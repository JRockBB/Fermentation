module CocoaBenchmarkData

using SHA
using TOML

const ROOT = normpath(joinpath(@__DIR__, ".."))
const RAW = joinpath("data", "raw", "cocoa_benchmark")
const ARCHIVE_DIR = "mmorenozam-ecbfm-5ced89b"
const BENCHMARK_STATES = (:G, :F, :E, :L, :A, :Xy, :Xl, :Xa)
# Main paper Table 1, in printed order. The other nine trials remain available.
const BENCHMARK_TRIALS = (:ghhp1, :dowb1, :ghhp2, :ghhp3, :ghhp4, :ghhp5,
    :brwb1, :brwb2, :ecpt1, :ecpt2, :ecwb1, :ecwb2, :brpb1, :brst1,
    :brwb3, :brwb4, :brst2, :brwb7, :mywb3, :hnwb1, :ecpb1, :niwb1, :niwb2)
const HIDDEN_LACTATE_TRIALS = (:ghhp2, :ghhp3, :ghhp4, :ghhp5)
const MECHANISMS = (
    (id=1, added=3, name="First-order loss of ethanol; lactate; acetate", location="Supplement S3.1 equations 5-8"),
    (id=2, added=7, name="LAB fructose uptake and additional metabolite yields", location="Supplement S3.2 equations 9-14"),
    (id=3, added=2, name="Yeast acetate production from glucose and fructose", location="Supplement S3.3 equation 15"),
    (id=4, added=4, name="Yeast lactate uptake and ethanol production", location="Supplement S3.4 equations 16-19"),
    (id=5, added=3, name="AAB acetate uptake and growth", location="Supplement S3.5 equations 20-22"),
)

digest(path) = bytes2hex(open(sha256, path))

"""Verify downloaded originals and every byte-preserved archive member."""
function verify_benchmark(; root=ROOT)
    raw = joinpath(root, RAW)
    for (file, hash, source, license) in (
        ("ecbfm-v1.0.zip", "735255d908338acd573d993d6ab99963a567a3a81773bd910156c2a9754364c4",
         "doi:10.5281/zenodo.5510370", "other-open (Zenodo; no license file in archive)"),
        ("rsos210274_si_001.pdf", "54fbaa2b5942f2b4464fade8a3f9987476da1038c06eac10a6c35458048b6f92",
         "doi:10.6084/m9.figshare.19149492.v1", "CC-BY-4.0"))
        path = joinpath(raw, file)
        digest(path) == hash || error("Benchmark checksum differs: $file")
        p = TOML.parsefile(path * ".provenance.toml")
        p["sha256"] == hash && p["bytes"] == filesize(path) &&
            p["source"] == source && p["license"] == license || error("Benchmark provenance differs: $file")
    end
    expected = Set{String}()
    for line in readlines(joinpath(root, "docs", "cocoa_benchmark_archive_files.csv"))[2:end]
        file, bytes, hash = split(line, ',')
        path = normpath(joinpath(raw, file))
        startswith(path, normpath(joinpath(raw, ARCHIVE_DIR)) * string(Base.Filesystem.path_separator)) ||
            error("Archive manifest path escapes extraction directory")
        filesize(path) == parse(Int, bytes) && digest(path) == hash || error("Archive member differs: $file")
        push!(expected, path)
    end
    actual = Set(normpath(joinpath(dir, file))
        for (dir, _, files) in walkdir(joinpath(raw, ARCHIVE_DIR)) for file in files)
    actual == expected && length(expected) == 78 || error("Archive member set differs")
    return length(expected)
end

# Strict literal parsing: never evaluate downloaded R code or expressions.
function literals(text)
    tokens = strip.(split(text, ','))
    all(t -> occursin(r"^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$", t), tokens) ||
        error("Expected finite numeric R literals only")
    values = parse.(Float64, tokens)
    all(isfinite, values) || error("Nonfinite source literal")
    return values
end

function one_match(pattern, text)
    matches = collect(eachmatch(pattern, text))
    length(matches) == 1 || error("Expected exactly one source field: $pattern")
    return only(matches)
end

"""Parse only the 28 eight-series or four seven-series dataset blocks."""
function parse_driver(text; hidden_lactate=false)
    source = replace(text, r"#[^\n]*" => "")
    blocks = collect(eachmatch(r"(?s)if\s*\(args\[1\]\s*==\s*(\d+)\)\s*\{(.*?)\}", source))
    expected = hidden_lactate ? 4 : 28
    length(blocks) == expected || error("Unexpected number of source datasets")
    result = NamedTuple[]
    for (index, block) in enumerate(blocks)
        parse(Int, block[1]) == index || error("Dataset branch order changed")
        text = block[2]
        trial = Symbol(one_match(r"lbl\s*<-\s*'([a-z0-9]+)_'", text)[1])
        rows = parse(Int, one_match(r"\bT\s*=\s*(\d+)\b", text)[1])
        t0 = only(literals(one_match(r"\bt0\s*=\s*([^\r\n]+)", text)[1]))
        timepoints = literals(one_match(r"\bts\s*=\s*c\(([^)]*)\)", text)[1])
        initial_pattern = hidden_lactate ? r"\bxi_r\s*=\s*c\(([^)]*)\)" : r"\bx0\s*=\s*c\(([^)]*)\)"
        initial_source = literals(one_match(initial_pattern, text)[1])
        m = one_match(r"(?s)\bx\s*=\s*structure\(c\((.*?)\),\s*\.Dim\s*=\s*c\((\d+),\s*(\d+)\)\)", text)
        n, k = parse.(Int, (m[2], m[3]))
        n == rows == length(timepoints) && k == length(initial_source) == (hidden_lactate ? 7 : 8) ||
            error("Source dimensions differ")
        numbers = literals(m[1])
        length(numbers) == n*k && all(>=(0), [initial_source; numbers]) ||
            error("Source concentrations differ")
        all(>(0), diff([t0; timepoints])) || error("Source times must increase")
        matrix = reshape(numbers, n, k)
        # Source Camu-2008 order G,F,E,A,Xy,Xl,Xa; hidden L is last in Stan.
        columns = hidden_lactate ? [1,2,3,5,6,7,8] : collect(1:8)
        physical_initial = Vector{Union{Missing,Float64}}(missing, 8)
        observations = Matrix{Union{Missing,Float64}}(missing, n, 8)
        scales = Vector{Union{Missing,Float64}}(missing, 8)
        physical_initial[columns] = initial_source
        observations[:,columns] = matrix
        maxima = vec(maximum(vcat(initial_source', matrix); dims=1))
        all(>(0), maxima) || error("Each observed series must have a positive maximum")
        scales[columns] = maxima
        initial = zeros(8)
        initial[columns] = initial_source ./ maxima
        # post_sims2.R:139 appends a fixed zero for the unobserved normalized LA.
        # It supplies neither a measured LA value nor a physical LA scale.
        normalized = Matrix{Union{Missing,Float64}}(missing, n, 8)
        normalized[:,columns] = matrix ./ maxima'
        data = (; t0, times=timepoints, initial, scales, values=vec(normalized), sigma=:sigma)
        push!(result, (; trial, branch=index, hidden_lactate, physical_initial,
            observations, data, in_paper_benchmark=trial in BENCHMARK_TRIALS))
    end
    length(unique(r.trial for r in result)) == expected || error("Duplicate trial identifiers")
    hidden_lactate && Set(r.trial for r in result) != Set(HIDDEN_LACTATE_TRIALS) &&
        error("Unexpected hidden-lactate trial")
    return result
end

function benchmark_trials(; root=ROOT, scope=:benchmark, verify=true)
    scope in (:benchmark, :all) || throw(ArgumentError("scope must be :benchmark or :all"))
    verify && verify_benchmark(; root)
    rows = NamedTuple[]
    for (file, hidden_lactate) in (("post_sims1.R", false), ("post_sims2.R", true))
        source_file = replace(joinpath(RAW, ARCHIVE_DIR, "posteriors", file), '\\'=>'/')
        path = joinpath(root, source_file)
        append!(rows, [merge(row, (; source_file, source_sha256=digest(path)))
            for row in parse_driver(read(path, String); hidden_lactate)])
    end
    length(unique(r.trial for r in rows)) == 32 || error("Expected 32 unique trials")
    Set(r.trial for r in rows if r.in_paper_benchmark) == Set(BENCHMARK_TRIALS) ||
        error("Published benchmark trial set differs")
    scope == :all && return rows
    return [only(filter(r -> r.trial == id, rows)) for id in BENCHMARK_TRIALS]
end

"""List all 32 combinations and verify their parameter counts in both Stan families."""
function mechanism_variants(; root=ROOT)
    rows = NamedTuple[]
    for mask in 0:31
        ids = [i for i in 1:5 if !iszero(mask & (1 << (i-1)))]
        code = isempty(ids) ? "0" : join(ids)
        kinetics = 24 + sum((MECHANISMS[i].added for i in ids); init=0)
        for (folder, prefix) in (("", "mm"), ("camu2008", "mmc"))
            path = joinpath(root, RAW, ARCHIVE_DIR, "Stan_models", folder, "$prefix$code.stan")
            text = replace(read(path, String), r"//[^\n]*" => "")
            block = one_match(r"(?s)\nparameters\s*\{(.*?)\}", text)[1]
            count = length(collect(eachmatch(r"real\s*<lower=0>\s*\w+\s*;", block)))
            count == kinetics + 1 || error("Parameter count differs for $prefix$code")
        end
        push!(rows, (; code, mechanisms=join(ids, ';'), kinetic_parameters=kinetics,
            fitted_parameters=kinetics+1, complete_source="Stan_models/mm$code.stan",
            hidden_source="Stan_models/camu2008/mmc$code.stan"))
    end
    return rows
end

"""Printed baseline diagnostics, with independently recorded transcription corrections."""
function published_baseline(; root=ROOT)
    raw = joinpath(root, "data", "raw", "transcribed", "cocoa2022")
    for (file, hash) in (
        ("baseline_diagnostics.csv", "62537f45a3332b2e0eae59f98c2270474d70d337f5de90b3158db8797dd38479"),
        ("baseline_diagnostics_corrections.csv", "8550c9c43c800e14a7f28ba4096dc35deca2fec17c65ce4c0d421658fb317bb3"))
        digest(joinpath(raw,file)) == hash || error("Source transcription differs: $file")
    end
    lines = readlines(joinpath(raw, "baseline_diagnostics.csv"))
    header = split(first(lines), ',')
    header == ["trial","published_baseline_success","max_rhat","min_bulk_ess",
        "min_tail_ess","mean_ess_fraction","bma_weight"] || error("Unexpected diagnostics header")
    entries = [Dict(zip(header, split(line, ','; keepempty=true))) for line in lines[2:end]]
    for line in readlines(joinpath(raw, "baseline_diagnostics_corrections.csv"))[2:end]
        id, field, old, new, _, _ = split(line, ',')
        row = only(filter(r -> r["trial"] == id, entries))
        row[field] == old || error("Transcription correction precondition failed")
        row[field] = new
    end
    [Symbol(r["trial"]) for r in entries] == collect(BENCHMARK_TRIALS) ||
        error("Published diagnostic trial order differs")
    rows = [merge((; trial=Symbol(r["trial"]), success=parse(Bool,r["published_baseline_success"])),
        NamedTuple{Tuple(Symbol.(header[3:end]))}(Tuple(isempty(r[h]) ? missing : parse(Float64,r[h])
            for h in header[3:end]))) for r in entries]
    all(r -> r.success ? all(x -> !ismissing(x) && isfinite(x), Tuple(r)[3:end]) :
        all(ismissing, Tuple(r)[3:end]), rows) || error("Inconsistent published success and diagnostics")
    return rows
end

csv_value(x) = ismissing(x) ? "" : string(x)
function write_csv(path, header, rows)
    open(path, "w") do io
        println(io, join(header, ','))
        for row in rows
            fields = csv_value.(row)
            any(s -> occursin(r"[,\r\n\"]", s), fields) && error("Unexpected CSV delimiter")
            println(io, join(fields, ','))
        end
    end
end

function prepare_benchmark(; root=ROOT,
        output=joinpath(root, "data", "processed", "cocoa_benchmark"))
    trials = benchmark_trials(; root, scope=:all)
    variants = mechanism_variants(; root)
    references = published_baseline(; root)
    mkpath(output)
    for trial in trials
        rows = Tuple[]
        for (j, state) in enumerate(BENCHMARK_STATES)
            push!(rows, (trial.trial, trial.data.t0, state, trial.physical_initial[j],
                ismissing(trial.physical_initial[j]) ? "unobserved_initial" : "initial",
                "mg/g pulp", trial.data.scales[j], trial.data.initial[j]))
            for (i, t) in enumerate(trial.data.times)
                value = trial.observations[i,j]
                push!(rows, (trial.trial, t, state, value, ismissing(value) ? "unobserved" : "observation",
                    "mg/g pulp", trial.data.scales[j], trial.data.values[(j-1)*length(trial.data.times)+i]))
            end
        end
        write_csv(joinpath(output, "$(trial.trial).csv"),
            ("trial","time_h","variable","value","role","unit","scale","normalized_value"), rows)
    end
    write_csv(joinpath(output,"datasets.csv"),
        ("trial","in_paper_benchmark","observed_states","fitted_observations","timepoints",
         "last_time_h","source_file","source_branch","source_sha256"),
        [(t.trial,t.in_paper_benchmark,t.hidden_lactate ? 7 : 8,count(!ismissing,t.data.values),
          length(t.data.times),last(t.data.times),t.source_file,t.branch,t.source_sha256) for t in trials])
    write_csv(joinpath(output,"mechanism_variants.csv"), keys(first(variants)), Tuple.(variants))
    write_csv(joinpath(output,"mechanisms.csv"), keys(first(MECHANISMS)), Tuple.(MECHANISMS))
    write_csv(joinpath(output,"published_baseline.csv"), keys(first(references)), Tuple.(references))
    paths = ["scripts/prepare_cocoa_benchmark.jl", "docs/cocoa_benchmark_archive_files.csv",
        joinpath(RAW,"ecbfm-v1.0.zip"),joinpath(RAW,"rsos210274_si_001.pdf"),
        "data/raw/transcribed/cocoa2022/baseline_diagnostics.csv",
        "data/raw/transcribed/cocoa2022/baseline_diagnostics_corrections.csv"]
    provenance = Dict("source"=>"doi:10.1098/rsos.210274; doi:10.5281/zenodo.5510370",
        "time_units"=>"h", "observed_state_units"=>"mg/g pulp", "all_trials"=>32,
        "paper_benchmark_trials"=>23, "state_order"=>String.(collect(BENCHMARK_STATES)),
        "transformations"=>"Parse numeric literals; reorder Camu-2008 states; preserve missing LA; normalize observed series by maxima including initial row",
        "hidden_lactate"=>"Fixed normalized initial zero from post_sims2.R:139; physical initial and scale unknown; all LA observations missing",
        "archive_license"=>"Zenodo other-open; no explicit archive license file; not assumed CC0",
        "inputs_sha256"=>Dict(replace(p,'\\'=>'/')=>digest(joinpath(root,p)) for p in paths))
    open(joinpath(output,"provenance.toml"),"w") do io
        TOML.print(io, provenance; sorted=true)
    end
    println("Prepared 32 trials (23 in published benchmark), 32 variants, and baseline references.")
    return (; trials, variants, references)
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    CocoaBenchmarkData.prepare_benchmark()
end
