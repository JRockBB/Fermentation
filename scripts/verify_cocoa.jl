# P2.1: verify immutable archive members and inspect numeric R literals.
# This script never evaluates R code and never writes raw or processed data.
using SHA
using TOML

function numeric_literals(text)
    tokens = strip.(split(text,','))
    pattern = r"^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$"
    all(token -> occursin(pattern,token),tokens) || error("Expected only numeric R literals")
    return parse.(Float64,tokens)
end

function inspect_dataset(path)
    source = replace(read(path,String),r"#[^\n]*"=>"")
    initial_match = match(r"(?m)^x0\s*=\s*c\(([^)]*)\)",source)
    times_match = match(r"(?m)^ts\s*=\s*c\(([^)]*)\)",source)
    matrix_match = match(r"(?s)x\s*=\s*structure\(c\((.*?)\),\s*\.Dim\s*=\s*c\((\d+),\s*(\d+)\)\)",source)
    any(isnothing,(initial_match,times_match,matrix_match)) && error("Unexpected data layout in $path")
    initial = numeric_literals(initial_match[1])
    times = numeric_literals(times_match[1])
    numbers = numeric_literals(matrix_match[1])
    rows, columns = parse(Int,matrix_match[2]), parse(Int,matrix_match[3])
    columns == length(initial) == 8 || error("Expected eight source states")
    rows == length(times) && length(numbers) == rows*columns || error("Source data dimensions differ")
    all(>(0),diff([0.;times])) && all(>=(0),[initial;numbers]) || error("Unexpected source time or concentration")
    # R and Julia both store matrices column-major.
    observations = reshape(numbers,rows,columns)
    scales = vec(maximum(vcat(initial',observations);dims=1))
    return (;initial,times,observations,scales)
end

function verify_cocoa(root=joinpath(@__DIR__,".."))
    raw = joinpath(root,"data","raw","cocoa")
    archive = joinpath(raw,"CocoaModel.zip")
    expected = "9cc6c82ae32c7360408d0bef072e80462df882a004d7aa8975c52f301568c40a"
    bytes2hex(open(sha256,archive)) == expected || error("Archive SHA-256 mismatch")
    provenance = TOML.parsefile(archive*".provenance.toml")
    provenance["sha256"] == expected && provenance["bytes"] == filesize(archive) || error("Archive provenance mismatch")
    entries = readlines(joinpath(root,"docs","cocoa_archive_files.csv"))[2:end]
    expected_paths = String[]
    for entry in entries
        relative, size, digest = split(entry,',')
        path = joinpath(raw,split(relative,'/')...)
        isfile(path) || error("Missing archive member: $relative")
        filesize(path) == parse(Int,size) || error("File length differs: $relative")
        bytes2hex(open(sha256,path)) == digest || error("File SHA-256 differs: $relative")
        push!(expected_paths,normpath(path))
    end
    actual_paths = [normpath(joinpath(dir,file)) for (dir,_,files) in walkdir(joinpath(raw,"CocoaModel-master")) for file in files]
    Set(actual_paths) == Set(expected_paths) || error("Unexpected files in the extracted archive")
    println("Verified archive and ",length(entries)," byte-identical members.")
    supplement = joinpath(raw,"rsos180964_si_001.pdf")
    supplement_hash = "10a337b696e50ca417e915ac90e22064e04a65533b89fb3300fc5b7542084a85"
    bytes2hex(open(sha256,supplement)) == supplement_hash || error("Supplement SHA-256 mismatch")
    supplement_metadata = TOML.parsefile(supplement*".provenance.toml")
    supplement_metadata["sha256"] == supplement_hash &&
        supplement_metadata["bytes"] == filesize(supplement) &&
        supplement_metadata["source"] == "doi:10.6084/m9.figshare.7152614.v2" &&
        supplement_metadata["license"] == "CC-BY-4.0" || error("Supplement provenance mismatch")
    println("Verified separately acquired publisher supplement (CC BY 4.0).")
    total = 0
    for (name,expected_rows) in (("camu.R",16),("p_box1.R",13),("p_box2.R",13))
        data = inspect_dataset(joinpath(raw,"CocoaModel-master","simulation",name))
        size(data.observations) == (expected_rows,8) || error("Dataset dimensions changed")
        total += length(data.observations)
        println(name,": ",expected_rows," fitted times plus t0; ",length(data.observations)," likelihood observations; last time ",last(data.times))
        println("  source scaling maxima [Glc,Fru,EtOH,LA,Ac,Y,LAB,AAB]: ",data.scales)
    end
    println("Total likelihood observations: ",total," (24 initial-state values are fixed separately).")
    return total
end

if abspath(PROGRAM_FILE) == @__FILE__
    verify_cocoa()
end
