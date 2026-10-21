module CocoaBayesAudit

using FermentModels, Statistics, MCMCDiagnosticTools, TOML
using Distributions: Normal, logpdf
include("diagnose_cocoa_bayes.jl")
const D = CocoaBayesDiagnostics
const C = D.C
const B = D.B

mutable struct Checks
    count::Int
end
function require!(checks,condition,message)
    condition === true || throw(ArgumentError(message))
    checks.count += 1
    return nothing
end

# Diagnostic NaNs are legitimate for failed chains; they must agree on both sides.
matches(a::Real,b::Real) = (isnan(a) && isnan(b)) || isapprox(a,b;rtol=1e-10,atol=1e-12)
matches(a::AbstractArray,b::AbstractArray) = size(a)==size(b) && all(matches.(a,b))

"""Strict reader for the unquoted numeric/identifier tables emitted by this runner."""
function table(path,header)
    lines = readlines(path)
    !isempty(lines) && split(first(lines),',')==collect(header) ||
        throw(ArgumentError("header differs: $path"))
    rows = [split(line,',';keepempty=true) for line in lines[2:end]]
    all(length(row)==length(header) for row in rows) ||
        throw(ArgumentError("row width differs: $path"))
    return rows
end

function audit_summary!(checks,saved,path)
    samples,chains = saved.settings["samples_per_chain"],saved.settings["chains"]
    rows = table(joinpath(path,"parameters.csv"),
        ("parameter","mean","median","lower_95","upper_95","rhat","ess_bulk","ess_tail"))
    require!(checks,getindex.(rows,1)==String.(collect(D.NAMES)),"parameter summary order differs")
    rhats,bulks,tails = Float64[],Float64[],Float64[]
    for (name,row) in zip(D.NAMES,rows)
        draws = reshape([r.parameters[name] for r in saved.rows],samples,chains)
        q = quantile(vec(draws),[.025,.5,.975])
        r,b,t = only(rhat(draws;kind=:rank)),only(ess(draws;kind=:bulk)),only(ess(draws;kind=:tail))
        expected = [mean(draws),q[2],q[1],q[3],r,b,t]
        require!(checks,matches(parse.(Float64,row[2:end]),expected),"summary differs for $name")
        push!(rhats,r); push!(bulks,b); push!(tails,t)
    end
    passed = chains>=2 && saved.report["divergences"]==0 && saved.report["depth_hits"]==0 &&
        all(isfinite,rhats) && maximum(rhats)<1.01 &&
        all(isfinite,bulks) && minimum(bulks)>=100chains &&
        all(isfinite,tails) && minimum(tails)>=100chains
    require!(checks,saved.report["diagnostics_passed"]===passed,"diagnostic decision differs")
    require!(checks,saved.report["status"]==(passed ? "posterior_diagnostics_passed" :
        "posterior_diagnostics_failed"),"diagnostic status differs")
    require!(checks,saved.report["rhat_threshold"]==1.01 &&
        saved.report["minimum_ess"]==100chains,"diagnostic threshold differs")
    for (name,value) in (("max_rhat",maximum(rhats)),("min_ess_bulk",minimum(bulks)),
            ("min_ess_tail",minimum(tails)))
        require!(checks,matches(saved.report[name],value),"report $name differs")
    end
    return passed
end

function audit_likelihood!(checks,saved,trial,model,path,replay_per_chain)
    observations = table(joinpath(path,"likelihood_observations.csv"),
        ("column","flat_index","time_h","variable","observed_normalized"))
    indices = findall(!ismissing,trial.data.values)
    nt = length(trial.data.times)
    expected = [("obs_$i",i,trial.data.times[mod1(i,nt)],String(COCOA_STATES[cld(i,nt)]),
        trial.data.values[i]) for i in indices]
    require!(checks,length(observations)==length(expected),"observation map length differs")
    for (row,e) in zip(observations,expected)
        require!(checks,(row[1],parse(Int,row[2]),parse(Float64,row[3]),row[4],parse(Float64,row[5]))==e,
            "observation map differs")
    end
    samples,chains = saved.settings["samples_per_chain"],saved.settings["chains"]
    n = min(replay_per_chain,samples)
    grid = n==1 ? [cld(samples,2)] : round.(Int,range(1,samples;length=n))
    replays = 0
    open(joinpath(path,"log_likelihood.csv")) do io
        eof(io) && throw(ArgumentError("empty likelihood table"))
        require!(checks,split(readline(io),',')==["chain","draw",("obs_$i" for i in indices)...],
            "likelihood header differs")
        count = 0
        for line in eachline(io)
            count += 1
            require!(checks,count<=samples*chains,"extra likelihood draw at row $count")
            row = split(line,',';keepempty=true)
            require!(checks,length(row)==2+length(indices),"likelihood width differs at row $count")
            chain,draw = cld(count,samples),mod1(count,samples)
            require!(checks,parse(Int,row[1])==chain && parse(Int,row[2])==draw,
                "likelihood identity differs at row $count; expected chain $chain draw $draw")
            values = parse.(Float64,row[3:end])
            require!(checks,all(isfinite,values),"nonfinite likelihood at chain $chain draw $draw")
            if draw in grid
                parameters = saved.rows[count].parameters
                prediction = model.predict(parameters,trial.data)
                require!(checks,!isnothing(prediction) && length(prediction)==length(trial.data.values) &&
                    all(i->isfinite(prediction[i]),indices),"likelihood replay prediction failed")
                expected_values = [logpdf(Normal(prediction[i],parameters.sigma),trial.data.values[i]) for i in indices]
                require!(checks,matches(values,expected_values),"replayed likelihood differs at chain $chain draw $draw")
                replays += 1
            end
        end
        require!(checks,count==samples*chains==saved.report["likelihood_draws"],
            "likelihood draw count differs")
    end
    return replays
end

function audit_predictive!(checks,saved,trial,model,path)
    selection = table(joinpath(path,"predictive_selection.csv"),("column","chain","draw"))
    samples,chains = saved.settings["samples_per_chain"],saved.settings["chains"]
    require!(checks,length(selection)==saved.settings["predictive_draws"]>0,
        "predictive selection count differs")
    selected = Int[]
    for (i,row) in enumerate(selection)
        column,chain,draw = parse.(Int,row)
        require!(checks,column==i && 1<=chain<=chains && 1<=draw<=samples,
            "predictive selection index differs")
        push!(selected,(chain-1)*samples+draw)
    end
    require!(checks,length(unique(selected))==length(selected),"duplicate predictive selection")
    columns = map(selected) do i
        prediction = model.predict(saved.rows[i].parameters,trial.data)
        require!(checks,!isnothing(prediction) && length(prediction)==length(trial.data.values) &&
            all(isfinite,prediction),"predictive replay failed")
        prediction
    end
    latent = hcat(columns...)
    rows = table(joinpath(path,"predictive.csv"),
        ("time_h","variable","observed_normalized","latent_mean_normalized","latent_median_normalized",
         "latent_lower_normalized","latent_upper_normalized","observation_lower_normalized",
         "observation_upper_normalized","source_scale_mg_g","observed_mg_g","latent_mean_mg_g",
         "latent_median_mg_g","latent_lower_mg_g","latent_upper_mg_g","observation_lower_mg_g",
         "observation_upper_mg_g"))
    require!(checks,length(rows)==length(trial.data.values),"predictive row count differs")
    nt = length(trial.data.times)
    for (i,row) in enumerate(rows)
        state,time = cld(i,nt),mod1(i,nt)
        require!(checks,parse(Float64,row[1])==trial.data.times[time] &&
            row[2]==String(COCOA_STATES[state]),"predictive state/time order differs")
        q = quantile(latent[i,:],[.025,.5,.975])
        require!(checks,matches(parse.(Float64,row[4:7]),[mean(latent[i,:]),q[2],q[1],q[3]]),
            "latent predictive bands differ at state $(COCOA_STATES[state]), time $(trial.data.times[time]) h")
        scale = trial.data.scales[state]
        if ismissing(scale)
            require!(checks,isempty(row[3]) && all(isempty,row[8:end]),
                "unobserved state has measurement or physical fields")
        else
            observed = trial.data.values[i]
            require!(checks,ismissing(observed) ? isempty(row[3]) : parse(Float64,row[3])==observed,
                "normalized observation differs")
            physical = trial.observations[time,state]
            require!(checks,ismissing(physical) ? isempty(row[11]) : parse(Float64,row[11])==physical,
                "physical observation differs")
            require!(checks,parse(Float64,row[10])==scale,"source scale differs")
            normalized,physical_bands = parse.(Float64,row[4:9]),parse.(Float64,row[12:17])
            require!(checks,all(isfinite,normalized) && all(isfinite,physical_bands) &&
                normalized[5]<=normalized[6] && matches(physical_bands,normalized.*scale),
                "measurement bands or physical conversion differ")
        end
    end
    return length(selected)
end

function audit_comparison!(checks,directory,saved,trial)
    header = ("trial","in_paper_benchmark","observations","status","max_rhat","min_ess_bulk","min_ess_tail",
        "divergences","depth_hits","retained_draws","elapsed_seconds","published_baseline_success",
        "published_max_rhat","published_min_bulk_ess","published_min_tail_ess","published_bma_weight",
        "psis_loo","published_psis_loo")
    rows = table(joinpath(directory,"comparison.csv"),header)
    require!(checks,getindex.(rows,1)==saved.settings["trials"],"comparison trial order differs")
    row = only(filter(r->r[1]==string(trial.trial),rows))
    columns = Dict(zip(header,row))
    require!(checks,parse(Bool,columns["in_paper_benchmark"])==trial.in_paper_benchmark &&
        columns["status"]==saved.report["status"],"comparison status or membership differs")
    for key in ("observations","max_rhat","min_ess_bulk","min_ess_tail","divergences","depth_hits",
            "retained_draws","elapsed_seconds")
        require!(checks,matches(parse(Float64,columns[key]),saved.report[key]),"comparison $key differs")
    end
    references = filter(r->r.trial==trial.trial,B.published_baseline())
    reference = isempty(references) ? nothing : only(references)
    require!(checks,parse(Bool,columns["published_baseline_success"])==
        (!isnothing(reference) && reference.success),"published success differs")
    for (key,field) in (("published_max_rhat",:max_rhat),("published_min_bulk_ess",:min_bulk_ess),
            ("published_min_tail_ess",:min_tail_ess),("published_bma_weight",:bma_weight))
        value = isnothing(reference) ? missing : getproperty(reference,field)
        require!(checks,ismissing(value) ? isempty(columns[key]) : matches(parse(Float64,columns[key]),value),
            "published $key differs")
    end
    require!(checks,isempty(columns["psis_loo"]) && isempty(columns["published_psis_loo"]) &&
        saved.report["direct_metric_comparison_available"]===false,"unexpected direct PSIS comparison")
end

"""
Audit completed Bayesian benchmark artifacts without changing their acceptance.
All likelihood rows are checked for structure and finiteness; a deterministic
grid per chain is replayed. Every selected latent predictive draw is replayed.
Observation-noise draws were not saved, so their bands receive structural and
unit checks only. An integrity pass does not imply passing MCMC diagnostics.
"""
function audit_bayes(directory;output,trials=nothing,replay_per_chain=3)
    ispath(output) && throw(ArgumentError("use a new audit output file"))
    isfile(joinpath(directory,"aborted_run.toml")) &&
        throw(ArgumentError("an interrupted experiment cannot receive a completed-fit audit"))
    replay_per_chain isa Integer && !(replay_per_chain isa Bool) && replay_per_chain>0 ||
        throw(ArgumentError("replay_per_chain must be a positive integer"))
    settings = TOML.parsefile(joinpath(directory,"settings.toml"))
    available = filter(t->string(t.trial) in settings["trials"],B.benchmark_trials(;scope=Symbol(settings["scope"])))
    string.([t.trial for t in available])==settings["trials"] ||
        throw(ArgumentError("saved trial list differs from scope order"))
    selected = C.F.select_trials(available,trials)
    isempty(selected) && throw(ArgumentError("no selected trials"))
    checks = Checks(0)
    require!(checks,get(settings,"time_unit",nothing)=="h" &&
        get(settings,"source_concentration_unit",nothing)=="mg/g pulp","source unit labels differ")
    files = ("scripts/fit_cocoa_benchmark_bayes.jl","scripts/fit_cocoa_benchmark.jl",
        "scripts/cocoa_fit_output.jl","src/core/bayes.jl",C.F.LIKELIHOOD_IMPLEMENTATION...)
    for file in files
        require!(checks,get(settings["implementation_sha256"],file,nothing)==B.digest(joinpath(B.ROOT,file)),
            "saved implementation differs: $file")
    end
    parent_hashes = Dict("settings.toml"=>B.digest(joinpath(directory,"settings.toml")),
        "comparison.csv"=>B.digest(joinpath(directory,"comparison.csv")))
    reports = Dict{String,Any}()
    for trial in selected
        id = string(trial.trial)
        saved = D.saved_draws(directory,trial)
        path = joinpath(directory,id)
        for (file,hash) in saved.provenance["parent_files_sha256"]
            parent_hashes["$id/$file"] = hash
        end
        for file in ("parameters.csv","log_likelihood.csv","likelihood_observations.csv",
                "predictive.csv","predictive_selection.csv")
            parent_hashes["$id/$file"] = B.digest(joinpath(path,file))
        end
        before = checks.count
        passed = audit_summary!(checks,saved,path)
        initial = NamedTuple{D.NAMES}(ntuple(_->.5,length(D.NAMES)))
        lower = NamedTuple{D.NAMES}(ntuple(_->1e-5,length(D.NAMES)))
        upper = NamedTuple{D.NAMES}(ntuple(_->5.,length(D.NAMES)))
        model = cocoa_fit_model(;initial,lower,upper,solver=Symbol(settings["solver"]),
            abstol=settings["abstol"],reltol=settings["reltol"],maxiters=saved.input["solver_maxiters"])
        replayed = audit_likelihood!(checks,saved,trial,model,path,replay_per_chain)
        predictive = audit_predictive!(checks,saved,trial,model,path)
        audit_comparison!(checks,directory,saved,trial)
        reports[id] = Dict("integrity_passed"=>true,"checks"=>checks.count-before,
            "sampling_diagnostics_passed"=>passed,"status"=>saved.report["status"],
            "likelihood_rows_checked"=>length(saved.rows),"likelihood_draws_replayed"=>replayed,
            "likelihood_replay_complete"=>replayed==length(saved.rows),
            "likelihood_values_checked"=>length(saved.rows)*saved.report["observations"],
            "latent_predictive_draws_replayed"=>predictive)
    end
    # Detect modification during the audit before publishing a successful result.
    for (file,hash) in parent_hashes
        require!(checks,B.digest(joinpath(directory,file))==hash,"parent changed during audit: $file")
    end
    result = Dict("integrity_passed"=>true,"checks"=>checks.count,"trials"=>reports,
        "julia_version"=>string(VERSION),"replay_per_chain"=>replay_per_chain,
        "solver_budget_scope"=>"One complete Julia ODE solve; saveat does not reset maxiters. This differs from Stan's between-output-time step budget.",
        "stan_budget_reference"=>"https://mc-stan.org/docs/2_24/stan-users-guide/control-parameters-for-ode-solving.html",
        "parent_directory"=>replace(relpath(abspath(directory),B.ROOT),'\\'=>'/'),
        "parent_sha256"=>parent_hashes,
        "limitations"=>"Integrity is separate from MCMC acceptance and biological validity. Likelihood replay covers the recorded count per trial; other rows receive structural checks only. Latent bands are replayed completely. Observation-noise draws are unavailable; measurement bands receive only structural and unit checks. No PSIS-LOO score is computed.",
        "implementation_sha256"=>Dict(file=>B.digest(joinpath(B.ROOT,file)) for file in
            ("scripts/audit_cocoa_benchmark_bayes.jl","scripts/diagnose_cocoa_bayes.jl",files...)))
    C.F.write_toml(output,result)
    return result
end

end

if abspath(PROGRAM_FILE)==@__FILE__
    2<=length(ARGS)<=3 || error("Usage: audit_cocoa_benchmark_bayes.jl FIT_DIRECTORY NEW_AUDIT.toml [TRIAL,TRIAL,...]")
    result = CocoaBayesAudit.audit_bayes(ARGS[1];output=ARGS[2],
        trials=length(ARGS)==3 ? Symbol.(split(ARGS[3],',')) : nothing)
    println("Artifact integrity passed: ",result["checks"]," checks. See ",ARGS[2]," for separate sampling status.")
end
