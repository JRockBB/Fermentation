# P1.7: a deliberately confounded synthetic observation model.
using FermentModels
using Random
using TOML

function identifiability_demo(output=joinpath(@__DIR__,"output","identifiability"))
    model = FitModel((p,d)->fill(p.a*p.b,length(d.values));
        initial=(a=1.,b=1.),lower=(a=0.1,b=0.1),upper=(a=10.,b=10.))
    data = (values=fill(2.,12),sigma=0.1)
    profile = profile_likelihood(model,data,:a,range(0.5,4.;length=15);
        starts=3,rng=MersenneTwister(67))
    profile.reference.converged || error("No reference fit converged")
    sensitivity = local_sensitivities(model,data,profile.reference.parameters)
    mkpath(output)
    open(joinpath(output,"profile.csv"),"w") do io
        println(io,"fixed_a,fitted_b,negative_log_likelihood,deviance,converged")
        for point in profile.points
            b = isnothing(point.parameters) ? missing : point.parameters.b
            println(io,join((point.value,b,point.objective,point.deviance,point.converged),','))
        end
    end
    open(joinpath(output,"sensitivities.csv"),"w") do io
        println(io,"observation,d_prediction_d_log_a_over_sigma,d_prediction_d_log_b_over_sigma")
        for i in axes(sensitivity.jacobian,1)
            println(io,join((i,sensitivity.jacobian[i,1],sensitivity.jacobian[i,2]),','))
        end
    end
    report = Dict("synthetic"=>true,"script"=>"examples/identifiability.jl",
        "seed"=>67,"profile_status"=>string(profile.status),"sensitivity_status"=>string(sensitivity.status),
        "rank"=>sensitivity.rank,"parameters"=>2,"singular_values"=>sensitivity.singular_values,
        "interpretation"=>"Only the product a*b is constrained by these synthetic observations.")
    open(joinpath(output,"report.toml"),"w") do io
        TOML.print(io,report;sorted=true)
    end
    println("Local rank: ",sensitivity.rank," / 2; profile status: ",profile.status)
    println("Results: ",abspath(output))
    return (;profile,sensitivity)
end

if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) ? identifiability_demo() : identifiability_demo(only(ARGS))
end
