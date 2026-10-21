include(joinpath(@__DIR__,"..","scripts","cocoa_fit_output.jl"))

@testset "Cocoa checkpoint columns preserve parameter identities" begin
    # Loading TOML sorts names; the ODE and sampler use equation order.
    loaded = load_params(joinpath(@__DIR__,"..","params","cocoa_camu_scaled.toml")).values
    names = (COCOA_PARAMETERS...,:sigma)
    canonical = NamedTuple{names}(Tuple(loaded[name] for name in names))
    @test keys(loaded) != keys(canonical)
    for columns in (keys(loaded),keys(canonical)), parameters in (loaded,canonical)
        io = IOBuffer()
        write_cocoa_draw_header(io,columns)
        write_cocoa_draw(io,2,7,parameters,(divergent=false,tree_depth=9),columns)
        header,row = split.(split(chomp(String(take!(io))),'\n'),',')
        record = Dict(zip(header,row))
        @test length(header) == length(row) == 29
        @test record["chain"] == "2" && record["draw"] == "7"
        @test all(parse(Float64,record[String(name)]) == loaded[name] for name in names)
        @test record["divergent"] == "false" && record["tree_depth"] == "9"
    end
end

@testset "Cocoa MLE restart parameters remain named and bounded" begin
    loaded = load_params(joinpath(@__DIR__,"..","params","cocoa_camu_scaled.toml")).values
    mktempdir() do directory
        path = joinpath(directory,"checkpoint.toml")
        tentative = merge(loaded,(mu1=loaded.mu1*1.1,)) # synthetic saved-start fixture
        open(path,"w") do io
            TOML.print(io,Dict("converged"=>false,"parameters"=>Dict(string(k)=>v for (k,v) in pairs(tentative))))
        end
        restored = cocoa_fit_initial(loaded,path)
        @test restored == tentative
        @test keys(restored) == keys(loaded)
        lower = NamedTuple{keys(loaded)}(ntuple(_ -> 1e-5,25))
        upper = NamedTuple{keys(loaded)}(ntuple(_ -> 5.,25))
        @test cocoa_fit_model(;initial=restored,lower,upper).initial.mu1 == tentative.mu1
        # A start is computational input, never an accepted estimate or relaxed bound.
        outside = merge(restored,(mu1=5.,))
        @test_throws DomainError cocoa_fit_model(;initial=outside,lower,upper)
        for bad in (Dict(),Dict("parameters"=>Dict()),Dict("parameters"=>Dict("unknown"=>1.)))
            open(path,"w") do io
                TOML.print(io,bad)
            end
            @test_throws ArgumentError cocoa_fit_initial(loaded,path)
        end
    end
end
