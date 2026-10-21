function with_parameter_file(f, entries)
    mktempdir() do directory
        path = joinpath(directory, "parameters with spaces.toml")
        open(path, "w") do io
            TOML.print(io, entries)
        end
        return f(path)
    end
end

function synthetic_entry()
    return Dict{String, Any}(
        "value" => 0.25,
        "units" => "1/h",
        "source" => "synthetic:test/parameters.jl",
        "location" => "Synthetic parser input",
        "confidence" => "high",
    )
end

@testset "Parameter and provenance loader" begin
    @testset "Values and provenance stay aligned" begin
        path = joinpath(@__DIR__, "fixtures", "params", "valid.toml")
        result = load_params(path)
        @test keys(result) == (:values, :provenance)
        values, provenance = result
        @test values === (a_rate=0.25, z_integer=2)
        @test [row.name for row in provenance] == collect(keys(values))
        @test provenance[1] == (
            name=:a_rate,
            units="1/h",
            source="synthetic:test/fixtures/params/valid.toml",
            location="Floating-point parser fixture",
            notes="Synthetic input; not suitable for a fermentation model.",
            confidence="high",
        )
        @test provenance[2].notes == ""
        @test provenance[2].units == "1"
        @test provenance[2].confidence == "medium"

        for confidence in ("high", "medium", "low")
            entry = merge(synthetic_entry(), Dict("confidence" => confidence))
            with_parameter_file(Dict("μmax" => entry)) do path
                loaded = cd(dirname(path)) do
                    load_params(basename(path))
                end
                @test loaded.values.μmax === 0.25
                @test only(loaded.provenance).name == :μmax
                @test only(loaded.provenance).confidence == confidence
            end
        end

        with_parameter_file(Dict{String, Any}()) do path
            loaded = load_params(path)
            @test loaded.values === NamedTuple()
            @test isempty(loaded.provenance)
        end
    end

    @testset "Explicit placeholders" begin
        placeholder = Dict{String, Any}(
            "units" => "1/h",
            "confidence" => "placeholder",
            "notes" => "TODO: obtain a measured value and its citation.",
        )
        # Both an omitted value and TOML's explicit string sentinel mean missing.
        for entry in (placeholder, merge(placeholder, Dict("value" => "missing")))
            with_parameter_file(Dict("rate" => entry)) do path
                loaded = load_params(path)
                @test ismissing(loaded.values.rate)
                row = only(loaded.provenance)
                @test ismissing(row.source)
                @test ismissing(row.location)
                @test row.confidence == "placeholder"
            end
        end

        sourced = merge(placeholder, Dict("source" => "synthetic:test/parameters.jl",
                                         "location" => "Value awaiting extraction"))
        with_parameter_file(Dict("rate" => sourced)) do path
            row = only(load_params(path).provenance)
            @test row.source == sourced["source"]
            @test row.location == sourced["location"]
        end

        # Source-free values must retain an explicit placeholder confidence tag.
        numeric = merge(placeholder, Dict("value" => 0.5, "notes" => "Synthetic only."))
        with_parameter_file(Dict("rate" => numeric)) do path
            loaded = load_params(path)
            @test loaded.values.rate === 0.5
            @test ismissing(only(loaded.provenance).source)
            @test only(loaded.provenance).confidence == "placeholder"
        end

        with_parameter_file(Dict("known" => synthetic_entry(), "unknown" => placeholder)) do path
            loaded = load_params(path)
            @test loaded.values.known === 0.25
            @test ismissing(loaded.values.unknown)
            @test [row.name for row in loaded.provenance] == [:known, :unknown]
        end

        for change in (Dict("units" => ""), Dict("value" => Inf), Dict("source" => 42))
            with_parameter_file(Dict("rate" => merge(placeholder, change))) do path
                @test_throws ArgumentError load_params(path)
            end
        end

        for notes in ("", "Obtain a measured value.", "NOTODO")
            entry = merge(placeholder, Dict("notes" => notes))
            with_parameter_file(Dict("rate" => entry)) do path
                @test_throws ArgumentError load_params(path)
            end
        end
        without_notes = copy(placeholder)
        delete!(without_notes, "notes")
        with_parameter_file(Dict("rate" => without_notes)) do path
            @test_throws ArgumentError load_params(path)
        end

        template = load_params(joinpath(@__DIR__, "..", "params", "template.toml"))
        @test ismissing(template.values.lab_mu_max)
        @test only(template.provenance).confidence == "placeholder"
    end

    @testset "Required provenance and units" begin
        for field in ("units", "source", "location", "confidence")
            @testset "$field" begin
                entry = synthetic_entry()
                delete!(entry, field)
                with_parameter_file(Dict("rate" => entry)) do path
                    @test_throws ArgumentError load_params(path)
                end
                for bad_text in ("", " \t ", 42, ["1/h"])
                    entry = merge(synthetic_entry(), Dict(field => bad_text))
                    with_parameter_file(Dict("rate" => entry)) do path
                        @test_throws ArgumentError load_params(path)
                    end
                end
            end
        end
        for confidence in ("medium", "low")
            entry = merge(synthetic_entry(), Dict("confidence" => confidence))
            delete!(entry, "source")
            with_parameter_file(Dict("rate" => entry)) do path
                @test_throws ArgumentError load_params(path)
            end
        end
        for change in (Dict("confidence" => "certain"), Dict("notes" => 42))
            with_parameter_file(Dict("rate" => merge(synthetic_entry(), change))) do path
                @test_throws ArgumentError load_params(path)
            end
        end

        # Errors identify the offending file, parameter, and field.
        entry = synthetic_entry()
        delete!(entry, "source")
        with_parameter_file(Dict("rate" => entry)) do path
            error = try
                load_params(path)
                nothing
            catch exception
                exception
            end
            @test error isa ArgumentError
            message = sprint(showerror, error)
            @test occursin(repr(path), message)
            @test occursin("rate", message)
            @test occursin("source", message)
        end
    end

    @testset "Scalar values and schema validation" begin
        for bad_value in (true, false, Inf, -Inf, NaN, "0.25", "missing", [1, 2],
                          Dict("nested" => 1), TOML.parse("value = 2000-01-01")["value"])
            with_parameter_file(Dict("rate" => merge(synthetic_entry(), Dict("value" => bad_value)))) do path
                @test_throws ArgumentError load_params(path)
            end
        end
        entry = synthetic_entry()
        delete!(entry, "value")
        with_parameter_file(Dict("rate" => entry)) do path
            @test_throws ArgumentError load_params(path)
        end

        # Bounds are model-specific: finite negative and zero values are retained.
        for value in (0, -2.5)
            with_parameter_file(Dict("rate" => merge(synthetic_entry(), Dict("value" => value)))) do path
                @test load_params(path).values.rate === value
            end
        end

        for entries in (Dict("rate" => 1), Dict("" => synthetic_entry()),
                        Dict("rate" => merge(synthetic_entry(), Dict("unit" => "1/h"))))
            with_parameter_file(entries) do path
                @test_throws ArgumentError load_params(path)
            end
        end
    end

    @testset "File and TOML errors" begin
        mktempdir() do directory
            @test_throws ErrorException load_params(joinpath(directory, "absent.toml"))
            path = joinpath(directory, "invalid.toml")
            for text in ("[unclosed", "[rate]\nvalue = 1\nvalue = 2\n", "[rate]\n[rate]\n")
                write(path, text)
                @test_throws TOML.ParserError load_params(path)
            end
        end
    end
end
