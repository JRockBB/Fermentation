include(joinpath(@__DIR__, "..", "scripts", "fetch", "download.jl"))

@testset "Reproducible downloads" begin
    mktempdir() do directory
        original = joinpath(directory, "source.bin")
        bytes = UInt8[0x00, 0xff, 0x0a, 0x0d, 0x42]
        write(original, bytes)
        url = "file://" * (Sys.iswindows() ? "/" : "") * replace(abspath(original), '\\' => '/', " " => "%20")
        dest = joinpath(directory, "raw", "copy.bin")
        digest = FetchData.file_sha256(original)
        args = (; expected_sha256=digest, source="synthetic fixture", license="test")
        result = FetchData.fetch_file(url, dest; args...)
        @test !result.reused
        @test read(dest) == bytes
        metadata = TOML.parsefile(result.provenance)
        @test metadata["url"] == url
        @test metadata["sha256"] == digest
        @test metadata["bytes"] == length(bytes)
        @test endswith(metadata["retrieved_at_utc"], "Z")
        original_metadata = read(result.provenance)
        @test FetchData.fetch_file(url, dest; args...).reused
        @test read(result.provenance) == original_metadata
        @test_throws ErrorException FetchData.fetch_file(url * "?changed", dest; args...)
        @test_throws ErrorException FetchData.fetch_file(url, dest; source="different", license="test")
        @test_throws ErrorException FetchData.fetch_file(url, dest; args..., expected_sha256="0"^64)
        failed = joinpath(directory, "bad.bin")
        @test_throws ErrorException FetchData.fetch_file(url, failed; expected_sha256="0"^64)
        @test !ispath(failed)
        @test !ispath(failed * ".provenance.toml")
        @test_throws ArgumentError FetchData.fetch_file(url, failed; expected_sha256="bad")
        # Corruption is detected and the existing bytes are left untouched.
        write(dest, "changed fixture")
        @test_throws ErrorException FetchData.fetch_file(url, dest; args...)
        @test read(dest, String) == "changed fixture"
        orphan = joinpath(directory, "orphan.bin")
        write(orphan, bytes)
        @test_throws ErrorException FetchData.fetch_file(url, orphan)
    end
end
