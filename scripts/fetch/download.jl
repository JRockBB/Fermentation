module FetchData

using Dates
using Downloads
using SHA
using TOML

file_sha256(path) = open(bytes2hex ∘ sha256, path)

"""
    fetch_file(url, destination; expected_sha256=nothing, source="", license="")

Download original bytes and write a `.provenance.toml` sidecar containing the
URL, UTC retrieval time, byte count, SHA-256, source, and license. Repeated calls
verify and reuse an existing pair; they never replace original data or metadata.
An incomplete pair or changed contents raises an error. Downloads are staged in
the destination directory and checked before publication. Concurrent writers to
the same destination are not supported.
"""
function fetch_file(url::AbstractString, destination::AbstractString;
                    expected_sha256=nothing, source="", license="")
    isempty(strip(url)) && throw(ArgumentError("URL must not be blank"))
    expected = if isnothing(expected_sha256)
        nothing
    elseif expected_sha256 isa AbstractString && occursin(r"^[0-9a-fA-F]{64}$", expected_sha256)
        lowercase(expected_sha256)
    else
        throw(ArgumentError("expected_sha256 must contain 64 hexadecimal characters"))
    end
    path = abspath(destination)
    sidecar = path * ".provenance.toml"
    if ispath(path) || ispath(sidecar)
        isfile(path) && isfile(sidecar) ||
            error("Incomplete download/provenance pair at $path; refusing to overwrite")
        metadata = TOML.parsefile(sidecar)
        get(metadata, "url", nothing) == url || error("Recorded URL differs for $path")
        get(metadata, "source", nothing) == source || error("Recorded source differs for $path")
        get(metadata, "license", nothing) == license || error("Recorded license differs for $path")
        digest = file_sha256(path)
        get(metadata, "sha256", nothing) == digest || error("Recorded SHA-256 mismatch for $path")
        get(metadata, "bytes", nothing) == filesize(path) || error("Recorded size mismatch for $path")
        isnothing(expected) || digest == expected || error("Expected SHA-256 mismatch for $path")
        return (; path, provenance=sidecar, sha256=digest, reused=true)
    end

    mkpath(dirname(path))
    temporary, handle = mktemp(dirname(path))
    close(handle)
    temporary_metadata, handle = mktemp(dirname(path))
    close(handle)
    try
        Downloads.download(url, temporary)
        digest = file_sha256(temporary)
        isnothing(expected) || digest == expected || error("Downloaded SHA-256 mismatch for $url")
        metadata = Dict(
            "url" => String(url),
            "retrieved_at_utc" => Dates.format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SS.sss") * "Z",
            "sha256" => digest,
            "bytes" => filesize(temporary),
            "source" => String(source),
            "license" => String(license),
        )
        open(temporary_metadata, "w") do io
            TOML.print(io, metadata; sorted=true)
        end
        # force=false preserves any existing raw material.
        mv(temporary, path; force=false)
        mv(temporary_metadata, sidecar; force=false)
        return (; path, provenance=sidecar, sha256=digest, reused=false)
    finally
        isfile(temporary) && rm(temporary)
        isfile(temporary_metadata) && rm(temporary_metadata)
    end
end

end
