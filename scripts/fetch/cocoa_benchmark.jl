# P2.5: acquire the archived 2022 benchmark without executing its R/Stan code.
include("download.jl")

raw = joinpath(@__DIR__, "..", "..", "data", "raw", "cocoa_benchmark")
const HASHES = Dict(
    "ecbfm-v1.0.zip" => "735255d908338acd573d993d6ab99963a567a3a81773bd910156c2a9754364c4",
    "rsos210274_si_001.pdf" => "54fbaa2b5942f2b4464fade8a3f9987476da1038c06eac10a6c35458048b6f92",
)
for (url, file, source, license) in (
    ("https://zenodo.org/api/records/5510370", "zenodo_5510370.json",
        "doi:10.5281/zenodo.5510370", "repository metadata"),
    ("https://zenodo.org/api/records/5510370/files/mmorenozam/ecbfm-v1.0.zip/content",
        "ecbfm-v1.0.zip", "doi:10.5281/zenodo.5510370", "other-open (Zenodo; no license file in archive)"),
    ("https://api.figshare.com/v2/articles/19149492/versions/1", "figshare_19149492_v1.json",
        "doi:10.6084/m9.figshare.19149492.v1", "repository metadata"),
    ("https://ndownloader.figshare.com/files/34022499", "rsos210274_si_001.pdf",
        "doi:10.6084/m9.figshare.19149492.v1", "CC-BY-4.0"),
    ("https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8847890/fullTextXML",
        "PMC8847890.xml", "doi:10.1098/rsos.210274", "CC-BY-4.0"),
)
    result = FetchData.fetch_file(url, joinpath(raw, file);
        expected_sha256=get(HASHES, file, nothing), source, license)
    println(result.reused ? "Verified: " : "Downloaded: ", file, " ", result.sha256)
end
