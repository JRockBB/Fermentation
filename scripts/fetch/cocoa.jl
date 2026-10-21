include("download.jl")

result = FetchData.fetch_file(
    "https://zenodo.org/records/4935403/files/CocoaModel.zip?download=1",
    joinpath(@__DIR__, "..", "..", "data", "raw", "cocoa", "CocoaModel.zip");
    # SHA-256 of the archive whose MD5 matches Zenodo's 39dc2c73f73301be8ec48096fb25e319.
    expected_sha256="9cc6c82ae32c7360408d0bef072e80462df882a004d7aa8975c52f301568c40a",
    source="doi:10.5061/dryad.321d33v",
    license="CC0-1.0",
)
println(result.reused ? "Verified existing archive: " : "Downloaded archive: ", result.path)
println("SHA-256: ", result.sha256)
