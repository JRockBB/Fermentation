include("download.jl")

result = FetchData.fetch_file(
    "https://ndownloader.figshare.com/files/13237406",
    joinpath(@__DIR__,"..","..","data","raw","cocoa","rsos180964_si_001.pdf");
    # Figshare article 7152614 v2; publisher MD5 d2a400dadd0768b34afe64b6b9d05023.
    expected_sha256="10a337b696e50ca417e915ac90e22064e04a65533b89fb3300fc5b7542084a85",
    source="doi:10.6084/m9.figshare.7152614.v2",
    license="CC-BY-4.0",
)
println(result.reused ? "Verified existing supplement: " : "Downloaded supplement: ",result.path)
println("SHA-256: ",result.sha256)
