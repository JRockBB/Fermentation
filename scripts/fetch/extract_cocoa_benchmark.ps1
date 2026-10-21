# P2.5. Extract original bytes once; later calls verify every member.
$ErrorActionPreference = 'Stop'
Add-Type -AssemblyName System.IO.Compression.FileSystem
$repo = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot '../..'))
$raw = Join-Path $repo 'data/raw/cocoa_benchmark'
$archive = Join-Path $raw 'ecbfm-v1.0.zip'
$expected = '735255d908338acd573d993d6ab99963a567a3a81773bd910156c2a9754364c4'
if ((Get-FileHash -LiteralPath $archive -Algorithm SHA256).Hash.ToLowerInvariant() -ne $expected) {
    throw 'Benchmark archive checksum differs'
}
$zip = [IO.Compression.ZipFile]::OpenRead($archive)
try {
    $records = @()
    foreach ($entry in $zip.Entries) {
        $target = [IO.Path]::GetFullPath((Join-Path $raw $entry.FullName))
        $prefix = [IO.Path]::GetFullPath($raw) + [IO.Path]::DirectorySeparatorChar
        if (-not $target.StartsWith($prefix, [StringComparison]::OrdinalIgnoreCase)) {
            throw 'Archive member escapes destination'
        }
        if ($entry.FullName.EndsWith('/')) { continue }
        $stream = $entry.Open()
        $sha = [Security.Cryptography.SHA256]::Create()
        try { $hash = [BitConverter]::ToString($sha.ComputeHash($stream)).Replace('-', '').ToLowerInvariant() }
        finally { $stream.Dispose(); $sha.Dispose() }
        if (Test-Path -LiteralPath $target) {
            if ((Get-FileHash -LiteralPath $target -Algorithm SHA256).Hash.ToLowerInvariant() -ne $hash) {
                throw "Existing member differs: $target"
            }
        } else {
            [IO.Directory]::CreateDirectory([IO.Path]::GetDirectoryName($target)) | Out-Null
            [IO.Compression.ZipFileExtensions]::ExtractToFile($entry, $target, $false)
            if ((Get-FileHash -LiteralPath $target -Algorithm SHA256).Hash.ToLowerInvariant() -ne $hash) {
                throw "Extracted member differs: $target"
            }
        }
        $records += "$($entry.FullName),$($entry.Length),$hash"
    }
    $manifest = Join-Path $repo 'docs/cocoa_benchmark_archive_files.csv'
    $text = (@('path,bytes,sha256') + $records) -join "`n"
    if (Test-Path -LiteralPath $manifest) {
        if ([IO.File]::ReadAllText($manifest).TrimEnd() -ne $text) { throw 'Member manifest differs' }
    } else {
        [IO.File]::WriteAllText($manifest, $text + "`n", [Text.UTF8Encoding]::new($false))
    }
    Write-Output "Verified $($records.Count) archive members."
} finally { $zip.Dispose() }
