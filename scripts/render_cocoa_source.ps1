param(
    [Parameter(Mandatory=$true)][string]$InputPdf,
    [Parameter(Mandatory=$true)][string]$OutputDirectory,
    [int]$Width = 6000
)
$ErrorActionPreference = 'Stop'
$Pages = @(44,45,46)
$expectedHash = 'D5E63C6ECA4C2A82FEA9E13F3BB79D456681DFEA8A567637632E5613370A6A99'
if ((Get-FileHash -LiteralPath $InputPdf -Algorithm SHA256).Hash -ne $expectedHash) { throw 'Unexpected thesis PDF checksum' }
$outputRoot = $ExecutionContext.SessionState.Path.GetUnresolvedProviderPathFromPSPath($OutputDirectory)
$null = [System.IO.Directory]::CreateDirectory($outputRoot)
# Windows display scaling can affect the rendered size. The digitizer requires
# 9000 x 12729; adjust Width proportionally if another display scale differs.
Add-Type -AssemblyName System.Runtime.WindowsRuntime
$null = [Windows.Data.Pdf.PdfDocument, Windows.Data.Pdf, ContentType=WindowsRuntime]
$null = [Windows.Storage.StorageFile, Windows.Storage, ContentType=WindowsRuntime]
$null = [Windows.Storage.Streams.InMemoryRandomAccessStream, Windows.Storage.Streams, ContentType=WindowsRuntime]
$null = [Windows.Data.Pdf.PdfPageRenderOptions, Windows.Data.Pdf, ContentType=WindowsRuntime]
$asTaskMethod = [System.WindowsRuntimeSystemExtensions].GetMethods() | Where-Object {
    $_.Name -eq 'AsTask' -and $_.IsGenericMethod -and $_.GetParameters().Count -eq 1 -and
    $_.GetParameters()[0].ParameterType.Name -eq 'IAsyncOperation`1'
} | Select-Object -First 1
$asActionMethod = [System.WindowsRuntimeSystemExtensions].GetMethods() | Where-Object {
    $_.Name -eq 'AsTask' -and -not $_.IsGenericMethod -and $_.GetParameters().Count -eq 1
} | Select-Object -First 1
$asStreamMethod = [System.IO.WindowsRuntimeStreamExtensions].GetMethods() | Where-Object {
    $_.Name -eq 'AsStreamForRead' -and $_.GetParameters().Count -eq 1
} | Select-Object -First 1
function Await-WinRT($Operation, [Type]$ResultType) {
    $task = $asTaskMethod.MakeGenericMethod($ResultType).Invoke($null,@($Operation))
    return $task.GetAwaiter().GetResult()
}
$inputPath = (Resolve-Path -LiteralPath $InputPdf).Path
$file = Await-WinRT ([Windows.Storage.StorageFile]::GetFileFromPathAsync($inputPath)) ([Windows.Storage.StorageFile])
$pdf = Await-WinRT ([Windows.Data.Pdf.PdfDocument]::LoadFromFileAsync($file)) ([Windows.Data.Pdf.PdfDocument])
Write-Output "PDF pages: $($pdf.PageCount)"
foreach ($number in $Pages) {
    if ($number -ge $pdf.PageCount) { continue }
    $page = $pdf.GetPage([uint32]$number)
    $stream = [Windows.Storage.Streams.InMemoryRandomAccessStream]::new()
    $options = [Windows.Data.Pdf.PdfPageRenderOptions]::new()
    $options.DestinationWidth = $Width
    try {
        $renderTask = $asActionMethod.Invoke($null,@($page.RenderToStreamAsync($stream,$options)))
        $renderTask.GetAwaiter().GetResult()
        $reader = $asStreamMethod.Invoke($null,@($stream.GetInputStreamAt(0)))
        $outputPath = Join-Path $outputRoot ('cocoa-thesis-page-' + ($number+1) + '.png')
        $writer = [System.IO.File]::Open($outputPath,[System.IO.FileMode]::CreateNew)
        try { $reader.CopyTo($writer) } finally { $writer.Dispose(); $reader.Dispose() }
        Write-Output $outputPath
    } finally { $stream.Dispose(); $page.Dispose() }
}
