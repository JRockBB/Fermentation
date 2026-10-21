$ErrorActionPreference = 'Stop'
$repoRoot = (Get-Location).Path
$fitHistory = Join-Path $repoRoot 'data/processed/cocoa_benchmark/fits'
$receiptPath = Join-Path $fitHistory 'bayes_brwb3_depth_pilot1.process.json'
$stdoutPath = Join-Path $fitHistory 'bayes_brwb3_depth_pilot1.stdout.log'
$stderrPath = Join-Path $fitHistory 'bayes_brwb3_depth_pilot1.stderr.log'
foreach ($newPath in @($receiptPath, $stdoutPath, $stderrPath, (Join-Path $fitHistory 'bayes_brwb3_depth_pilot1'))) {
    if (Test-Path -LiteralPath $newPath) { throw "Refusing to overwrite $newPath" }
}
$startedUtc = (Get-Date).ToUniversalTime().ToString('o')
$pilotProcess = Start-Process -FilePath 'C:/Users/jmurp/.julia/juliaup/julia-1.10.11+0.x64.w64.mingw32/bin/julia.exe' -ArgumentList @('--startup-file=no', '--check-bounds=yes', '--threads=2', '--project', '.julia-depot/p25_brwb3_depth_pilot.jl') -WorkingDirectory $repoRoot -RedirectStandardOutput $stdoutPath -RedirectStandardError $stderrPath -WindowStyle Hidden -PassThru
$retainedHandle = $pilotProcess.Handle
$receipt = [ordered]@{
    status = 'running'
    process_id = $pilotProcess.Id
    started_utc = $startedUtc
    watchdog_seconds = 3600
    watchdog_scope = 'Only the Julia process started by this wrapper; includes package loading and postprocessing.'
    exit_code_capture = 'Native process handle retained before WaitForExit; missing exit code fails the launcher.'
    stdout = 'bayes_brwb3_depth_pilot1.stdout.log'
    stderr = 'bayes_brwb3_depth_pilot1.stderr.log'
}
$receipt | ConvertTo-Json -Depth 4 | Set-Content -LiteralPath $receiptPath -Encoding UTF8
try {
    $completed = $pilotProcess.WaitForExit(3600000)
    if (-not $completed) {
        $receipt['status'] = 'watchdog_terminated'
        $pilotProcess.Kill()
        $pilotProcess.WaitForExit()
        $exitStatus = 124
    } else {
        $pilotProcess.WaitForExit()
        $exitStatus = $pilotProcess.ExitCode
        $receipt['status'] = if ($null -eq $exitStatus) { 'exit_code_unavailable' } else { 'exited' }
    }
    $receipt['exit_code'] = $exitStatus
    $receipt['native_exit_code'] = $pilotProcess.ExitCode
    $receipt['finished_utc'] = (Get-Date).ToUniversalTime().ToString('o')
    $receipt['stdout_sha256'] = (Get-FileHash -LiteralPath $stdoutPath -Algorithm SHA256).Hash.ToLowerInvariant()
    $receipt['stderr_sha256'] = (Get-FileHash -LiteralPath $stderrPath -Algorithm SHA256).Hash.ToLowerInvariant()
    $receipt | ConvertTo-Json -Depth 4 | Set-Content -LiteralPath $receiptPath -Encoding UTF8
    if ($null -eq $exitStatus) { throw 'Native exit code unavailable; receipt preserved without inferring success' }
} finally {
    $pilotProcess.Dispose()
}
exit $exitStatus
