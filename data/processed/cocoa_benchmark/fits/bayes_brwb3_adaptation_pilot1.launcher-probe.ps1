# Follow-up bookkeeping probe. Does not alter or rerun the scientific pilot.
$ErrorActionPreference = 'Stop'
$destination = Join-Path $PSScriptRoot 'bayes_brwb3_adaptation_pilot1.launcher-probe.json'
if (Test-Path -LiteralPath $destination) { throw 'Refusing to overwrite the probe report' }
$results = @()
foreach ($expectedExitCode in @(0, 7)) {
    $arguments = '-NoProfile -NonInteractive -Command "Start-Sleep -Milliseconds 200; exit ' + $expectedExitCode + '"'
    $child = Start-Process -FilePath (Join-Path $PSHOME 'powershell.exe') -ArgumentList $arguments -WindowStyle Hidden -PassThru
    # Retain the native handle before waiting. Fail closed if no code is returned.
    $retainedHandle = $child.Handle
    try {
        if (-not $child.WaitForExit(10000)) {
            $child.Kill()
            $child.WaitForExit()
            throw 'The probe child exceeded its deadline'
        }
        $child.WaitForExit()
        $actualExitCode = $child.ExitCode
        if ($null -eq $actualExitCode -or $actualExitCode -ne $expectedExitCode) {
            throw "Expected exit $expectedExitCode; observed $actualExitCode"
        }
        $results += [ordered]@{ expected_exit_code = $expectedExitCode; observed_exit_code = $actualExitCode; passed = $true }
    } finally {
        $child.Dispose()
    }
}
$report = [ordered]@{
    purpose = 'Verify retaining a process handle before waiting captures both zero and nonzero child exit codes for future launchers.'
    native_exit_code_of_original_pilot = $null
    original_pilot_receipt_modified = $false
    script_sha256 = (Get-FileHash -LiteralPath $PSCommandPath -Algorithm SHA256).Hash.ToLowerInvariant()
    finished_utc = (Get-Date).ToUniversalTime().ToString('o')
    cases = $results
}
$report | ConvertTo-Json -Depth 5 | Set-Content -LiteralPath $destination -Encoding UTF8
$report | ConvertTo-Json -Depth 5
