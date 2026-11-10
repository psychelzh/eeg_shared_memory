# This local cleanup uses the reviewed 2026-11-10 inventory. Preview is the default.
# Run from PowerShell with -Apply only after reviewing logs/cache-cleanup-20261110.
[CmdletBinding()]
param([switch]$Apply)

$ErrorActionPreference = 'Stop'
$repoPath = [System.IO.Path]::GetFullPath((Split-Path -Parent $PSScriptRoot))
$auditRoot = Join-Path $repoPath 'logs/cache-cleanup-20261110'
$objectRoot = Join-Path $repoPath '_targets/objects'
$userRoot = Join-Path $repoPath '_targets/user'
$metaPath = Join-Path $repoPath '_targets/meta/meta'
$snapshot = Get-Content -LiteralPath (Join-Path $auditRoot 'audit-fingerprints.json') -Raw | ConvertFrom-Json

foreach ($entry in @($snapshot.Files)) {
    $actual = (Get-FileHash -LiteralPath $entry.Path -Algorithm SHA256).Hash
    if ($actual -ne $entry.Hash) {
        throw "The audited source or inventory changed; regenerate the audit before cleanup: $($entry.Path)"
    }
}

$keepNames = [System.Collections.Generic.HashSet[string]]::new([System.StringComparer]::Ordinal)
Import-Csv -LiteralPath (Join-Path $auditRoot 'keep_names.csv') | ForEach-Object { [void]$keepNames.Add($_.name) }
$references = [System.Collections.Generic.HashSet[string]]::new([System.StringComparer]::OrdinalIgnoreCase)
Import-Csv -LiteralPath (Join-Path $auditRoot 'referenced_files.csv') | ForEach-Object {
    [void]$references.Add([System.IO.Path]::GetFullPath($_.local))
}
$plan = @(Import-Csv -LiteralPath (Join-Path $auditRoot 'deletion-plan.csv'))
$totalBytes = ($plan | Measure-Object -Property Bytes -Sum).Sum
Write-Output ("Audited unused files: {0}; total: {1:N3} GiB" -f $plan.Count, ($totalBytes / 1GB))
if (-not $Apply) {
    $plan | Group-Object Kind | ForEach-Object {
        [pscustomobject]@{
            Kind = $_.Name
            Files = $_.Count
            GiB = [math]::Round(($_.Group | Measure-Object Bytes -Sum).Sum / 1GB, 3)
        }
    } | Format-Table
    Write-Output 'Preview only. No cache files or metadata were changed.'
    return
}

# Validate every absolute path and its current size before deleting any file.
foreach ($item in $plan) {
    $fullPath = [System.IO.Path]::GetFullPath($item.Path)
    $allowedRoot = switch ($item.Kind) {
        'object' { $objectRoot }
        'generated' { $userRoot }
        default { throw 'Unknown cleanup category.' }
    }
    $prefix = [System.IO.Path]::GetFullPath($allowedRoot) + [System.IO.Path]::DirectorySeparatorChar
    if (-not $fullPath.StartsWith($prefix, [System.StringComparison]::OrdinalIgnoreCase)) {
        throw "Cleanup path is outside the audited cache directory: $fullPath"
    }
    if ($item.Kind -eq 'object' -and $keepNames.Contains([System.IO.Path]::GetFileName($fullPath))) {
        throw "The plan includes a retained target: $fullPath"
    }
    if ($references.Contains($fullPath)) { throw "The plan includes a referenced file: $fullPath" }
    $current = Get-Item -LiteralPath $fullPath
    if ($current.PSIsContainer -or $current.Length -ne [long]$item.Bytes) {
        throw "The audited cache entry changed: $fullPath"
    }
}
foreach ($path in $references) {
    if ($path.StartsWith($userRoot + [System.IO.Path]::DirectorySeparatorChar, [System.StringComparison]::OrdinalIgnoreCase) -and -not (Test-Path -LiteralPath $path -PathType Leaf)) {
        throw "A retained generated file is missing: $path"
    }
}

$rscriptCommand = Get-Command Rscript -ErrorAction SilentlyContinue
$rscriptPath = if ($rscriptCommand) { $rscriptCommand.Source } else { 'C:/Program Files/R/R-4.6.2/bin/Rscript.exe' }
if (-not (Test-Path -LiteralPath $rscriptPath -PathType Leaf)) { throw 'Rscript was not found.' }
Copy-Item -LiteralPath $metaPath -Destination (Join-Path $auditRoot 'meta.before-prune')
$freeBefore = (Get-PSDrive -Name D).Free
$removed = 0
foreach ($item in $plan) {
    Remove-Item -LiteralPath $item.Path -Force
    $removed++
    if ($removed % 2000 -eq 0) { Write-Output "Removed $removed / $($plan.Count) verified unused files" }
}

Push-Location -LiteralPath $repoPath
try {
    $env:RENV_CONFIG_SANDBOX_ENABLED = 'FALSE'
    & $rscriptPath -e 'targets::tar_prune(cloud = FALSE, callr_function = NULL)'
    if ($LASTEXITCODE -ne 0) { throw 'Metadata pruning failed; the metadata backup is in the audit directory.' }
} finally {
    Pop-Location
}

# Check all retained data files against the pre-cleanup inventory.
foreach ($filename in @('retained-objects.csv', 'retained-generated.csv')) {
    foreach ($entry in (Import-Csv -LiteralPath (Join-Path $auditRoot $filename))) {
        $current = Get-Item -LiteralPath $entry.FullName
        if ($current.Length -ne [long]$entry.Length) { throw "A retained file changed: $($entry.FullName)" }
    }
}
[pscustomobject]@{
    Removed = $removed
    PlannedGiB = [math]::Round($totalBytes / 1GB, 3)
    FreeBeforeGiB = [math]::Round($freeBefore / 1GB, 3)
    FreeAfterGiB = [math]::Round((Get-PSDrive -Name D).Free / 1GB, 3)
    FinishedUtc = [DateTime]::UtcNow.ToString('o')
} | ConvertTo-Json | Set-Content -LiteralPath (Join-Path $auditRoot 'cleanup-completed.json') -Encoding utf8
Write-Output 'Cleanup and retained-file verification completed.'
