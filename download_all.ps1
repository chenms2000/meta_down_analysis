param(
    [switch]$Run,
    [switch]$ResolveIndexes,
    [switch]$IncludeDisabled,
    [switch]$AcceptLicensed,
    [string]$Sources = "",
    [string]$Layers = "open_core",
    [string]$OutputRoot = "raw_lake",
    [string]$ReleaseId = "",
    [int]$Jobs = 2,
    [int]$LimitPerSource = 0
)

$ErrorActionPreference = "Stop"

$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$argsList = @(
    "$scriptDir\scripts\download_databases.py",
    "--catalog", "$scriptDir\config\source_catalog.toml",
    "--output-root", "$scriptDir\$OutputRoot",
    "--manifest-dir", "$scriptDir\manifests",
    "--layers", $Layers,
    "--jobs", "$Jobs"
)

if ($ReleaseId -ne "") {
    $argsList += @("--release-id", $ReleaseId)
}
if ($Sources -ne "") {
    $argsList += @("--sources", $Sources)
}
if ($ResolveIndexes) {
    $argsList += "--resolve-indexes"
}
if ($IncludeDisabled) {
    $argsList += "--include-disabled"
}
if ($AcceptLicensed) {
    $argsList += "--accept-licensed"
}
if ($LimitPerSource -gt 0) {
    $argsList += @("--limit-per-source", "$LimitPerSource")
}

if ($Run) {
    $argsList += "--yes"
} else {
    $argsList += "--dry-run"
}

python @argsList
