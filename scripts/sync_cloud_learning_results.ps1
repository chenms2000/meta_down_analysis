param(
    [string]$RemoteUser = "ms",
    [string]$RemoteHost = "10.64.2.166",
    [string]$RemoteWorkspace = "/disk3/ms/meta_down_analysis_learning_scibio",
    [string]$RunId = "learn_scibio_20260605T165940Z",
    [string]$LocalWorkspace = (Resolve-Path ".").Path,
    [string]$LocalResultsRoot = "",
    [switch]$Overwrite,
    [switch]$KeepArchive,
    [switch]$LegacyScp
)

$ErrorActionPreference = "Stop"

function Require-Command {
    param([string]$Name)
    if (-not (Get-Command $Name -ErrorAction SilentlyContinue)) {
        throw "Required command '$Name' was not found. Install OpenSSH Client or add it to PATH."
    }
}

Require-Command ssh
Require-Command scp
Require-Command tar

function Invoke-Native {
    param(
        [string]$Command,
        [string[]]$Arguments
    )
    & $Command @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "Command failed with exit code ${LASTEXITCODE}: $Command $($Arguments -join ' ')"
    }
}

if ([string]::IsNullOrWhiteSpace($LocalResultsRoot)) {
    $LocalResultsRoot = Join-Path $LocalWorkspace "learning_runs"
}

$Remote = "$RemoteUser@$RemoteHost"
$RemoteRunPath = "$RemoteWorkspace/learning_runs/$RunId"
$RemoteArchive = "/tmp/${RunId}.tgz"
$LocalResultsRoot = [System.IO.Path]::GetFullPath($LocalResultsRoot)
$LocalRunPath = Join-Path $LocalResultsRoot $RunId
$LocalTmp = Join-Path $LocalWorkspace ".tmp\cloud-downloads"
$LocalArchive = Join-Path $LocalTmp "${RunId}.tgz"
$ManifestDir = Join-Path $LocalWorkspace "manifests\cloud_learning_sync"
$ManifestPath = Join-Path $ManifestDir "${RunId}.sync.json"

if ((Test-Path $LocalRunPath) -and -not $Overwrite) {
    throw "Local run already exists: $LocalRunPath. Re-run with -Overwrite to replace the extracted copy."
}

New-Item -ItemType Directory -Force -Path $LocalResultsRoot, $LocalTmp, $ManifestDir | Out-Null

Write-Host "Checking remote run: ${Remote}:$RemoteRunPath"
Invoke-Native ssh @($Remote, "test -d '$RemoteRunPath' && test -f '$RemoteRunPath/reports/priority_ranker_report.json'")

Write-Host "Creating remote archive: $RemoteArchive"
Invoke-Native ssh @($Remote, "rm -f '$RemoteArchive' && tar -C '$RemoteWorkspace/learning_runs' -czf '$RemoteArchive' '$RunId'")

Write-Host "Downloading archive to: $LocalArchive"
$ScpArgs = @()
if ($LegacyScp) {
    $ScpArgs += "-O"
}
$ScpArgs += "${Remote}:$RemoteArchive"
$ScpArgs += "$LocalArchive"
Invoke-Native scp $ScpArgs

if ((Test-Path $LocalRunPath) -and $Overwrite) {
    Remove-Item -LiteralPath $LocalRunPath -Recurse -Force
}

Write-Host "Extracting into: $LocalResultsRoot"
Invoke-Native tar @("-xzf", "$LocalArchive", "-C", "$LocalResultsRoot")

$ExpectedFiles = @(
    "rankings\drug_priorities.parquet",
    "rankings\pathway_priorities.parquet",
    "rankings\relation_priorities.parquet",
    "rankings\target_priorities.parquet",
    "rankings\entity_embeddings.parquet",
    "rankings\entity_neighbors.parquet",
    "rankings\entity_clusters.parquet",
    "rankings\entity_cluster_assignments.parquet",
    "weak_labels\weak_labels.parquet",
    "models\priority_ranker.joblib",
    "reports\priority_ranker_report.json",
    "reports\global_learning_report.md",
    "reports\global_top_pathways.csv",
    "reports\global_top_targets.csv",
    "reports\global_top_drugs.csv",
    "reports\global_top_relations.csv"
)

$Missing = @()
foreach ($RelativePath in $ExpectedFiles) {
    $Path = Join-Path $LocalRunPath $RelativePath
    if (-not (Test-Path $Path)) {
        $Missing += $RelativePath
    }
}

$ArchiveHash = Get-FileHash -Algorithm SHA256 -Path $LocalArchive
$Manifest = [ordered]@{
    run_id = $RunId
    remote = "${Remote}:$RemoteRunPath"
    local_run_path = $LocalRunPath
    downloaded_archive = $LocalArchive
    archive_sha256 = $ArchiveHash.Hash
    synced_at_utc = (Get-Date).ToUniversalTime().ToString("o")
    expected_files = $ExpectedFiles
    missing_expected_files = $Missing
    note = "Learning outputs are research-priority artifacts. They are ignored by Git and should not be interpreted as clinical recommendations."
}
$Manifest | ConvertTo-Json -Depth 5 | Set-Content -Encoding UTF8 -Path $ManifestPath

if (-not $KeepArchive) {
    Remove-Item -LiteralPath $LocalArchive -Force
    $Manifest["downloaded_archive"] = $null
    $Manifest | ConvertTo-Json -Depth 5 | Set-Content -Encoding UTF8 -Path $ManifestPath
}

Write-Host ""
Write-Host "Downloaded run: $LocalRunPath"
Write-Host "Sync manifest: $ManifestPath"

if ($Missing.Count -gt 0) {
    Write-Warning ("Missing expected files: " + ($Missing -join ", "))
    exit 2
}

Write-Host "All expected learning outputs are present."
