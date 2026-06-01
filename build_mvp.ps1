param(
    [string]$ReleaseId = "",
    [switch]$SkipArticles,
    [ValidateSet("none", "abstract", "full")]
    [string]$ArticleSentenceScope = "abstract",
    [int]$MaxArticles = 0,
    [switch]$PubChemEnrich
)

$ErrorActionPreference = "Stop"

$argsList = @(
    ".\scripts\build_normalized_store.py",
    "--article-sentence-scope", $ArticleSentenceScope
)

if ($ReleaseId) {
    $argsList += @("--release-id", $ReleaseId)
}

if ($SkipArticles) {
    $argsList += "--skip-articles"
}

if ($MaxArticles -gt 0) {
    $argsList += @("--max-articles", [string]$MaxArticles)
}

if ($PubChemEnrich) {
    $argsList += "--pubchem-enrich"
}

python @argsList
