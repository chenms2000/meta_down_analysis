param(
    [string]$ReleaseId = "",
    [switch]$SkipArticles,
    [ValidateSet("none", "abstract", "full")]
    [string]$ArticleSentenceScope = "abstract",
    [ValidateSet("rules", "scispacy", "hybrid")]
    [string]$SentenceParser = "rules",
    [string]$ScispacyModel = "en_core_sci_sm",
    [int]$MaxArticles = 0,
    [switch]$PubChemEnrich
)

$ErrorActionPreference = "Stop"

$argsList = @(
    ".\scripts\build_normalized_store.py",
    "--article-sentence-scope", $ArticleSentenceScope,
    "--sentence-parser", $SentenceParser,
    "--scispacy-model", $ScispacyModel
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
