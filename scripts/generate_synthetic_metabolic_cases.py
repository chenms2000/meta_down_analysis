"""Generate synthetic metabolite tables for rule and calibration stress tests.

Synthetic cases are pretraining, workflow, robustness, and coverage fixtures only.
They are not biological truth labels.
"""

from __future__ import annotations

import argparse
import json
import random
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


THEME_SEEDS: dict[str, list[tuple[str, str]]] = {
    "renal_organic_anion_uremic_toxin_handling": [
        ("indoxyl sulfate", "up"),
        ("p-cresol sulfate", "up"),
        ("hippuric acid", "up"),
        ("phenylacetylglutamine", "up"),
        ("TMAO", "up"),
        ("uric acid", "up"),
        ("creatinine", "up"),
    ],
    "acylcarnitine_fatty_acid_oxidation_pressure": [
        ("carnitine", "down"),
        ("acetylcarnitine", "down"),
        ("propionylcarnitine", "up"),
        ("butyrylcarnitine", "up"),
        ("palmitoylcarnitine", "up"),
        ("stearoylcarnitine", "up"),
    ],
    "arginine_no_metabolism": [("ADMA", "up"), ("arginine", "down"), ("nitrite", "down"), ("citrulline", "up"), ("ornithine", "up")],
    "purine_degradation_nucleotide_stress": [("hypoxanthine", "up"), ("xanthine", "up"), ("uric acid", "up"), ("adenosine", "down"), ("AMP", "down")],
    "idh_like_metabolic_pressure": [("2-hydroxyglutarate", "up"), ("alpha-ketoglutarate", "down")],
    "glial_glutamate_gaba_metabolism": [("glutamate", "up"), ("glutamine", "down"), ("GABA", "down"), ("N-acetylaspartate", "down"), ("myo-inositol", "up")],
    "glutathione_redox_stress": [("reduced glutathione", "down"), ("oxidized glutathione", "up"), ("cysteine", "down")],
    "nucleotide_sugar_glycosylation": [("UDP-glucose", "up"), ("UDP-galactose", "up"), ("UDP-N-acetylglucosamine", "up")],
    "phospholipid_membrane_remodeling": [("phosphatidylcholine", "down"), ("lysophosphatidylcholine", "up"), ("choline", "up")],
}

CONTEXTS = {
    "renal_tubular": {"organ": "kidney", "cell_type": "renal proximal tubular epithelial cell", "sample_type": "tissue"},
    "glioma_brain": {"disease": "glioma", "organ": "brain", "cell_type": "glial tumor cell", "sample_type": "tumor tissue"},
    "epithelial": {"cell_type": "epithelial cell", "experimental_context": "barrier injury and membrane repair"},
    "generalized": {},
}

NOISE = ["glucose", "lactate", "succinate", "taurine", "betaine", "tryptophan", "kynurenine", "sphingosine"]
ALIASES = {"p-cresol sulfate": ["p cresol sulfate", "p-cresyl sulfate"], "alpha-ketoglutarate": ["2-oxoglutarate", "AKG"], "TMAO": ["trimethylamine N-oxide"]}


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def perturb_name(name: str, rng: random.Random) -> str:
    choices = ALIASES.get(name, [])
    if choices and rng.random() < 0.35:
        return rng.choice(choices)
    if rng.random() < 0.12:
        return name.lower()
    return name


def direction_to_fc(direction: str, rng: random.Random) -> float:
    mag = rng.uniform(0.35, 1.8)
    return round(mag if direction == "up" else -mag, 4)


def build_case(index: int, rng: random.Random, generalized_ratio: float) -> dict[str, Any]:
    context_id = "generalized" if rng.random() < generalized_ratio else rng.choice(["renal_tubular", "glioma_brain", "epithelial"])
    if context_id == "renal_tubular":
        theme_pool = ["renal_organic_anion_uremic_toxin_handling", "acylcarnitine_fatty_acid_oxidation_pressure", "arginine_no_metabolism", "purine_degradation_nucleotide_stress"]
    elif context_id == "glioma_brain":
        theme_pool = ["glial_glutamate_gaba_metabolism", "idh_like_metabolic_pressure", "glutathione_redox_stress"]
    elif context_id == "epithelial":
        theme_pool = ["nucleotide_sugar_glycosylation", "phospholipid_membrane_remodeling"]
    else:
        theme_pool = list(THEME_SEEDS)
    themes = sorted(set(rng.sample(theme_pool, k=min(len(theme_pool), rng.randint(1, 3)))))
    records = []
    for theme in themes:
        seeds = THEME_SEEDS[theme]
        for name, direction in rng.sample(seeds, k=min(len(seeds), rng.randint(2, min(5, len(seeds))))):
            padj = rng.choice([0.005, 0.01, 0.03, 0.08])
            records.append({"metabolite": perturb_name(name, rng), "log2FC": direction_to_fc(direction, rng), "padj": padj})
    for name in rng.sample(NOISE, k=rng.randint(1, 4)):
        records.append({"metabolite": perturb_name(name, rng), "log2FC": direction_to_fc(rng.choice(["up", "down"]), rng), "padj": rng.choice([0.04, 0.2, 0.6])})
    rng.shuffle(records)
    return {
        "case_id": f"synthetic_{index:05d}",
        "synthetic_use": ["pretraining", "pipeline_test", "robustness_test", "coverage_test"],
        "not_truth_label": True,
        "mode": "generalized" if context_id == "generalized" else "context_aware",
        "context": CONTEXTS[context_id],
        "expected_themes": themes,
        "records": records,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate synthetic metabolic analysis cases.")
    parser.add_argument("--output", required=True)
    parser.add_argument("--cases", type=int, default=200)
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument("--generalized-ratio", type=float, default=0.35)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    rng = random.Random(args.seed)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8") as handle:
        for index in range(args.cases):
            handle.write(json.dumps(build_case(index, rng, args.generalized_ratio), ensure_ascii=False, sort_keys=True) + "\n")
    manifest = {
        "created_at_utc": utc_now(),
        "output": str(output),
        "cases": args.cases,
        "seed": args.seed,
        "note": "Synthetic cases are weak training and stress-test inputs, not final biological labels.",
    }
    (output.with_suffix(output.suffix + ".manifest.json")).write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
