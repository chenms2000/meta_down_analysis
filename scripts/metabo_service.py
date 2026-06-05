"""Read-only MVP service layer over frozen metabolism-oncology releases.

The service intentionally keeps the online layer deterministic: all responses are
computed from parquet artifacts in normalized_store, graph_projection, and the
PubChem CID cache for a selected release.
"""

from __future__ import annotations

import argparse
import csv
import copy
import gzip
import hashlib
import heapq
import io
import json
import math
import re
import sys
import time
from collections import Counter, defaultdict
from dataclasses import dataclass
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any, Iterable
from urllib.parse import parse_qs, unquote, urlparse

try:
    import pyarrow as pa
    import pyarrow.dataset as ds
    import pyarrow.parquet as pq
except Exception:  # pragma: no cover - handled by require_arrow.
    pa = None
    ds = None
    pq = None

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from build_graph_projection import (  # noqa: E402
    latest_normalized_release_id,
    normalize_id_token,
    normalize_lookup_key,
)
from build_compound_match_index import (  # noqa: E402
    DEFAULT_OUTPUT_ROOT as DEFAULT_COMPOUND_ROOT,
    inchikey_connectivity,
    normalize_formula,
)
from resolve_entity import candidate_queries  # noqa: E402
from llm_safe_adapter import (  # noqa: E402
    ADAPTER_VERSION as LOCAL_LLM_ADAPTER_VERSION,
    EXTERNAL_ADAPTER_VERSION,
    ExternalLLMConfig,
    ExternalLLMChunkedTextNarratorBackend,
    ExternalLLMNarratorBackend,
    ExternalLLMTextNarratorBackend,
    INPUT_CONTRACT_VERSION as LLM_ADAPTER_INPUT_CONTRACT_VERSION,
    LLMAdapterError,
    OUTPUT_CONTRACT_VERSION as LLM_ADAPTER_OUTPUT_CONTRACT_VERSION,
    build_local_adapter_output,
    default_openai_compatible_transport,
    extract_chat_completion_content,
    guard_adapter_output,
    guard_policy_hash,
)


DEFAULT_NORMALIZED_ROOT = "normalized_store"
DEFAULT_GRAPH_ROOT = "graph_projection"
DEFAULT_PUBCHEM_ROOT = "pubchem_cid_cache"
DEFAULT_COMPOUND_MATCH_ROOT = DEFAULT_COMPOUND_ROOT
DEFAULT_LITERATURE_ROOT = "literature_evidence"
DEFAULT_PREDICTION_OVERLAY_ROOT = "manual_sources/prediction_overlays"
DEFAULT_DATABASE_ACCURACY_ROOT = "database_accuracy_store"
DEFAULT_GOLD_STANDARD_PATHS = (
    "config/gold_standard_conclusions.json",
    "manual_sources/gold_standard_conclusions.json",
    "manual_sources/gold_standard_conclusions.csv",
)
EXTERNAL_LLM_CONFIG_ERROR_CODES = {
    "external_llm_disabled",
    "invalid_llm_api_key_characters",
    "missing_llm_api_key",
    "missing_llm_endpoint",
    "missing_llm_model",
    "unsupported_llm_provider",
}
DEFAULT_MIN_MATCH_SCORE = 60.0
DEFAULT_MIN_MATCH_MARGIN = 15.0
EXPANDED_SEED_WEIGHTS = {
    "strict_identity": 1.0,
    "soft_identity": 0.5,
    "analog_candidate": 0.3,
    "class_or_pool": 0.2,
    "ratio_component": 0.25,
    "unresolved": 0.0,
}
DEFAULT_MAX_CANDIDATES = 10
DEFAULT_MAX_PATHS = 25
DEFAULT_MAX_HOPS = 4
DEFAULT_MAX_EDGES_PER_NODE = 40
DEFAULT_PROPAGATION_MAX_ADAPTIVE_EDGES_PER_NODE = 160
DEFAULT_MAX_SUBGRAPH_EDGES = 500
DEFAULT_MAX_EVIDENCE_ITEMS = 25
DEFAULT_PPM_TOLERANCE = 10.0
DEFAULT_RT_TOLERANCE = 0.25
DEFAULT_MS2_MZ_TOLERANCE = 0.02
DEFAULT_PROPAGATION_ALPHA = 0.85
DEFAULT_PROPAGATION_ITERATIONS = 30
DEFAULT_PROPAGATION_TOLERANCE = 1e-9
DEFAULT_MAX_PROPAGATION_NODES = 1000
DEFAULT_MAX_PROPAGATION_EDGES = 5000
DEFAULT_PROPAGATION_BEAM_PER_TYPE = 40
DEFAULT_PROPAGATION_RETAINED_MASS = 0.95
DEFAULT_MIN_LITERATURE_OVERLAY_PROB = 0.65
DEFAULT_MAX_PREDICTION_ROWS = 10
ANALYSIS_PACK_CONTRACT_VERSION = "analysis_pack.v1"
DATABASE_ACCURACY_CONTRACT_VERSION = "database_accuracy_store.v2"
PREDICTION_PACK_CONTRACT_VERSION = "prediction_pack.v1"
INTERPRETATION_REPORT_CONTRACT_VERSION = "interpretation_report.v1"
CONCLUSION_EVALUATION_CONTRACT_VERSION = "conclusion_evaluation.v1"
EXPLAIN_CONTRACT_VERSION = "explain_endpoint.v1"
CHAT_CONTRACT_VERSION = "chat_endpoint.v1"
ANALYSIS_PACK_RANKING_LIMIT = 25
ANALYSIS_PACK_PATH_LIMIT = 10
CONCLUSION_OBSERVATION_CANDIDATE_LIMIT = 40
GOLD_GENERIC_TERMS = {
    "acid",
    "alpha",
    "axis",
    "class",
    "compound",
    "control",
    "coa",
    "cycle",
    "exact",
    "fatty",
    "identity",
    "metabolism",
    "metabolic",
    "negative",
    "positive",
    "seed",
    "trait",
}
BROAD_CONTEXT_TERMS = {
    "adjacent",
    "cancer",
    "carcinoma",
    "cell",
    "cells",
    "control",
    "noncancerous",
    "normal",
    "tumor",
    "tumour",
}

COMMON_METABOLITE_NAME_ALIASES = {
    "glucose": ["D-glucose", "D-glucopyranose", "alpha-D-glucose", "beta-D-glucose"],
    "glucose 6-phosphate": ["D-glucopyranose 6-phosphate", "glucose-6-phosphate", "glucose 6 phosphate"],
    "fructose 6-phosphate": ["fructose-6-phosphate", "fructose 6 phosphate"],
    "fructose 1,6-bisphosphate": ["D-fructose 1,6-bisphosphate", "fructose-1,6-bisphosphate"],
    "lactate": ["L-lactic acid", "L-lactate", "(S)-lactic acid"],
    "lactic acid": ["L-lactic acid", "L-lactate"],
    "glutamine": ["L-glutamine"],
    "glutamate": ["L-glutamic acid", "L-glutamate"],
    "pyruvate": ["pyruvic acid"],
    "citrate": ["citric acid", "citrate(3-)"],
    "succinate": ["succinic acid", "succinate(2-)"],
    "malate": ["L-malic acid", "malic acid"],
    "fumarate": ["fumaric acid"],
    "alanine": ["L-alanine"],
    "serine": ["L-serine"],
    "glycine": ["glycine"],
    "aspartate": ["L-aspartic acid", "L-aspartate"],
    "asparagine": ["L-asparagine"],
    "cysteine": ["L-cysteine"],
    "oxidized glutathione": ["glutathione disulfide", "GSSG"],
    "methionine sulfoxide": ["methionine S-oxide"],
    "sorbitol": ["D-glucitol", "D-sorbitol"],
    "free carnitine": ["L-carnitine", "carnitine", "(R)-carnitine"],
    "udp-glucose": ["Uridine 5'-diphosphoglucose", "UDP-D-glucose", "UDP-glucose"],
    "udp-galactose": ["UDP-D-galactose", "UDP-galactose"],
    "udp-n-acetylglucosamine": ["UDP-N-acetyl-alpha-D-glucosamine", "UDP-GlcNAc"],
    "uridine monophosphate": ["uridine 5'-monophosphate", "UMP"],
    "nad+": ["NAD(+)", "nicotinamide adenine dinucleotide"],
    "acetylcarnitine": ["O-acetyl-L-carnitine", "acetyl-L-carnitine"],
    "butyrylcarnitine": ["O-butanoyl-L-carnitine", "butyryl-L-carnitine"],
    "palmitoylcarnitine": ["O-palmitoyl-L-carnitine", "palmitoyl-L-carnitine"],
    "ceramide d18:1/16:0": ["N-hexadecanoylsphingosine", "C16 ceramide"],
    "c16 ceramide": ["N-hexadecanoylsphingosine"],
    "asymmetric dimethylarginine": ["N(omega),N(omega)-dimethyl-L-arginine", "ADMA"],
    "adma": ["N(omega),N(omega)-dimethyl-L-arginine"],
    "xanthine": ["7H-xanthine"],
}

COMMON_BIOCHEMICAL_ROOT_TERMS = frozenset(
    {
        "acetate",
        "acetyl",
        "acetylcarnitine",
        "acetylglucosamine",
        "adenine",
        "adenosine",
        "adenosylhomocysteine",
        "adenosylmethionine",
        "adp",
        "alanine",
        "amp",
        "arachidonic",
        "arginine",
        "asparagine",
        "aspartate",
        "atp",
        "betaine",
        "carnitine",
        "ceramide",
        "choline",
        "citrate",
        "citrulline",
        "coa",
        "creatine",
        "creatinine",
        "cysteine",
        "cystine",
        "cytidine",
        "ctp",
        "cdp",
        "cmp",
        "ethanolamine",
        "fructose",
        "fucose",
        "fumarate",
        "galactose",
        "glcnac",
        "glucosamine",
        "glucose",
        "glutamate",
        "glutamine",
        "glutathione",
        "glycerol",
        "glycerophosphocholine",
        "glycine",
        "gsh",
        "gssg",
        "gtp",
        "gdp",
        "gmp",
        "guanosine",
        "histidine",
        "homocysteine",
        "hydroxyproline",
        "hypotaurine",
        "inosine",
        "inositol",
        "isoleucine",
        "isocitrate",
        "ketoglutarate",
        "kynurenine",
        "lactate",
        "leukotriene",
        "leucine",
        "lysine",
        "malate",
        "methionine",
        "myo-inositol",
        "nad",
        "nadh",
        "nadp",
        "nadph",
        "neuraminic",
        "nicotinamide",
        "ornithine",
        "oxaloacetate",
        "palmitoylcarnitine",
        "phenylalanine",
        "phosphocholine",
        "phosphoethanolamine",
        "proline",
        "prostaglandin",
        "prpp",
        "putrescine",
        "pyruvate",
        "ribose",
        "sam",
        "sah",
        "serine",
        "sphinganine",
        "sphingomyelin",
        "sphingosine",
        "spermidine",
        "spermine",
        "succinate",
        "taurine",
        "threonine",
        "tryptophan",
        "tyrosine",
        "udp",
        "ump",
        "uridine",
        "utp",
        "valine",
    }
)

COMMON_BIOCHEMICAL_MODIFIER_TERMS = frozenset(
    {
        "acid",
        "acyl",
        "adenine",
        "adenosine",
        "carnitine",
        "coa",
        "diphosphate",
        "glucosamine",
        "glutathione",
        "monophosphate",
        "nucleotide",
        "phosphate",
        "phospho",
        "prostaglandin",
        "sulfate",
        "sulfoxide",
        "triphosphate",
        "uridine",
    }
)

PATHWAY_THEME_TERMS = {
    "glucose": ["glycolysis", "glucose", "gluconeogenesis", "hexokinase", "glut"],
    "d-glucose": ["glycolysis", "glucose", "gluconeogenesis", "hexokinase", "glut"],
    "d-glucopyranose": ["glycolysis", "glucose", "gluconeogenesis", "hexokinase", "glut"],
    "lactate": ["lactate", "lactic", "pyruvate", "monocarboxylate", "slc16"],
    "l-lactate": ["lactate", "lactic", "pyruvate", "monocarboxylate", "slc16"],
    "lactic acid": ["lactate", "lactic", "pyruvate", "monocarboxylate", "slc16"],
    "l-lactic acid": ["lactate", "lactic", "pyruvate", "monocarboxylate", "slc16"],
    "glutamine": ["glutamine", "glutamate", "glutaminolysis", "gfpt"],
    "l-glutamine": ["glutamine", "glutamate", "glutaminolysis", "gfpt"],
    "glutamine zwitterion": ["glutamine", "glutamate", "glutaminolysis", "gfpt"],
}

BIOLOGICAL_ENTITY_POOLS = {
    "succinate_pool": {
        "label": "Succinate / dicarboxylate pool",
        "terms": ("succinate", "succinic acid", "succinate 2", "fumarate", "fumaric acid", "malate", "malic acid"),
    },
    "glutamate_glutamine_pool": {
        "label": "Glutamate / glutamine pool",
        "terms": ("glutamate", "glutamic acid", "glutamine", "glutaminolysis"),
    },
    "cysteine_homocysteine_pool": {
        "label": "Cysteine / homocysteine sulfur amino-acid pool",
        "terms": (
            "cysteine",
            "cystine",
            "cystathionine",
            "homocysteine",
            "methionine",
            "adenosylhomocysteine",
            "adenosylmethionine",
            "s-adenosylhomocysteine",
            "s-adenosylmethionine",
        ),
    },
    "glutathione_redox_pair": {
        "label": "Glutathione redox pair",
        "terms": ("glutathione", "oxidized glutathione", "reduced glutathione", "gsh", "gssg"),
    },
    "udp_sugar_pool": {
        "label": "UDP-sugar pool",
        "terms": (
            "udp glucose",
            "udp-glucose",
            "udp galactose",
            "udp-galactose",
            "udp n acetylglucosamine",
            "udp-n-acetylglucosamine",
            "udp glcnac",
            "udp-glcnac",
            "udp glucuronic",
            "udp-glucuronic",
        ),
    },
    "phospholipid_membrane_pool": {
        "label": "Phospholipid membrane-remodeling pool",
        "terms": (
            "phosphatidylcholine",
            "lysophosphatidylcholine",
            "phosphatidylethanolamine",
            "lysophosphatidylethanolamine",
            "phosphatidylserine",
            "phosphatidylinositol",
            "phosphocholine",
            "glycerophosphocholine",
            "choline",
            "ethanolamine",
            "plasmalogen",
        ),
    },
    "sphingolipid_membrane_pool": {
        "label": "Sphingolipid membrane pool",
        "terms": ("ceramide", "sphingosine", "sphinganine", "sphingomyelin", "glycosphingolipid"),
    },
    "glycolysis_glucose_pool": {
        "label": "Glucose utilization / glycolysis pool",
        "terms": (
            "glucose",
            "glucose 6 phosphate",
            "glucose-6-phosphate",
            "fructose 6 phosphate",
            "fructose-6-phosphate",
            "pyruvate",
            "pyruvic acid",
            "lactate",
            "lactic acid",
        ),
    },
    "osmolyte_epithelial_pool": {
        "label": "Osmolyte / epithelial homeostasis pool",
        "terms": ("taurine", "hypotaurine", "betaine", "myo inositol", "myo-inositol", "sorbitol", "glycerophosphocholine"),
    },
    "epithelial_barrier_repair_pool": {
        "label": "Epithelial barrier and membrane-repair pool",
        "terms": (
            "udp n acetylglucosamine",
            "udp-n-acetylglucosamine",
            "n acetylglucosamine",
            "n-acetylglucosamine",
            "glucosamine",
            "choline",
            "phosphocholine",
            "glycerophosphocholine",
            "sphingosine",
            "ceramide",
            "proline",
            "hydroxyproline",
        ),
    },
    "renal_organic_anion_uremic_toxin_pool": {
        "label": "Renal organic-anion / uremic-toxin handling pool",
        "terms": (
            "indoxyl sulfate",
            "indoxyl sulphate",
            "p cresol sulfate",
            "p-cresol sulfate",
            "p-cresyl sulfate",
            "hippuric acid",
            "hippurate",
            "phenylacetylglutamine",
            "tmao",
            "trimethylamine n oxide",
            "trimethylamine-n-oxide",
            "uric acid",
            "urate",
            "creatinine",
        ),
    },
    "acylcarnitine_fatty_acid_oxidation_pool": {
        "label": "Carnitine / acylcarnitine fatty-acid oxidation pool",
        "terms": (
            "carnitine",
            "acetylcarnitine",
            "propionylcarnitine",
            "butyrylcarnitine",
            "palmitoylcarnitine",
            "stearoylcarnitine",
            "long chain acylcarnitine",
            "long-chain acylcarnitine",
        ),
    },
    "arginine_no_pool": {
        "label": "Arginine / nitric-oxide bioavailability pool",
        "terms": ("arginine", "citrulline", "ornithine", "adma", "asymmetric dimethylarginine", "nitrite", "nitrate"),
    },
    "purine_degradation_pool": {
        "label": "Purine degradation / nucleotide stress pool",
        "terms": ("hypoxanthine", "xanthine", "uric acid", "urate", "adenosine", "amp", "adenosine monophosphate"),
    },
    "idh_like_pool": {
        "label": "2-hydroxyglutarate / alpha-ketoglutarate IDH-like pressure pool",
        "terms": (
            "2 hydroxyglutarate",
            "2-hydroxyglutarate",
            "d 2 hydroxyglutarate",
            "l 2 hydroxyglutarate",
            "alpha ketoglutarate",
            "alpha-ketoglutarate",
            "akg",
            "2 oxoglutarate",
            "2-oxoglutarate",
        ),
    },
    "glial_glutamate_gaba_pool": {
        "label": "Glial glutamate / GABA / osmolyte pool",
        "terms": ("glutamate", "glutamine", "gaba", "gamma aminobutyric", "n acetylaspartate", "n-acetylaspartate", "myo inositol", "myo-inositol"),
    },
}

BIOCHEMICAL_THEME_SPECS = {
    "amino_acid_redox_stress": {
        "input_terms": (
            "glutamate",
            "glutamine",
            "cysteine",
            "cystine",
            "methionine",
            "homocysteine",
            "taurine",
            "glutathione",
        ),
        "pool_ids": ("glutamate_glutamine_pool", "cysteine_homocysteine_pool", "glutathione_redox_pair"),
        "pathway_terms": ("amino acid", "glutamine", "glutamate", "cysteine", "redox", "oxidative stress"),
    },
    "sulfur_one_carbon_metabolism": {
        "input_terms": (
            "cystathionine",
            "homocysteine",
            "methionine",
            "adenosylhomocysteine",
            "adenosylmethionine",
            "taurine",
            "hypotaurine",
            "methylthioadenosine",
        ),
        "pool_ids": ("cysteine_homocysteine_pool", "osmolyte_epithelial_pool"),
        "pathway_terms": ("methionine", "cysteine", "homocysteine", "sulfur amino", "taurine", "one carbon", "transsulfuration"),
    },
    "nucleotide_sugar_glycosylation": {
        "input_terms": (
            "udp glucose",
            "udp-glucose",
            "udp galactose",
            "udp-galactose",
            "udp n acetylglucosamine",
            "udp-n-acetylglucosamine",
            "udp glcnac",
            "udp-glcnac",
            "glycosylation",
            "glycan",
        ),
        "pool_ids": ("udp_sugar_pool", "epithelial_barrier_repair_pool"),
        "pathway_terms": ("glycosylation", "glycan", "hexosamine", "udp sugar", "nucleotide sugar", "o-linked", "n-linked"),
    },
    "phospholipid_membrane_remodeling": {
        "input_terms": (
            "phosphatidylcholine",
            "lysophosphatidylcholine",
            "phosphatidylethanolamine",
            "lysophosphatidylethanolamine",
            "phosphocholine",
            "glycerophosphocholine",
            "choline",
            "ethanolamine",
            "sphingomyelin",
            "ceramide",
        ),
        "pool_ids": ("phospholipid_membrane_pool", "sphingolipid_membrane_pool"),
        "pathway_terms": ("phospholipid", "membrane", "choline", "ethanolamine", "sphingolipid", "lipid remodeling"),
    },
    "glutathione_redox_stress": {
        "input_terms": ("glutathione", "oxidized glutathione", "reduced glutathione", "gsh", "gssg", "cysteine", "cystine"),
        "pool_ids": ("glutathione_redox_pair", "cysteine_homocysteine_pool"),
        "pathway_terms": ("glutathione", "redox", "oxidative stress", "reactive oxygen", "cysteine"),
    },
    "glycolysis_glucose_utilization": {
        "input_terms": (
            "glucose",
            "glucose 6 phosphate",
            "fructose 6 phosphate",
            "lactate",
            "pyruvate",
            "glycolysis",
        ),
        "pool_ids": ("glycolysis_glucose_pool",),
        "pathway_terms": ("glycolysis", "glucose", "gluconeogenesis", "lactate", "pyruvate", "hexokinase"),
    },
    "tca_mitochondrial_pressure": {
        "input_terms": ("citrate", "isocitrate", "alpha-ketoglutarate", "succinate", "fumarate", "malate", "oxaloacetate", "tca"),
        "pool_ids": ("succinate_pool",),
        "pathway_terms": ("tca", "citric acid", "tricarboxylic", "succinate", "fumarate", "malate", "mitochondrial"),
    },
    "mitochondrial_tca_pressure": {
        "input_terms": ("citrate", "isocitrate", "alpha-ketoglutarate", "succinate", "fumarate", "malate", "oxaloacetate", "tca"),
        "pool_ids": ("succinate_pool",),
        "pathway_terms": ("tca", "citric acid", "tricarboxylic", "succinate", "fumarate", "malate", "mitochondrial"),
    },
    "osmolyte_epithelial_homeostasis": {
        "input_terms": ("taurine", "hypotaurine", "betaine", "myo inositol", "myo-inositol", "sorbitol", "glycerophosphocholine"),
        "pool_ids": ("osmolyte_epithelial_pool",),
        "pathway_terms": ("osmolyte", "osmoregulation", "taurine", "inositol", "epithelial", "transport"),
    },
    "epithelial_barrier_membrane_repair": {
        "input_terms": (
            "udp n acetylglucosamine",
            "udp-n-acetylglucosamine",
            "glucosamine",
            "choline",
            "phosphocholine",
            "glycerophosphocholine",
            "ceramide",
            "sphingosine",
            "proline",
            "hydroxyproline",
        ),
        "pool_ids": ("epithelial_barrier_repair_pool", "udp_sugar_pool", "phospholipid_membrane_pool", "sphingolipid_membrane_pool"),
        "pathway_terms": ("epithelial", "barrier", "tight junction", "membrane repair", "glycosylation", "phospholipid"),
    },
    "bile_acid": {
        "input_terms": (
            "taurocholic",
            "glycocholic",
            "glycochenodeoxycholic",
            "chenodeoxycholic",
            "deoxycholic",
            "lithocholic",
            "ursodeoxycholic",
            "cholic acid",
            "bile acid",
            "bile salt",
        ),
        "pathway_terms": ("bile acid", "bile salt", "enterohepatic", "biliary transport"),
    },
    "tryptophan_kynurenine": {
        "input_terms": (
            "tryptophan",
            "kynurenine",
            "kynurenic",
            "xanthurenic",
            "anthranilic",
            "quinolinic",
            "indole",
            "nicotinamide",
        ),
        "pathway_terms": ("tryptophan", "kynurenine", "nad biosynthesis", "indole"),
    },
    "microbial_uremic": {
        "input_terms": (
            "indoxyl",
            "p-cresol",
            "cresol sulfate",
            "phenylacetylglutamine",
            "trimethylamine",
            "tmao",
        ),
        "pathway_terms": ("gut liver", "indole", "xenobiotic", "uremic", "trimethylamine"),
    },
    "sulfur_one_carbon": {
        "input_terms": (
            "cystathionine",
            "homocysteine",
            "methionine",
            "adenosylhomocysteine",
            "adenosylmethionine",
            "taurine",
            "hypotaurine",
            "methylthioadenosine",
        ),
        "pool_ids": ("cysteine_homocysteine_pool", "osmolyte_epithelial_pool"),
        "pathway_terms": ("methionine", "cysteine", "homocysteine", "sulfur amino", "taurine", "transsulfuration"),
    },
    "polyamine": {
        "input_terms": ("putrescine", "spermidine", "spermine", "acetylputrescine", "methylthioadenosine"),
        "pathway_terms": ("polyamine", "putrescine", "spermidine", "spermine"),
    },
    "sphingolipid": {
        "input_terms": ("sphingosine", "sphinganine", "ceramide", "sphingolipid"),
        "pool_ids": ("sphingolipid_membrane_pool",),
        "pathway_terms": ("sphingolipid", "ceramide", "sphingosine", "glycosphingolipid"),
    },
    "heme_bilirubin": {
        "input_terms": ("bilirubin", "biliverdin", "heme"),
        "pathway_terms": ("heme", "bilirubin", "porphyria"),
    },
    "carnitine_fatty_acid": {
        "input_terms": ("carnitine", "propionylcarnitine", "palmitoylcarnitine"),
        "pathway_terms": ("carnitine", "fatty acid", "beta-oxidation"),
    },
    "renal_organic_anion_uremic_toxin_handling": {
        "input_terms": (
            "indoxyl sulfate",
            "p cresol sulfate",
            "p-cresol sulfate",
            "hippuric acid",
            "hippurate",
            "phenylacetylglutamine",
            "tmao",
            "trimethylamine n oxide",
            "uric acid",
            "urate",
            "creatinine",
        ),
        "pool_ids": ("renal_organic_anion_uremic_toxin_pool",),
        "pathway_terms": ("organic anion", "uremic", "indoxyl", "urate", "creatinine", "renal transport", "slc22", "oat"),
    },
    "renal_osmolyte_homeostasis": {
        "input_terms": ("taurine", "hypotaurine", "betaine", "myo inositol", "myo-inositol", "sorbitol", "glycerophosphocholine"),
        "pool_ids": ("osmolyte_epithelial_pool",),
        "pathway_terms": ("osmolyte", "osmoregulation", "taurine", "inositol", "renal", "kidney", "solute transport"),
    },
    "acylcarnitine_fatty_acid_oxidation_pressure": {
        "input_terms": (
            "carnitine",
            "acetylcarnitine",
            "propionylcarnitine",
            "butyrylcarnitine",
            "palmitoylcarnitine",
            "stearoylcarnitine",
            "acylcarnitine",
        ),
        "pool_ids": ("acylcarnitine_fatty_acid_oxidation_pool",),
        "pathway_terms": ("carnitine", "acylcarnitine", "fatty acid oxidation", "beta oxidation", "mitochondrial"),
    },
    "arginine_no_metabolism": {
        "input_terms": ("arginine", "citrulline", "ornithine", "adma", "asymmetric dimethylarginine", "nitrite", "nitrate", "nitric oxide"),
        "pool_ids": ("arginine_no_pool",),
        "pathway_terms": ("arginine", "nitric oxide", "nitrite", "urea cycle", "citrulline", "nos", "no synthesis"),
    },
    "purine_degradation_nucleotide_stress": {
        "input_terms": ("hypoxanthine", "xanthine", "uric acid", "urate", "adenosine", "amp", "adenosine monophosphate"),
        "pool_ids": ("purine_degradation_pool",),
        "pathway_terms": ("purine", "nucleotide", "adenosine", "xanthine", "urate", "uric acid"),
    },
    "idh_like_metabolic_pressure": {
        "input_terms": (
            "2 hydroxyglutarate",
            "2-hydroxyglutarate",
            "hydroxyglutarate",
            "alpha ketoglutarate",
            "alpha-ketoglutarate",
            "akg",
            "2 oxoglutarate",
            "2-oxoglutarate",
        ),
        "pool_ids": ("idh_like_pool",),
        "pathway_terms": ("idh", "isocitrate dehydrogenase", "2 hydroxyglutarate", "2-hydroxyglutarate", "alpha ketoglutarate"),
    },
    "glial_glutamate_gaba_metabolism": {
        "input_terms": ("glutamate", "glutamine", "gaba", "gamma aminobutyric", "n acetylaspartate", "n-acetylaspartate", "myo inositol", "myo-inositol"),
        "pool_ids": ("glial_glutamate_gaba_pool", "glutamate_glutamine_pool", "osmolyte_epithelial_pool"),
        "pathway_terms": ("glutamate", "glutamine", "gaba", "neurotransmitter", "astrocy", "glial", "myo-inositol"),
    },
}

THEME_ID_ALIASES = {
    "sulfur_one_carbon": "sulfur_one_carbon_metabolism",
    "sphingolipid": "phospholipid_membrane_remodeling",
    "carnitine_fatty_acid": "acylcarnitine_fatty_acid_oxidation_pressure",
    "microbial_uremic": "renal_organic_anion_uremic_toxin_handling",
    "tca_mitochondrial_pressure": "mitochondrial_tca_pressure",
}

CONTEXT_THEME_PRIORS = {
    "renal_organic_anion_uremic_toxin_handling": {
        "renal_organic_anion_uremic_toxin_handling": 0.16,
        "renal_osmolyte_homeostasis": 0.08,
        "acylcarnitine_fatty_acid_oxidation_pressure": 0.08,
        "arginine_no_metabolism": 0.06,
        "purine_degradation_nucleotide_stress": 0.06,
    },
    "renal_osmolyte_homeostasis": {"renal_osmolyte_homeostasis": 0.14},
    "fatty_acid_oxidation_pressure": {"acylcarnitine_fatty_acid_oxidation_pressure": 0.12},
    "proximal_tubule_solute_transport_pressure": {
        "renal_organic_anion_uremic_toxin_handling": 0.12,
        "renal_osmolyte_homeostasis": 0.08,
    },
    "glial_glutamate_gaba_metabolism": {"glial_glutamate_gaba_metabolism": 0.16},
    "idh_like_metabolic_pressure": {"idh_like_metabolic_pressure": 0.14},
    "hypoxia_glycolytic_shift": {"glycolysis_glucose_utilization": 0.08, "hypoxia_glycolytic_shift": 0.1},
    "epithelial_barrier_membrane_repair": {
        "epithelial_barrier_membrane_repair": 0.14,
        "nucleotide_sugar_glycosylation": 0.08,
        "phospholipid_membrane_remodeling": 0.08,
        "renal_osmolyte_homeostasis": 0.04,
    },
    "nucleotide_sugar_glycosylation": {"nucleotide_sugar_glycosylation": 0.1},
    "phospholipid_membrane_remodeling": {"phospholipid_membrane_remodeling": 0.1},
}

CONTEXT_MISMATCH_GROUPS = {
    "astrocytic_glial": ("astrocy", "glial", "oligodendrocyte", "microglia"),
    "neuronal": ("neuronal", "neuron", "synaptic", "axon", "dendrite", "neurotransmitter"),
    "taste_perception": ("taste", "gustatory", "olfactory", "smell perception"),
    "immune_cell_specific": ("immune cell", "lymphocyte", "t cell", "b cell", "macrophage", "monocyte", "neutrophil", "dendritic cell"),
    "hematopoietic": ("hematopoietic", "erythrocyte", "platelet", "myeloid", "lymphoid"),
}

CANCER_CONTEXT_GROUPS = {
    "colorectal": ("colorectal", "colon cancer", "rectal cancer", "crc", "coad", "read"),
    "melanoma": ("melanoma",),
    "breast": ("breast cancer", "brca", "tnbc", "triple negative breast"),
    "lung": ("lung cancer", "lung adenocarcinoma", "luad", "lusc"),
    "liver": ("hepatocellular carcinoma", "hcc", "liver cancer"),
    "skin_squamous": ("cutaneous squamous cell carcinoma", "cscc", "squamous cell carcinoma"),
}

EPITHELIAL_CONTEXT_TERMS = ("epithelial", "epithelium", "barrier", "mucosal", "intestinal", "colon", "airway")

DRUG_OVERLAY_DOWNGRADE_TERMS = (
    "sulfuric acid",
    "generic class",
    "drug class",
    "chemical class",
    "non therapeutic",
    "non-therapeutic",
    "reagent",
    "vehicle",
    "solvent",
    "chemical entity",
    "pan inhibitor",
    "broad inhibitor",
)

GENERIC_PATHWAY_TERMS = (
    "transcription",
    "translation",
    "rrna",
    "ribosome",
    "dna damage",
    "dna repair",
    "nuclear envelope",
    "signaling",
    "smad",
    "foxo",
    "pten",
    "axin",
)

BIOCHEMICAL_TOKEN_STOPWORDS = {
    "acid",
    "acids",
    "alpha",
    "beta",
    "delta",
    "gamma",
    "d",
    "l",
    "n",
    "o",
    "s",
    "r",
    "cis",
    "trans",
    "and",
    "or",
    "the",
    "with",
    "from",
    "into",
    "via",
    "type",
    "family",
    "member",
    "protein",
    "proteins",
    "gene",
    "genes",
    "regulation",
    "regulates",
    "signaling",
    "pathway",
    "pathways",
    "metabolism",
    "metabolic",
    "biosynthesis",
    "degradation",
    "synthesis",
    "transport",
    "disease",
    "disorder",
    "defective",
    "deficiency",
}

NODE_COLUMNS = [
    "node_uid",
    "node_idx",
    "node_type",
    "canonical_name",
    "display_name",
    "primary_external_id",
    "external_xrefs",
    "source_priority",
    "source_release",
    "license_id",
]

EDGE_COLUMNS = [
    "edge_uid",
    "edge_type",
    "subject_uid",
    "subject_type",
    "predicate",
    "object_uid",
    "object_type",
    "weight",
    "source_table",
    "source_name",
    "source_record_id",
    "evidence_level",
    "source_release",
    "license_id",
    "parser_hash",
    "metadata_json",
]

ALLOWED_EXPLANATION_EDGE_TYPES = {
    "metabolite_same_as_metabolite",
    "metabolite_participates_in_pathway",
    "metabolite_participates_in_reaction",
    "metabolite_regulates_gene",
    "metabolite_changed_in_cancer",
    "metabolite_associated_with_disease",
    "protein_participates_in_reaction",
    "gene_participates_in_reaction",
    "reaction_in_pathway",
    "gene_involved_in_pathway",
    "gene_regulates_metabolic_process",
    "target_maps_to_gene",
    "target_associated_with_disease",
    "pathway_associated_with_disease",
    "pathway_parent_of_pathway",
    "disease_is_a_disease",
}

TARGET_ANALYSIS_TYPES = {"pathway", "reaction", "gene", "protein", "target", "disease"}

PREDICTION_TYPE_WEIGHTS = {
    "metabolic": 1.0,
    "transport": 1.0,
    "biochemical_process": 0.8,
    "signaling": 0.5,
    "regulatory": 0.5,
    "immune_inflammation": 0.4,
    "development_cell_cycle": 0.4,
    "transcription_translation": 0.4,
    "disease": 0.3,
    "target": 0.3,
    "drug": 0.2,
    "drug_metabolism": 0.2,
    "literature_hotspot": 0.1,
    "cell_context_overlay": 0.2,
    "unknown_mixed": 0.5,
}

PREDICTION_DISTANCE_WEIGHTS = {0: 1.0, 1: 0.85, 2: 0.55}

CONFIDENCE_TIERS = (
    ("high", 0.75),
    ("medium", 0.45),
    ("exploratory", 0.18),
    ("low", 0.0),
)

METABOLIC_STATE_SIGNATURES = {
    "hypoxia_glycolytic_shift": {
        "label": "Hypoxia or glycolytic shift",
        "terms": ("lactate", "pyruvate", "glucose", "glycolysis", "monocarboxylate"),
        "theme_ids": ("glycolysis_glucose_utilization",),
    },
    "mitochondrial_pressure": {
        "label": "Mitochondrial or TCA-cycle pressure",
        "terms": ("citrate", "succinate", "fumarate", "malate", "tca", "oxidative phosphorylation"),
        "theme_ids": ("mitochondrial_tca_pressure",),
    },
    "oxidative_stress": {
        "label": "Oxidative stress or redox imbalance",
        "terms": ("glutathione", "cysteine", "methionine", "nad", "redox", "oxidative"),
        "theme_ids": ("glutathione_redox_stress", "amino_acid_redox_stress"),
    },
    "one_carbon_pressure": {
        "label": "One-carbon or methylation pressure",
        "terms": ("serine", "glycine", "folate", "methionine", "homocysteine", "one carbon"),
        "theme_ids": ("sulfur_one_carbon_metabolism",),
    },
    "lipid_remodeling": {
        "label": "Lipid remodeling",
        "terms": ("fatty acid", "carnitine", "ceramide", "sphingosine", "cholesterol", "lipid"),
        "theme_ids": ("phospholipid_membrane_remodeling",),
    },
    "microbial_or_xenobiotic_signal": {
        "label": "Microbial co-metabolism or xenobiotic signal",
        "terms": ("indole", "cresol", "tmao", "xenobiotic", "drug", "bile acid"),
    },
}

PATHWAY_TYPE_KEYWORDS = {
    "drug_metabolism": (
        "drug",
        "xenobiotic",
        "cytochrome p450",
        "phase i",
        "phase ii",
        "adme",
        "detoxification",
    ),
    "disease": (
        "disease",
        "cancer",
        "carcinoma",
        "tumor",
        "neoplasm",
        "syndrome",
        "defective",
        "deficiency",
        "disorder",
        "infection",
        "viral",
    ),
    "transport": (
        "transport",
        "transporter",
        "transmembrane",
        "export",
        "import",
        "uptake",
        "secretion",
        "carrier",
        "slc",
        "abc transporter",
    ),
    "metabolic": (
        "metabolism",
        "metabolic",
        "biosynthesis",
        "catabolism",
        "degradation",
        "glycolysis",
        "gluconeogenesis",
        "tca",
        "citric acid",
        "fatty acid",
        "lipid",
        "amino acid",
        "nucleotide",
        "bile acid",
        "bile salt",
        "one carbon",
        "folate",
        "urea cycle",
        "oxidation",
        "beta oxidation",
        "glutathione",
        "sphingolipid",
        "ceramide",
        "polyamine",
        "porphyrin",
        "heme",
    ),
    "immune_inflammation": (
        "immune",
        "inflammation",
        "inflammatory",
        "interleukin",
        "cytokine",
        "interferon",
        "t cell",
        "b cell",
        "complement",
    ),
    "transcription_translation": (
        "transcription",
        "translation",
        "ribosome",
        "rrna",
        "mrna",
        "splicing",
        "chromatin",
        "histone",
    ),
    "development_cell_cycle": (
        "cell cycle",
        "mitosis",
        "meiosis",
        "development",
        "differentiation",
        "apoptosis",
        "senescence",
    ),
    "signaling": (
        "signaling",
        "signal",
        "regulation",
        "regulates",
        "cascade",
        "receptor",
        "kinase",
        "phosphorylation",
        "pten",
        "mapk",
        "pi3k",
        "akt",
        "mtor",
        "wnt",
        "notch",
        "tgf",
        "smad",
    ),
}

CONTEXT_FIELD_ALIASES = {
    "cancer": "cancer_type",
    "cancer_type": "cancer_type",
    "disease": "cancer_type",
    "phenotype_context": "cancer_type",
    "tumor_type": "cancer_type",
    "oncotree": "cancer_type",
    "tissue": "tissue",
    "organ": "tissue",
    "site": "tissue",
    "body_site": "tissue",
    "cell_type": "cell_type",
    "celltype": "cell_type",
    "cell": "cell_type",
    "cell_state": "cell_state",
    "cellstate": "cell_state",
    "state": "cell_state",
    "phenotype": "cell_state",
    "cell_line": "cell_line",
    "cellline": "cell_line",
    "model": "cell_line",
    "context": "context_terms",
    "context_terms": "context_terms",
    "condition": "context_terms",
    "experimental_context": "experimental_context",
    "experiment": "experimental_context",
    "species": "species",
    "organism": "species",
    "platform": "platform",
    "sample_type": "sample_type",
    "sample": "sample_type",
}

PREDICTION_OVERLAY_FILES = {
    "drug_targets": ("drug_targets", "drug_target_edges", "drug_predictions"),
    "cell_type_signatures": ("cell_type_signatures", "cell_signatures", "cell_context_signatures"),
}

PROPAGATION_ALLOWED_TRANSITIONS = {
    ("metabolite", "metabolite_same_as_metabolite", "metabolite"),
    ("metabolite", "metabolite_participates_in_pathway", "pathway"),
    ("metabolite", "metabolite_participates_in_reaction", "reaction"),
    ("metabolite", "metabolite_regulates_gene", "gene"),
    ("metabolite", "metabolite_regulates_gene", "target"),
    ("metabolite", "metabolite_changed_in_cancer", "disease"),
    ("metabolite", "metabolite_associated_with_disease", "disease"),
    ("reaction", "reaction_in_pathway", "pathway"),
    ("reaction", "gene_participates_in_reaction", "gene"),
    ("reaction", "protein_participates_in_reaction", "protein"),
    ("pathway", "gene_involved_in_pathway", "gene"),
    ("pathway", "gene_regulates_metabolic_process", "gene"),
    ("pathway", "gene_regulates_metabolic_process", "target"),
    ("pathway", "pathway_associated_with_disease", "disease"),
    ("pathway", "pathway_parent_of_pathway", "pathway"),
    ("gene", "target_maps_to_gene", "target"),
    ("gene", "gene_involved_in_pathway", "pathway"),
    ("gene", "gene_regulates_metabolic_process", "pathway"),
    ("gene", "gene_participates_in_reaction", "reaction"),
    ("protein", "protein_participates_in_reaction", "reaction"),
    ("target", "metabolite_regulates_gene", "metabolite"),
    ("target", "gene_regulates_metabolic_process", "pathway"),
    ("target", "target_associated_with_disease", "disease"),
    ("target", "target_maps_to_gene", "gene"),
    ("disease", "disease_is_a_disease", "disease"),
}

PROPAGATION_EDGE_PRIORITY = {
    ("metabolite", "metabolite_same_as_metabolite", "metabolite"): 0,
    ("metabolite", "metabolite_participates_in_pathway", "pathway"): 0,
    ("metabolite", "metabolite_participates_in_reaction", "reaction"): 1,
    ("metabolite", "metabolite_regulates_gene", "gene"): 2,
    ("metabolite", "metabolite_regulates_gene", "target"): 2,
    ("metabolite", "metabolite_changed_in_cancer", "disease"): 3,
    ("metabolite", "metabolite_associated_with_disease", "disease"): 3,
    ("reaction", "reaction_in_pathway", "pathway"): 0,
    ("reaction", "gene_participates_in_reaction", "gene"): 1,
    ("reaction", "protein_participates_in_reaction", "protein"): 1,
    ("pathway", "gene_involved_in_pathway", "gene"): 0,
    ("pathway", "gene_regulates_metabolic_process", "gene"): 1,
    ("pathway", "gene_regulates_metabolic_process", "target"): 1,
    ("pathway", "pathway_associated_with_disease", "disease"): 2,
    ("pathway", "pathway_parent_of_pathway", "pathway"): 1,
    ("gene", "target_maps_to_gene", "target"): 0,
    ("gene", "gene_involved_in_pathway", "pathway"): 1,
    ("gene", "gene_regulates_metabolic_process", "pathway"): 1,
    ("target", "gene_regulates_metabolic_process", "pathway"): 1,
    ("target", "target_associated_with_disease", "disease"): 0,
    ("target", "target_maps_to_gene", "gene"): 1,
    ("target", "metabolite_regulates_gene", "metabolite"): 2,
    ("disease", "disease_is_a_disease", "disease"): 0,
}

GRAPH_EDGE_TO_LITERATURE_PREDICATE = {
    "target_associated_with_disease": "target_associated_with_disease",
    "gene_involved_in_pathway": "gene_regulates_metabolic_process",
}


def propagation_allows(current_type: str, edge_type: str, neighbor_type: str) -> bool:
    return (current_type, edge_type, neighbor_type) in PROPAGATION_ALLOWED_TRANSITIONS


def propagation_priority(current_type: str, edge_type: str, neighbor_type: str) -> int:
    return PROPAGATION_EDGE_PRIORITY.get((current_type, edge_type, neighbor_type), 99)


def require_arrow() -> None:
    if pa is None or ds is None or pq is None:
        raise RuntimeError("pyarrow is required for the metabolism service layer.")


def stable_json(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def content_hash(value: Any) -> str:
    return hashlib.sha256(stable_json(value).encode("utf-8")).hexdigest()


def clamp_probability(value: Any, default: float = 0.0) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        number = default
    if math.isnan(number) or math.isinf(number):
        number = default
    return min(0.999999, max(0.000001, number))


def read_json_file(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def latest_release_id(root: Path) -> str:
    return latest_normalized_release_id(root)


def table_dataset(path: Path) -> ds.Dataset:
    require_arrow()
    if not path.exists():
        raise FileNotFoundError(f"Missing parquet table: {path}")
    return ds.dataset(path, format="parquet")


def safe_metadata(metadata_json: str | None) -> dict[str, Any]:
    if not metadata_json:
        return {}
    try:
        parsed = json.loads(metadata_json)
    except (TypeError, ValueError):
        return {}
    return parsed if isinstance(parsed, dict) else {}


def parse_list(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, list):
        return [str(item) for item in value if item is not None and str(item) != ""]
    if isinstance(value, tuple):
        return [str(item) for item in value if item is not None and str(item) != ""]
    return [str(value)] if str(value) else []


def normalize_prediction_text(value: Any) -> str:
    text = str(value or "").casefold()
    text = re.sub(r"[^0-9a-z]+", " ", text).strip()
    return re.sub(r"\s+", " ", text)


def prediction_terms_from_value(value: Any) -> list[str]:
    terms: set[str] = set()

    def add(raw: Any) -> None:
        if raw is None:
            return
        if isinstance(raw, (list, tuple, set)):
            for item in raw:
                add(item)
            return
        if isinstance(raw, dict):
            for item in raw.values():
                add(item)
            return
        text = str(raw).strip()
        if not text:
            return
        parts = re.split(r"[;|,\n]+", text)
        for part in parts:
            normalized = normalize_prediction_text(part)
            if not normalized:
                continue
            terms.add(normalized)
            for token in normalized.split():
                if len(token) >= 3:
                    terms.add(token)

    add(value)
    return sorted(terms)


def gold_term_variants(term: str) -> set[str]:
    normalized = normalize_prediction_text(term)
    lookup_normalized = normalize_lookup_key(term)
    variants = {value for value in (normalized, lookup_normalized) if value}
    for canonical, aliases in COMMON_METABOLITE_NAME_ALIASES.items():
        canonical_norm = normalize_prediction_text(canonical)
        canonical_lookup = normalize_lookup_key(canonical)
        alias_norms = {normalize_prediction_text(alias) for alias in aliases}
        alias_lookup_norms = {normalize_lookup_key(alias) for alias in aliases}
        if normalized == canonical_norm or lookup_normalized == canonical_lookup:
            variants.update(alias_norms)
            variants.update(alias_lookup_norms)
        elif normalized in alias_norms or lookup_normalized in alias_lookup_norms:
            variants.add(canonical_norm)
            variants.add(canonical_lookup)
            variants.update(alias_norms)
            variants.update(alias_lookup_norms)
    return {variant for variant in variants if variant}


def normalized_text_has_term(normalized_text: str, term: str) -> bool:
    if not normalized_text or not term:
        return False
    for variant in gold_term_variants(term):
        if " " in variant:
            if variant in normalized_text:
                return True
        elif re.search(r"[^a-z0-9]", variant):
            pattern = rf"(?<![a-z0-9]){re.escape(variant)}(?![a-z0-9])"
            if re.search(pattern, normalized_text):
                return True
        elif variant in set(re.findall(r"[a-z][a-z0-9]+", normalized_text)):
            return True
    return False


def biochemical_rescue_key(value: Any) -> str:
    text = normalize_lookup_key(value)
    if not text:
        return ""
    text = text.replace("+", " plus ")
    text = re.sub(r"\((?:\d*[+-]|[+-])\)", " ", text)
    text = re.sub(r"\b(?:alpha|beta|gamma|delta|d|l|r|s|o|n)\b", " ", text)
    text = text.replace("5 diphospho", "diphospho").replace("5 monophosphate", "monophosphate")
    text = re.sub(r"[^a-z0-9+]+", " ", text)
    text = re.sub(r"\s+", " ", text).strip()
    replacements = {
        "glucose 6 phosphate": "glucose 6 phosphate",
        "fructose 6 phosphate": "fructose 6 phosphate",
        "uridine diphosphoglucose": "udp glucose",
        "udp glucose": "udp glucose",
        "udp galactose": "udp galactose",
        "udp n acetyl glucosamine": "udp n acetylglucosamine",
        "udp n acetylglucosamine": "udp n acetylglucosamine",
        "glucopyranose 6 phosphate": "glucose 6 phosphate",
        "fructose 1 6 diphosphate": "fructose 1 6 bisphosphate",
        "isocitric acid": "isocitrate",
        "phosphonatoenolpyruvate": "phosphoenolpyruvate",
        "glutathione disulfide": "oxidized glutathione",
        "gssg": "oxidized glutathione",
        "guanosine 5 monophosphate": "gmp",
        "guanosine monophosphate": "gmp",
        "guanosine 5 diphosphate": "gdp",
        "guanosine diphosphate": "gdp",
        "guanosine 5 triphosphate": "gtp",
        "guanosine triphosphate": "gtp",
        "cytidine 5 monophosphate": "cmp",
        "cytidine monophosphate": "cmp",
        "cytidine 5 diphosphate": "cdp",
        "cytidine diphosphate": "cdp",
        "cytidine 5 triphosphate": "ctp",
        "cytidine triphosphate": "ctp",
        "uridine 5 monophosphate": "ump",
        "uridine monophosphate": "ump",
        "uridine 5 diphosphate": "udp",
        "uridine diphosphate": "udp",
        "uridine 5 triphosphate": "utp",
        "uridine triphosphate": "utp",
        "adenosine 5 monophosphate": "amp",
        "adenosine monophosphate": "amp",
        "adenosine 5 diphosphate": "adp",
        "adenosine diphosphate": "adp",
        "adenosine 5 triphosphate": "atp",
        "adenosine triphosphate": "atp",
        "5 phosphono ribofuranosyl diphosphate": "prpp",
        "phosphono ribofuranosyl diphosphate": "prpp",
        "ribofuranosyl diphosphate": "prpp",
        "nad plus": "nad",
        "nicotinamide adenine dinucleotide": "nad",
        "o acetyl carnitine": "acetylcarnitine",
        "acetyl carnitine": "acetylcarnitine",
        "o palmitoyl carnitine": "palmitoylcarnitine",
        "palmitoyl carnitine": "palmitoylcarnitine",
        "methionine oxide": "methionine sulfoxide",
    }
    for source, target in replacements.items():
        if source in text:
            text = text.replace(source, target)
    return re.sub(r"\s+", " ", text).strip()


def biochemical_rescue_compatible(query: Any, candidate_name: Any) -> bool:
    query_key = biochemical_rescue_key(query)
    candidate_key = biochemical_rescue_key(candidate_name)
    if not query_key or not candidate_key:
        return False
    if query_key == candidate_key:
        return True
    compact_query = query_key.replace(" ", "")
    compact_candidate = candidate_key.replace(" ", "")
    if compact_query and compact_query == compact_candidate:
        return True
    if len(compact_query) >= 6 and compact_query in compact_candidate:
        return True
    aliases = COMMON_METABOLITE_NAME_ALIASES.get(normalize_lookup_key(query), [])
    return any(biochemical_rescue_key(alias) == candidate_key for alias in aliases)


def is_common_biochemical_rescue_query(query: Any) -> bool:
    normalized = normalize_lookup_key(query)
    rescue_key = biochemical_rescue_key(query)
    if not normalized or not rescue_key:
        return False
    compact = rescue_key.replace(" ", "")
    tokens = set(re.findall(r"[a-z][a-z0-9+]*", rescue_key))
    if normalized in COMMON_METABOLITE_NAME_ALIASES:
        return True
    if compact in {term.replace(" ", "").replace("-", "") for term in COMMON_BIOCHEMICAL_ROOT_TERMS}:
        return True
    if tokens & COMMON_BIOCHEMICAL_ROOT_TERMS:
        return True
    if any(token.endswith("carnitine") for token in tokens):
        return True
    if any(token.startswith("udp") or token in {"adp", "atp", "amp", "gdp", "gtp", "gmp", "udp", "utp", "ump", "ctp", "cmp"} for token in tokens):
        return True
    if any(token.startswith("nad") or token in {"fadh", "fad", "fmnh", "fmn"} for token in tokens):
        return True
    if any(term in rescue_key for term in COMMON_BIOCHEMICAL_MODIFIER_TERMS) and (
        tokens & COMMON_BIOCHEMICAL_ROOT_TERMS or len(tokens) <= 4
    ):
        return True
    lipid_headgroups = ("phosphatidyl", "lysophosphatidyl", "phosphocholine", "phosphoethanolamine", "sphingosine", "ceramide")
    if any(term in normalized for term in lipid_headgroups):
        return True
    return False


def common_biochemical_zero_score_rescue_rank(query: Any, candidate: dict[str, Any]) -> float:
    query_key = biochemical_rescue_key(query)
    display_key = biochemical_rescue_key(candidate.get("display_name", ""))
    if not query_key or not display_key or not biochemical_rescue_compatible(query, candidate.get("display_name", "")):
        return -1.0
    matches = candidate.get("matches") or []
    exact_blocked_matches = [
        match
        for match in matches
        if match.get("component") == "name"
        and match.get("triage_policy") == "auto_abstain_name_only"
        and biochemical_rescue_compatible(query, match.get("raw_value", ""))
    ]
    if not exact_blocked_matches:
        return -1.0
    display_text = display_key
    rank = 10.0
    compact_query = query_key.replace(" ", "")
    compact_display = display_key.replace(" ", "")
    if compact_query == compact_display:
        rank += 8.0
    aliases = [biochemical_rescue_key(alias) for alias in COMMON_METABOLITE_NAME_ALIASES.get(normalize_lookup_key(query), [])]
    if display_key in aliases:
        rank += 6.0
    if any(match.get("match_field") == "canonical_name" for match in exact_blocked_matches):
        rank += 4.0
    if any(match.get("match_field") in {"pubchem_title", "synonym", "pubchem_synonym"} for match in exact_blocked_matches):
        rank += 2.0
    if any(term in display_text for term in (" residue", "zwitterion", "conjugate")):
        rank -= 6.0
    if display_text.endswith("ate") and not query_key.endswith("ate"):
        rank -= 4.0
    return rank


def exact_prediction_terms_from_value(value: Any) -> list[str]:
    terms: set[str] = set()

    def add(raw: Any) -> None:
        if raw is None:
            return
        if isinstance(raw, (list, tuple, set)):
            for item in raw:
                add(item)
            return
        if isinstance(raw, dict):
            for item in raw.values():
                add(item)
            return
        text = str(raw).strip()
        if not text:
            return
        for part in re.split(r"[;|,\n]+", text):
            normalized = normalize_prediction_text(part)
            if normalized:
                terms.add(normalized)
            if ":" in part:
                suffix = normalize_prediction_text(part.split(":", 1)[1])
                if suffix:
                    terms.add(suffix)

    add(value)
    return sorted(terms)


def row_value(row: dict[str, Any], *keys: str) -> Any:
    lowered = {str(key).casefold(): value for key, value in row.items()}

    def present(value: Any) -> bool:
        return value is not None and value != ""

    for key in keys:
        if key in row and present(row[key]):
            return row[key]
        value = lowered.get(key.casefold())
        if present(value):
            return value
    return ""


def split_overlay_terms(row: dict[str, Any], *keys: str) -> list[str]:
    terms: set[str] = set()
    for key in keys:
        terms.update(prediction_terms_from_value(row_value(row, key)))
    return sorted(terms)


def split_overlay_exact_terms(row: dict[str, Any], *keys: str) -> list[str]:
    terms: set[str] = set()
    for key in keys:
        terms.update(exact_prediction_terms_from_value(row_value(row, key)))
    return sorted(terms)


def load_records_file(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    suffix = path.suffix.casefold()
    inner_suffix = path.with_suffix("").suffix.casefold() if suffix == ".gz" else suffix
    opener = gzip.open if suffix == ".gz" else open
    if suffix in {".json", ".jsonl"}:
        if suffix == ".jsonl":
            rows = []
            for line in path.read_text(encoding="utf-8").splitlines():
                if line.strip():
                    parsed = json.loads(line)
                    if isinstance(parsed, dict):
                        rows.append(parsed)
            return rows
        parsed = json.loads(path.read_text(encoding="utf-8"))
        if isinstance(parsed, dict):
            parsed = parsed.get("records", [])
        return [row for row in parsed if isinstance(row, dict)] if isinstance(parsed, list) else []
    delimiter = "\t" if inner_suffix == ".tsv" else ","
    with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter=delimiter)]


def normalize_query_notes(notes: Any) -> list[dict[str, Any]]:
    if notes is None or notes == "":
        return []
    raw_notes = notes if isinstance(notes, list) else [notes]
    normalized: list[dict[str, Any]] = []
    for note in raw_notes:
        if isinstance(note, dict):
            normalized.append(note)
        else:
            normalized.append({"mode": str(note), "active": False, "raw_note": str(note)})
    return normalized


def table_rows(path: Path, columns: list[str] | None = None, filter_expr: ds.Expression | None = None) -> list[dict[str, Any]]:
    dataset = table_dataset(path)
    table = dataset.to_table(columns=columns, filter=filter_expr)
    return table.to_pylist()


def score_resolver_row(row: dict[str, Any]) -> tuple[float, dict[str, float]]:
    lookup_rank = float(row.get("rank") or 0.0)
    query_priority = float(row.get("_query_priority") or 0.0)
    priority_score = min(100.0, max(0.0, query_priority / 120.0 * 100.0))
    score = min(100.0, (0.65 * lookup_rank) + (0.35 * priority_score))
    components = {
        "lookup_rank": round(lookup_rank, 6),
        "query_priority_score": round(priority_score, 6),
        "feature_support_score": 0.0,
        "ambiguity_penalty": 0.0,
    }
    return round(score, 6), components


def edge_probability(edge: dict[str, Any]) -> tuple[float, dict[str, float]]:
    edge_type = str(edge.get("edge_type") or "")
    evidence = str(edge.get("evidence_level") or "").casefold()
    source = str(edge.get("source_name") or "").casefold()
    source_table = str(edge.get("source_table") or "").casefold()
    metadata = safe_metadata(edge.get("metadata_json"))

    p_curated = 0.0
    support_rows = edge.get("literature_support") or []
    if support_rows:
        p_literature = independent_probability_union([float(row.get("p_literature") or 0.0) for row in support_rows])
    else:
        p_literature = clamp_probability(edge.get("p_literature"), 0.0)
    p_topology = 0.0
    p_user = 0.0

    if "literature_overlay" in evidence:
        p_topology = 0.0
    elif edge_type == "target_associated_with_disease":
        p_topology = clamp_probability(edge.get("weight"), 0.0)
    elif "curated" in evidence or source in {"reactome", "bridgedb"}:
        p_curated = 0.95
    elif edge_type in {"target_maps_to_gene", "pathway_parent_of_pathway", "disease_is_a_disease"}:
        p_curated = 0.85
    elif source_table in {"reactions", "reaction_participants"}:
        p_curated = 0.9
    else:
        p_topology = clamp_probability(edge.get("weight"), 0.5)

    if edge_type == "target_associated_with_disease":
        try:
            score_json = json.loads(metadata.get("score_components_json") or "{}")
        except (TypeError, ValueError):
            score_json = {}
        if isinstance(score_json, dict) and "associationScore" in score_json:
            p_topology = clamp_probability(score_json.get("associationScore"), p_topology)

    p_final = 1.0
    for component in (p_curated, p_literature, p_topology, p_user):
        p_final *= 1.0 - clamp_probability(component, 0.0)
    p_final = clamp_probability(1.0 - p_final, 0.0)
    return p_final, {
        "p_curated": round(p_curated, 6),
        "p_literature": round(p_literature, 6),
        "p_topology": round(p_topology, 6),
        "p_user": round(p_user, 6),
        "p_final": round(p_final, 6),
    }


def hypergeom_sf(x: int, population: int, success_population: int, draws: int) -> float:
    if x <= 0:
        return 1.0
    if population <= 0 or success_population <= 0 or draws <= 0:
        return 1.0
    max_k = min(success_population, draws)
    denominator = math.lgamma(population + 1) - math.lgamma(draws + 1) - math.lgamma(population - draws + 1)
    total = 0.0
    for k in range(x, max_k + 1):
        if draws - k > population - success_population:
            continue
        log_p = (
            math.lgamma(success_population + 1)
            - math.lgamma(k + 1)
            - math.lgamma(success_population - k + 1)
            + math.lgamma(population - success_population + 1)
            - math.lgamma(draws - k + 1)
            - math.lgamma(population - success_population - draws + k + 1)
            - denominator
        )
        total += math.exp(log_p)
    return min(1.0, max(0.0, total))


def benjamini_hochberg(rows: list[dict[str, Any]], p_value_key: str = "p_value") -> None:
    if not rows:
        return
    ordered = sorted(enumerate(rows), key=lambda item: (item[1].get(p_value_key, 1.0), item[0]))
    m = len(ordered)
    prev = 1.0
    for rank, (idx, row) in reversed(list(enumerate(ordered, start=1))):
        q = min(prev, float(row.get(p_value_key, 1.0)) * m / rank)
        rows[idx]["fdr"] = round(q, 12)
        prev = q


ADDUCTS = {
    "[M+H]+": {"delta": 1.007276466812, "charge": 1},
    "M+H": {"delta": 1.007276466812, "charge": 1},
    "[M+NA]+": {"delta": 22.989218, "charge": 1},
    "M+NA": {"delta": 22.989218, "charge": 1},
    "[M+K]+": {"delta": 38.963158, "charge": 1},
    "M+K": {"delta": 38.963158, "charge": 1},
    "[M+NH4]+": {"delta": 18.033823, "charge": 1},
    "M+NH4": {"delta": 18.033823, "charge": 1},
    "[M-H]-": {"delta": -1.007276466812, "charge": -1},
    "M-H": {"delta": -1.007276466812, "charge": -1},
    "[M+CL]-": {"delta": 34.969402, "charge": -1},
    "M+CL": {"delta": 34.969402, "charge": -1},
    "[M+HCOO]-": {"delta": 44.997655, "charge": -1},
    "M+HCOO": {"delta": 44.997655, "charge": -1},
    "[M+FA-H]-": {"delta": 44.997655, "charge": -1},
    "M+FA-H": {"delta": 44.997655, "charge": -1},
    "[M+2H]2+": {"delta": 2.014552933624, "charge": 2},
    "M+2H": {"delta": 2.014552933624, "charge": 2},
}


def normalize_adduct(value: Any) -> str:
    text = str(value or "").strip().upper().replace(" ", "")
    if not text:
        return ""
    text = text.replace("NA", "Na").replace("CL", "Cl")
    canonical = text.upper()
    return canonical


def maybe_number(value: Any) -> float | None:
    if value in {None, ""}:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if number == number and math.isfinite(number) else None


def clamp_unit(value: Any, default: float = 0.0) -> float:
    number = maybe_number(value)
    if number is None:
        number = default
    return min(1.0, max(0.0, number))


def independent_probability_union(values: list[float]) -> float:
    remaining = 1.0
    for value in values:
        remaining *= 1.0 - clamp_unit(value)
    return clamp_unit(1.0 - remaining)


def first_present(mapping: dict[str, Any], keys: tuple[str, ...]) -> Any:
    for key in keys:
        value = mapping.get(key)
        if value is not None and value != "":
            return value
    return None


def parse_log2_fold_change(record: dict[str, Any]) -> float | None:
    direct = first_present(record, ("log2fc", "log2_fc", "logfc", "log_fc", "log2_fold_change"))
    direct_number = maybe_number(direct)
    if direct_number is not None:
        return direct_number
    fold_change = maybe_number(first_present(record, ("fold_change", "foldchange", "fc")))
    if fold_change is None or fold_change <= 0:
        return None
    return math.log2(fold_change)


def p_value_support(value: Any) -> float:
    p_value = maybe_number(value)
    if p_value is None:
        return 0.0
    p_value = min(1.0, max(1e-300, p_value))
    return min(0.95, max(0.0, -math.log10(p_value) / 6.0))


def effect_support(log2_fold_change: float | None) -> float:
    if log2_fold_change is None:
        return 0.0
    return min(0.75, abs(log2_fold_change) / 4.0)


def direction_from_features(raw_direction: Any, log2_fold_change: float | None) -> str:
    text = str(raw_direction or "").strip().casefold()
    if text in {"up", "upregulated", "increase", "increased", "positive", "+", "1"}:
        return "up"
    if text in {"down", "downregulated", "decrease", "decreased", "negative", "-", "-1"}:
        return "down"
    if text in {"stable", "unchanged", "none", "0"}:
        return "unchanged"
    if log2_fold_change is None:
        return "unknown"
    if log2_fold_change > 0:
        return "up"
    if log2_fold_change < 0:
        return "down"
    return "unchanged"


def analysis_feature_from_record(record: Any, input_id: str = "") -> dict[str, Any]:
    normalized = {str(key).strip().casefold().replace(" ", "_"): value for key, value in record.items()} if isinstance(record, dict) else {}
    log2_fc = parse_log2_fold_change(normalized)
    p_value = maybe_number(first_present(normalized, ("p_value", "pvalue", "p", "p_ttest", "p_wilcoxon")))
    adjusted_p_value = maybe_number(first_present(normalized, ("adjusted_p_value", "padj", "fdr", "q_value", "qvalue", "q_ttest", "q_wilcoxon")))
    significance_value = adjusted_p_value if adjusted_p_value is not None else p_value
    direction = direction_from_features(first_present(normalized, ("direction", "regulation", "change_direction")), log2_fc)
    p_user = independent_probability_union([p_value_support(significance_value), effect_support(log2_fc)])
    significance_threshold = 0.05
    return {
        "input_id": input_id,
        "log2_fold_change": round(log2_fc, 6) if log2_fc is not None else None,
        "p_value": p_value,
        "adjusted_p_value": adjusted_p_value,
        "significance_value": significance_value,
        "significant": bool(significance_value is not None and significance_value <= significance_threshold),
        "direction": direction,
        "p_user": round(p_user, 6),
        "score_components": {
            "p_value_support": round(p_value_support(significance_value), 6),
            "effect_support": round(effect_support(log2_fc), 6),
            "p_user_formula": "1-(1-p_value_support)(1-effect_support)",
        },
    }


def normalize_method_key(value: Any) -> str:
    return normalize_lookup_key(str(value or "default"))


def parse_ms2_peaks(value: Any) -> list[tuple[float, float]]:
    peaks: list[tuple[float, float]] = []
    if value is None or value == "":
        return peaks
    if isinstance(value, list):
        for item in value:
            if isinstance(item, dict):
                mz = maybe_number(item.get("mz") or item.get("fragment_mz"))
                intensity = maybe_number(item.get("intensity") or item.get("fragment_intensity")) or 1.0
            elif isinstance(item, (list, tuple)) and item:
                mz = maybe_number(item[0])
                intensity = maybe_number(item[1]) if len(item) > 1 else 1.0
                intensity = intensity or 1.0
            else:
                mz = maybe_number(item)
                intensity = 1.0
            if mz and mz > 0:
                peaks.append((mz, max(0.0, intensity)))
        return peaks
    text = str(value).replace(";", " ").replace(",", " ")
    for token in text.split():
        if not token:
            continue
        if ":" in token:
            mz_text, intensity_text = token.split(":", 1)
        elif "|" in token:
            mz_text, intensity_text = token.split("|", 1)
        else:
            mz_text, intensity_text = token, "1"
        mz = maybe_number(mz_text)
        intensity = maybe_number(intensity_text) or 1.0
        if mz and mz > 0:
            peaks.append((mz, max(0.0, intensity)))
    return peaks


NUMERIC_CHAT_FIELDS = {
    "log2FC",
    "fold_change",
    "pvalue",
    "padj",
    "mz",
    "charge",
    "ppm_tolerance",
    "rt",
    "rt_tolerance",
    "ms2_mz_tolerance",
    "n_group1",
    "n_group2",
    "mean_group1",
    "mean_group2",
    "median_group1",
    "median_group2",
    "sd_group1",
    "sd_group2",
    "mean_diff",
    "median_diff",
    "pseudo_log2FC_shifted",
    "signed_log2FC_absmean",
    "cohen_d",
    "cliffs_delta",
    "t_welch",
    "z_welch",
    "u_stat",
    "z_wilcoxon",
}

SUPPORTED_CHAT_RECORD_FIELDS = {
    "input_id",
    "id",
    "name",
    "metabolite",
    "compound",
    "hmdb",
    "hmdb_id",
    "chebi",
    "chebi_id",
    "pubchem_cid",
    "cid",
    "pubchem",
    "kegg",
    "kegg_id",
    "inchikey",
    "inchi_key",
    "formula",
    "mz",
    "adduct",
    "charge",
    "ppm_tolerance",
    "rt",
    "rt_tolerance",
    "rt_method",
    "ms2_peaks",
    "msms_peaks",
    "fragments",
    "fragment_peaks",
    "ms2_mz_tolerance",
    "ion_mode",
    "polarity",
    "log2FC",
    "fold_change",
    "pvalue",
    "padj",
    "direction",
    "trait",
    "accession_id",
    "reported_trait",
    "summary_statistics_url",
    "pubmed_id",
    "paper_title",
    "journal",
    "publication_date",
    "efo_traits",
    "bg_traits",
    "initial_sample_description",
    "discovery_sample_ancestry",
    "cell_type",
    "celltype",
    "celltype_l",
    "group_col",
    "group1",
    "group2",
    "n_group1",
    "n_group2",
    "mean_group1",
    "mean_group2",
    "median_group1",
    "median_group2",
    "sd_group1",
    "sd_group2",
    "mean_diff",
    "median_diff",
    "pseudo_log2FC_shifted",
    "signed_log2FC_absmean",
    "cohen_d",
    "cliffs_delta",
}

CHAT_COLUMN_ALIASES = {
    "id": "input_id",
    "input": "input_id",
    "input_id": "input_id",
    "sample_id": "input_id",
    "row_id": "input_id",
    "name": "name",
    "names": "name",
    "metabolite": "name",
    "metabolite_name": "name",
    "compound": "name",
    "compound_name": "name",
    "feature_name": "name",
    "peak_name": "name",
    "代谢物": "name",
    "代谢物名称": "name",
    "名称": "name",
    "化合物": "name",
    "化合物名称": "name",
    "hmdb": "hmdb",
    "hmdb_id": "hmdb",
    "hmdbid": "hmdb",
    "chebi": "chebi",
    "chebi_id": "chebi",
    "chebiid": "chebi",
    "pubchem": "pubchem_cid",
    "pubchem_cid": "pubchem_cid",
    "pubchem_id": "pubchem_cid",
    "cid": "pubchem_cid",
    "kegg": "kegg",
    "kegg_id": "kegg",
    "keggid": "kegg",
    "inchikey": "inchikey",
    "inchi_key": "inchikey",
    "inchi_key_": "inchikey",
    "formula": "formula",
    "chemical_formula": "formula",
    "分子式": "formula",
    "mz": "mz",
    "m_z": "mz",
    "m_over_z": "mz",
    "mass_to_charge": "mz",
    "质荷比": "mz",
    "adduct": "adduct",
    "加合物": "adduct",
    "charge": "charge",
    "电荷": "charge",
    "ppm": "ppm_tolerance",
    "ppm_tolerance": "ppm_tolerance",
    "mass_error_ppm": "ppm_tolerance",
    "rt": "rt",
    "retention_time": "rt",
    "保留时间": "rt",
    "rt_tolerance": "rt_tolerance",
    "rt_method": "rt_method",
    "lc_method": "rt_method",
    "method": "rt_method",
    "ms2": "ms2_peaks",
    "ms2_peaks": "ms2_peaks",
    "msms": "ms2_peaks",
    "msms_peaks": "ms2_peaks",
    "fragments": "ms2_peaks",
    "fragment_peaks": "ms2_peaks",
    "ms2_mz_tolerance": "ms2_mz_tolerance",
    "fragment_mz_tolerance": "ms2_mz_tolerance",
    "ion_mode": "ion_mode",
    "polarity": "polarity",
    "log2fc": "log2FC",
    "log2_fc": "log2FC",
    "log2_fold_change": "log2FC",
    "logfc": "log2FC",
    "log_fc": "log2FC",
    "log2变化": "log2FC",
    "fold_change": "fold_change",
    "foldchange": "fold_change",
    "fc": "fold_change",
    "倍数变化": "fold_change",
    "p": "pvalue",
    "p_value": "pvalue",
    "pvalue": "pvalue",
    "pval": "pvalue",
    "p_val": "pvalue",
    "p_ttest": "pvalue",
    "p_wilcoxon": "pvalue",
    "p值": "pvalue",
    "padj": "padj",
    "adjusted_p_value": "padj",
    "q_value": "padj",
    "qvalue": "padj",
    "fdr": "padj",
    "q_ttest": "padj",
    "q_wilcoxon": "padj",
    "校正p值": "padj",
    "direction": "direction",
    "regulation": "direction",
    "change_direction": "direction",
    "trait": "trait",
    "gwas_trait": "trait",
    "accessionid": "accession_id",
    "accession_id": "accession_id",
    "gcst": "accession_id",
    "reportedtrait": "reported_trait",
    "reported_trait": "reported_trait",
    "reported_trait_name": "reported_trait",
    "summarystatistics": "summary_statistics_url",
    "summary_statistics": "summary_statistics_url",
    "summary_statistics_url": "summary_statistics_url",
    "summary_stats_url": "summary_statistics_url",
    "pubmedid": "pubmed_id",
    "pubmed_id": "pubmed_id",
    "pmid": "pubmed_id",
    "title": "paper_title",
    "paper_title": "paper_title",
    "journal": "journal",
    "publicationdate": "publication_date",
    "publication_date": "publication_date",
    "efotraits": "efo_traits",
    "efo_traits": "efo_traits",
    "bgtraits": "bg_traits",
    "bg_traits": "bg_traits",
    "initialsampledescription": "initial_sample_description",
    "initial_sample_description": "initial_sample_description",
    "discoverysampleancestry": "discovery_sample_ancestry",
    "discovery_sample_ancestry": "discovery_sample_ancestry",
    "celltype": "celltype",
    "cell_type": "cell_type",
    "celltype_l": "celltype_l",
    "cell_type_l": "celltype_l",
    "celltype_long": "celltype_l",
    "cell_type_long": "celltype_l",
    "group_col": "group_col",
    "group_column": "group_col",
    "group1": "group1",
    "group_1": "group1",
    "group2": "group2",
    "group_2": "group2",
    "n_group1": "n_group1",
    "n_group2": "n_group2",
    "mean_group1": "mean_group1",
    "mean_group2": "mean_group2",
    "median_group1": "median_group1",
    "median_group2": "median_group2",
    "sd_group1": "sd_group1",
    "sd_group2": "sd_group2",
    "mean_diff": "mean_diff",
    "median_diff": "median_diff",
    "pseudo_log2fc_shifted": "pseudo_log2FC_shifted",
    "pseudo_log2_fc_shifted": "pseudo_log2FC_shifted",
    "pseudo_logfc_shifted": "pseudo_log2FC_shifted",
    "signed_log2fc_absmean": "signed_log2FC_absmean",
    "signed_log2_fc_absmean": "signed_log2FC_absmean",
    "signed_logfc_absmean": "signed_log2FC_absmean",
    "cohen_d": "cohen_d",
    "cliffs_delta": "cliffs_delta",
    "t_welch": "t_welch",
    "z_welch": "z_welch",
    "u_stat": "u_stat",
    "z_wilcoxon": "z_wilcoxon",
    "趋势": "direction",
    "方向": "direction",
}

DIRECTION_CHAT_ALIASES = {
    "up": "up",
    "upregulated": "up",
    "increase": "up",
    "increased": "up",
    "higher": "up",
    "+": "up",
    "上调": "up",
    "升高": "up",
    "增加": "up",
    "down": "down",
    "downregulated": "down",
    "decrease": "down",
    "decreased": "down",
    "lower": "down",
    "-": "down",
    "下调": "down",
    "降低": "down",
    "减少": "down",
    "unchanged": "unchanged",
    "stable": "unchanged",
    "none": "unchanged",
    "不变": "unchanged",
    "无变化": "unchanged",
}


def canonical_chat_column(name: Any) -> str:
    raw = str(name or "").strip().lstrip("\ufeff")
    lowered = raw.casefold()
    if lowered in CHAT_COLUMN_ALIASES:
        return CHAT_COLUMN_ALIASES[lowered]
    token = re.sub(r"[^0-9a-z]+", "_", lowered).strip("_")
    return CHAT_COLUMN_ALIASES.get(token, token)


def normalize_chat_cell(field: str, value: Any) -> Any:
    if value is None:
        return ""
    if isinstance(value, str):
        text = value.strip()
        if text == "":
            return ""
        if field == "direction":
            return DIRECTION_CHAT_ALIASES.get(text.casefold(), DIRECTION_CHAT_ALIASES.get(text, text))
        if field in NUMERIC_CHAT_FIELDS:
            number = maybe_number(text)
            return number if number is not None else text
        return text
    if field == "direction":
        return DIRECTION_CHAT_ALIASES.get(str(value).strip().casefold(), value)
    return value


def normalize_chat_record(record: Any, row_index: int = 0) -> dict[str, Any]:
    if isinstance(record, str):
        text = record.strip()
        return {"input_id": f"row_{row_index}", "name": text} if text else {}
    if not isinstance(record, dict):
        return {}
    normalized: dict[str, Any] = {}
    for key, value in record.items():
        field = canonical_chat_column(key)
        if not field:
            continue
        cell = normalize_chat_cell(field, value)
        if cell == "" or cell is None:
            continue
        normalized[field] = cell
    normalized.setdefault("input_id", str(record.get("input_id") or record.get("id") or f"row_{row_index}"))
    return normalized


def normalize_gwas_accession(value: Any) -> str:
    token = re.sub(r"[^0-9A-Za-z]+", "", str(value or "")).upper()
    if token.startswith("CGST"):
        token = "GCST" + token[4:]
    return token


DIRECT_METABOLITE_RECORD_FIELDS = {
    "name",
    "metabolite",
    "compound",
    "hmdb",
    "hmdb_id",
    "chebi",
    "chebi_id",
    "pubchem",
    "pubchem_cid",
    "cid",
    "kegg",
    "kegg_id",
    "inchikey",
    "inchi_key",
    "formula",
    "mz",
    "ms2_peaks",
}

PRIMARY_COMPARISON_GROUP_TERMS = {
    "tumor",
    "tumour",
    "cancer",
    "malignant",
    "case",
    "disease",
    "lesion",
}

REFERENCE_COMPARISON_GROUP_TERMS = {
    "adjacent",
    "normal",
    "control",
    "healthy",
    "benign",
    "paratumor",
    "para_tumor",
    "non_tumor",
    "nontumor",
}


POOL_OR_CLASS_NAME_RE = re.compile(
    r"\b(?:lipid|sphingo|sphingosine|ceramide|acylcarnitine|carnitine|"
    r"phosphatidyl|lysophosphatidyl|plasmalogen|sterol|steroid|bile acid|"
    r"fatty acid|gpc|gpe|gpi|gps|gpg|d\d{1,2}:\d|\d{1,2}:\d|[a-z]{1,5}\(?\d{1,2}:\d)",
    flags=re.IGNORECASE,
)


def comparison_group_key(value: Any) -> str:
    return re.sub(r"[^0-9a-z]+", "_", str(value or "").strip().casefold()).strip("_")


def comparison_group_has_any(value: Any, terms: set[str]) -> bool:
    key = comparison_group_key(value)
    if not key:
        return False
    tokens = set(key.split("_"))
    return bool(tokens & terms or any(term in key for term in terms))


def normalized_record_lookup(record: Any) -> dict[str, Any]:
    if not isinstance(record, dict):
        return {}
    return {str(key).strip().casefold().replace(" ", "_"): value for key, value in record.items()}


def first_record_value(record: dict[str, Any], keys: Iterable[str]) -> Any:
    normalized = normalized_record_lookup(record)
    for key in keys:
        value = normalized.get(key)
        if value not in {None, ""}:
            return value
    return None


def has_direct_metabolite_identifier(record: dict[str, Any]) -> bool:
    normalized = normalized_record_lookup(record)
    return any(normalized.get(field) not in {None, ""} for field in DIRECT_METABOLITE_RECORD_FIELDS)


def trait_text_is_ratio_or_composite(value: Any) -> bool:
    text = clean_reported_trait_name(value)
    if not text:
        return False
    lowered = text.casefold()
    return bool(re.search(r"\bratio\b", lowered) or len(split_outside_parentheses(text)) > 1 or " + " in lowered)


def trait_text_is_pool_or_class(value: Any) -> bool:
    return bool(POOL_OR_CLASS_NAME_RE.search(str(value or "")))


def comparison_effect_display_label(effect_source: str) -> str:
    source = str(effect_source or "").strip()
    if source == "log2FC":
        return "log2FC"
    if source.startswith("pseudo_log2FC_shifted"):
        return "pseudo_log2FC_shifted"
    if source.startswith("signed_log2FC_absmean"):
        return "signed_log2FC_absmean"
    if source.startswith("mean_or_median_diff"):
        return "mean_or_median_diff"
    return source or "directional_effect"


def is_two_group_trait_comparison_record(record: Any) -> bool:
    if not isinstance(record, dict):
        return False
    normalized = normalized_record_lookup(record)
    if normalized.get("source_analysis_mode") == "two_group_trait_comparison_table":
        return False
    has_trait = any(normalized.get(field) not in {None, ""} for field in ("trait", "accession_id", "gcst", "gwas_trait"))
    has_groups = normalized.get("group1") not in {None, ""} and normalized.get("group2") not in {None, ""}
    has_stat = any(
        normalized.get(field) not in {None, ""}
        for field in (
            "log2fc",
            "log2_fc",
            "logfc",
            "fold_change",
            "pvalue",
            "p_value",
            "p_ttest",
            "p_wilcoxon",
            "padj",
            "q_ttest",
            "q_wilcoxon",
            "direction",
            "mean_diff",
            "median_diff",
        )
    )
    return bool(has_trait and has_groups and has_stat and not has_direct_metabolite_identifier(record))


def infer_primary_comparison_group(group1: Any, group2: Any) -> tuple[str, str]:
    group1_text = str(group1 or "").strip()
    group2_text = str(group2 or "").strip()
    if comparison_group_has_any(group1_text, PRIMARY_COMPARISON_GROUP_TERMS):
        return group1_text, "group1_primary_term"
    if comparison_group_has_any(group2_text, PRIMARY_COMPARISON_GROUP_TERMS):
        return group2_text, "group2_primary_term"
    if comparison_group_has_any(group1_text, REFERENCE_COMPARISON_GROUP_TERMS) and group2_text:
        return group2_text, "group1_reference_term"
    if comparison_group_has_any(group2_text, REFERENCE_COMPARISON_GROUP_TERMS) and group1_text:
        return group1_text, "group2_reference_term"
    return group1_text, "default_group1"


def direction_group_index(raw_direction: Any) -> int | None:
    text = comparison_group_key(raw_direction)
    if not text:
        return None
    group1_high = {
        "group1_hi",
        "group1_high",
        "group1_higher",
        "group_1_hi",
        "group_1_high",
        "g1_hi",
        "g1_high",
    }
    group2_high = {
        "group2_hi",
        "group2_high",
        "group2_higher",
        "group_2_hi",
        "group_2_high",
        "g2_hi",
        "g2_high",
    }
    group1_low = {"group1_low", "group1_lo", "group_1_low", "g1_low", "g1_lo"}
    group2_low = {"group2_low", "group2_lo", "group_2_low", "g2_low", "g2_lo"}
    if text in group1_high or text in group2_low:
        return 1
    if text in group2_high or text in group1_low:
        return 2
    return None


def normalize_comparison_p_value(value: Any) -> Any:
    number = maybe_number(value)
    if number is None:
        return value
    return 1e-300 if number == 0 else number


def comparison_effect_fallback(record: dict[str, Any], high_group_index: int | None, primary_is_group2: bool) -> tuple[float | None, str]:
    pseudo_value = maybe_number(first_present(record, ("pseudo_log2fc_shifted", "pseudo_log2_fc_shifted", "pseudo_logfc_shifted")))
    if pseudo_value is not None:
        magnitude = abs(pseudo_value)
        if high_group_index is not None:
            primary_high = (high_group_index == 2) if primary_is_group2 else (high_group_index == 1)
            return (magnitude if primary_high else -magnitude), "pseudo_log2FC_shifted_oriented_by_group_direction"
        return pseudo_value, "pseudo_log2FC_shifted"
    signed_value = maybe_number(first_present(record, ("signed_log2fc_absmean", "signed_log2_fc_absmean", "signed_logfc_absmean")))
    if signed_value is not None:
        return (-signed_value if primary_is_group2 else signed_value), "signed_log2FC_absmean_oriented"
    mean_diff = maybe_number(first_present(record, ("mean_diff", "median_diff")))
    if mean_diff is not None:
        return (-mean_diff if primary_is_group2 else mean_diff), "mean_or_median_diff_as_directional_effect"
    return None, ""


def convert_two_group_trait_comparison_record(record: dict[str, Any], row_index: int = 0) -> dict[str, Any]:
    normalized = normalized_record_lookup(record)
    group1 = first_present(normalized, ("group1", "group_1"))
    group2 = first_present(normalized, ("group2", "group_2"))
    primary_group, primary_reason = infer_primary_comparison_group(group1, group2)
    primary_is_group2 = comparison_group_key(primary_group) == comparison_group_key(group2)
    log2_fc = parse_log2_fold_change(normalized)
    raw_direction = first_present(normalized, ("direction", "regulation", "change_direction"))
    high_group_index = direction_group_index(raw_direction)
    effect_source = "log2FC"
    if log2_fc is None:
        log2_fc, effect_source = comparison_effect_fallback(normalized, high_group_index, primary_is_group2)
        oriented_log2_fc = log2_fc
    else:
        oriented_log2_fc = -log2_fc if primary_is_group2 else log2_fc
    if high_group_index is not None:
        high_group = group1 if high_group_index == 1 else group2
        direction = "up" if comparison_group_key(high_group) == comparison_group_key(primary_group) else "down"
    else:
        direction = direction_from_features(raw_direction, oriented_log2_fc)

    converted = dict(record)
    accession = normalize_gwas_accession(first_present(normalized, ("accession_id", "trait", "gcst", "gwas_trait")))
    if accession:
        converted["accession_id"] = accession
        converted.setdefault("trait", accession)
    converted["log2FC"] = round(oriented_log2_fc, 6) if oriented_log2_fc is not None else converted.get("log2FC")
    p_value = first_present(normalized, ("pvalue", "p_value", "p", "p_ttest", "p_wilcoxon"))
    adjusted_p = first_present(normalized, ("padj", "adjusted_p_value", "q_value", "qvalue", "fdr", "q_wilcoxon", "q_ttest"))
    if p_value not in {None, ""}:
        converted["pvalue"] = normalize_comparison_p_value(p_value)
    if adjusted_p not in {None, ""}:
        converted["padj"] = normalize_comparison_p_value(adjusted_p)
    converted["direction"] = direction
    if normalized.get("cell_type") in {None, ""}:
        cell_type = first_present(normalized, ("celltype_l", "cell_type_l", "celltype"))
        if cell_type not in {None, ""}:
            converted["cell_type"] = cell_type
    converted.setdefault("input_id", f"{accession or first_present(normalized, ('trait',)) or 'trait'}|{converted.get('cell_type') or normalized.get('celltype') or normalized.get('celltype_l') or 'cell'}|{group1}_vs_{group2}|row_{row_index}")
    converted["source_analysis_mode"] = "two_group_trait_comparison_table"
    converted["comparison_primary_group"] = primary_group
    converted["comparison_reference_group"] = str(group2 if comparison_group_key(primary_group) == comparison_group_key(group1) else group1 or "").strip()
    converted["comparison_primary_group_reason"] = primary_reason
    converted["comparison_original_direction"] = raw_direction or ""
    if log2_fc is not None:
        converted["comparison_original_log2FC"] = round(log2_fc, 6)
    converted["comparison_log2fc_orientation"] = "primary_vs_reference"
    converted["comparison_effect_source"] = effect_source
    converted["effect_label"] = comparison_effect_display_label(effect_source)
    if oriented_log2_fc is not None:
        converted["effect_value"] = round(oriented_log2_fc, 6)
    if converted["effect_label"] != "log2FC":
        converted["log2FC_semantics"] = "directional_effect_surrogate_not_abundance_log2fc"
    trait_label = first_present(normalized, ("reported_trait", "reportedtrait", "name", "metabolite", "compound"))
    if trait_text_is_ratio_or_composite(trait_label):
        converted["trait_value_type"] = "ratio_or_composite_trait"
        ratio_trait = ratio_trait_descriptor(
            trait_label,
            comparison_effect_display_label(effect_source),
            converted.get("effect_value", converted.get("log2FC")),
            converted.get("direction", ""),
        )
        if ratio_trait:
            converted["ratio_trait"] = ratio_trait
    elif trait_text_is_pool_or_class(trait_label):
        converted["trait_value_type"] = "class_or_pool_trait"
    else:
        converted["trait_value_type"] = "trait_score"
    return normalize_chat_record(converted, row_index)


def prepare_analysis_records(
    records: list[Any],
    source_metadata: dict[str, Any] | None = None,
) -> tuple[list[Any], dict[str, Any]]:
    metadata: dict[str, Any] = {
        "source": (source_metadata or {}).get("source", "records"),
        "record_count": 0,
        "columns": sorted({key for record in records if isinstance(record, dict) for key in record}),
        "warnings": list((source_metadata or {}).get("warnings", [])),
        "analysis_mode": "metabolite_table",
        "input_format_detected": "metabolite_table",
        "converted_record_count": 0,
    }
    if not records:
        return [], metadata
    converted_records: list[Any] = []
    comparison_indices = [index for index, record in enumerate(records) if is_two_group_trait_comparison_record(record)]
    existing_modes = {
        str(record.get("source_analysis_mode"))
        for record in records
        if isinstance(record, dict) and record.get("source_analysis_mode")
    }
    if comparison_indices:
        for index, record in enumerate(records):
            if index in comparison_indices and isinstance(record, dict):
                converted_records.append(convert_two_group_trait_comparison_record(record, index))
            else:
                converted_records.append(record)
        metadata["analysis_mode"] = "two_group_trait_comparison_table"
        metadata["input_format_detected"] = "two_group_trait_comparison_table"
        metadata["converted_record_count"] = len(comparison_indices)
        metadata["warnings"].append(
            {
                "code": "trait_score_comparison_not_direct_abundance",
                "severity": "warning",
                "message": (
                    "Two-group trait score comparison detected. Directional effects may be pseudo_log2FC_shifted, "
                    "mean_diff, or signed score differences, not direct LC-MS abundance log2FC."
                ),
            }
        )
        metadata["warnings"].append(
            {
                "code": "cell_level_trait_score_requires_patient_level_validation",
                "severity": "info",
                "message": (
                    "If rows were computed from single cells, prioritize effect sizes and validate with "
                    "patient-level pseudobulk or mixed models."
                ),
            }
        )
        primary_counts = Counter(
            str(record.get("comparison_primary_group") or "")
            for record in converted_records
            if isinstance(record, dict) and record.get("comparison_primary_group")
        )
        metadata["comparison_primary_group_counts"] = dict(sorted(primary_counts.items()))
        if len(comparison_indices) != len(records):
            metadata["warnings"].append(
                {
                    "code": "mixed_input_modes",
                    "message": "Some rows looked like two-group trait comparisons while others looked like direct metabolite records; both were retained.",
                    "comparison_rows": len(comparison_indices),
                    "total_rows": len(records),
                }
            )
    else:
        converted_records = list(records)
        if existing_modes:
            metadata["analysis_mode"] = sorted(existing_modes)[0]
            metadata["input_format_detected"] = sorted(existing_modes)[0]
    metadata["record_count"] = len(converted_records)
    metadata["columns"] = sorted({key for record in converted_records if isinstance(record, dict) for key in record})
    return converted_records, metadata


TRAIT_SCORE_Q_FIELDS = (
    "q_wilcoxon_global",
    "q_ttest_global",
    "q_wilcoxon",
    "q_ttest",
    "padj",
    "adjusted_p_value",
    "q_value",
    "qvalue",
    "fdr",
)
TRAIT_SCORE_EFFECT_FIELDS = (
    "cohen_d",
    "mean_diff",
    "median_diff",
    "z_wilcoxon",
    "z_welch_approx",
    "pseudo_log2fc_shifted",
    "signed_log2fc_absmean",
    "log2fc",
)


DIFFERENTIAL_Q_FIELDS = (
    "q_wilcoxon_global",
    "q_ttest_global",
    "padj",
    "adjusted_p_value",
    "q_value",
    "qvalue",
    "fdr",
    "q_wilcoxon",
    "q_ttest",
)
DIFFERENTIAL_P_FIELDS = (
    "pvalue",
    "p_value",
    "p",
    "p_wilcoxon",
    "p_ttest",
)
DIFFERENTIAL_EFFECT_FIELDS = (
    ("log2fc", "log2FC"),
    ("log2_fc", "log2FC"),
    ("logfc", "log2FC"),
    ("log_fc", "log2FC"),
    ("log2_fold_change", "log2FC"),
    ("fold_change", "fold_change"),
    ("foldchange", "fold_change"),
    ("fc", "fold_change"),
    ("pseudo_log2fc_shifted", "pseudo_log2FC_shifted"),
    ("signed_log2fc_absmean", "signed_log2FC_absmean"),
    ("mean_diff", "mean_diff"),
    ("median_diff", "median_diff"),
    ("cohen_d", "cohen_d"),
    ("cliffs_delta", "cliffs_delta"),
    ("z_wilcoxon", "z_wilcoxon"),
    ("z_welch", "z_welch"),
    ("z_welch_approx", "z_welch_approx"),
    ("t_welch", "t_welch"),
    ("statistic", "statistic"),
    ("score", "score"),
    ("diff", "diff"),
    ("delta", "delta"),
)
def trait_score_numeric_value(record: dict[str, Any], fields: Iterable[str]) -> tuple[float | None, str]:
    normalized = normalized_record_lookup(record)
    for field in fields:
        value = normalized.get(field)
        if value in {None, ""}:
            continue
        number = maybe_number(value)
        if number is not None:
            return number, field
    return None, ""


def differential_numeric_value(record: dict[str, Any], fields: Iterable[str]) -> tuple[float | None, str]:
    normalized = normalized_record_lookup(record)
    for field in fields:
        value = normalized.get(field)
        if value in {None, ""}:
            continue
        number = maybe_number(value)
        if number is not None:
            return number, field
    return None, ""


def differential_effect_value(record: dict[str, Any]) -> tuple[float | None, str, bool]:
    normalized = normalized_record_lookup(record)
    log2_fc = parse_log2_fold_change(normalized)
    if log2_fc is not None:
        return log2_fc, "log2FC", True
    for field, label in DIFFERENTIAL_EFFECT_FIELDS:
        value = normalized.get(field)
        if value in {None, ""}:
            continue
        number = maybe_number(value)
        if number is None:
            continue
        if label == "fold_change":
            if number <= 0:
                continue
            return math.log2(number), "fold_change_as_log2FC", True
        return number, label, False
    return None, "", False


def is_differential_table_record(record: Any) -> bool:
    if not isinstance(record, dict):
        return False
    normalized = normalized_record_lookup(record)
    has_effect = differential_effect_value(record)[0] is not None
    has_significance = any(normalized.get(field) not in {None, ""} for field in (*DIFFERENTIAL_Q_FIELDS, *DIFFERENTIAL_P_FIELDS))
    has_identifier = has_direct_metabolite_identifier(record) or any(
        normalized.get(field) not in {None, ""}
        for field in ("trait", "accession_id", "gcst", "gwas_trait", "reported_trait")
    )
    return bool(has_identifier and (has_effect or has_significance))


def normalize_differential_table_record(record: dict[str, Any], row_index: int = 0) -> dict[str, Any]:
    if is_two_group_trait_comparison_record(record):
        converted = convert_two_group_trait_comparison_record(record, row_index)
        converted["source_analysis_mode"] = "differential_table"
        converted["differential_input_subtype"] = "two_group_trait_comparison"
        return converted

    normalized = normalized_record_lookup(record)
    effect_value, effect_label, is_log2_effect = differential_effect_value(record)
    p_value, p_field = differential_numeric_value(record, DIFFERENTIAL_P_FIELDS)
    q_value, q_field = differential_numeric_value(record, DIFFERENTIAL_Q_FIELDS)
    converted = dict(record)
    if effect_value is not None:
        converted["log2FC"] = round(effect_value, 6)
        converted["effect_value"] = round(effect_value, 6)
        converted["effect_label"] = effect_label
        converted["comparison_effect_source"] = effect_label
        if not is_log2_effect:
            converted["log2FC_semantics"] = "directional_effect_surrogate_not_abundance_log2fc"
    if p_value is not None:
        converted["pvalue"] = normalize_comparison_p_value(p_value)
        converted["differential_p_field"] = p_field
    if q_value is not None:
        converted["padj"] = normalize_comparison_p_value(q_value)
        converted["differential_q_field"] = q_field
    converted["direction"] = direction_from_features(first_present(normalized, ("direction", "regulation", "change_direction")), effect_value)
    converted.setdefault("input_id", str(first_present(normalized, ("input_id", "id")) or f"differential_row_{row_index}"))
    converted["source_analysis_mode"] = "differential_table"
    converted["differential_input_subtype"] = "metabolite_differential_result"
    converted["differential_table_semantics"] = "precomputed_group_difference_not_raw_abundance"
    ratio_trait = ratio_trait_descriptor(
        first_present(normalized, ("reported_trait", "reportedtrait", "name", "metabolite", "compound")),
        str(converted.get("effect_label") or ""),
        converted.get("effect_value", converted.get("log2FC")),
        str(converted.get("direction") or ""),
    )
    if ratio_trait:
        converted["ratio_trait"] = ratio_trait
    return normalize_chat_record(converted, row_index)


def differential_group_key(record: dict[str, Any]) -> str:
    normalized = normalized_record_lookup(record)
    return str(
        first_present(normalized, ("celltype_l", "cell_type_l", "celltype", "cell_type", "group_col"))
        or "ALL_ROWS"
    ).strip() or "ALL_ROWS"


def select_differential_table_records(
    records: list[Any],
    *,
    q_threshold: float | None = 0.05,
    p_threshold: float | None = None,
    min_abs_effect: float = 0.0,
    top_per_group: int = 120,
    max_records: int = 300,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    eligible: list[dict[str, Any]] = []
    rejected = Counter()
    input_subtypes = Counter()
    for index, record in enumerate(records):
        if not isinstance(record, dict):
            rejected["not_object"] += 1
            continue
        if not is_differential_table_record(record):
            rejected["not_differential_record"] += 1
            continue
        normalized_record = normalize_differential_table_record(record, index)
        effect_value, effect_field = differential_numeric_value(normalized_record, ("effect_value", "log2fc", "log2_fc", "logfc"))
        if effect_value is None:
            rejected["missing_effect"] += 1
            continue
        q_value, q_field = differential_numeric_value(normalized_record, DIFFERENTIAL_Q_FIELDS)
        p_value, p_field = differential_numeric_value(normalized_record, DIFFERENTIAL_P_FIELDS)
        if q_threshold is not None and q_value is not None and q_value > q_threshold:
            rejected["q_above_threshold"] += 1
            continue
        if q_value is None and p_threshold is not None and p_value is not None and p_value > p_threshold:
            rejected["p_above_threshold"] += 1
            continue
        if abs(effect_value) < min_abs_effect:
            rejected["effect_below_threshold"] += 1
            continue
        row = dict(normalized_record)
        row["_differential_source_index"] = index
        row["_differential_q_value"] = q_value
        row["_differential_q_field"] = q_field
        row["_differential_p_value"] = p_value
        row["_differential_p_field"] = p_field
        row["_differential_effect_value"] = effect_value
        row["_differential_effect_field"] = effect_field
        input_subtypes[str(row.get("differential_input_subtype") or "unknown")] += 1
        eligible.append(row)

    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in eligible:
        grouped[differential_group_key(row)].append(row)

    def selection_sort_p(row: dict[str, Any]) -> float:
        q_value = row.get("_differential_q_value")
        p_value = row.get("_differential_p_value")
        if q_value is not None:
            return float(q_value)
        if p_value is not None:
            return float(p_value)
        return 1.0

    selected: list[dict[str, Any]] = []
    per_group_selected: dict[str, int] = {}
    for group, rows in sorted(grouped.items()):
        ranked = sorted(
            rows,
            key=lambda row: (
                selection_sort_p(row),
                -abs(float(row.get("_differential_effect_value") or 0.0)),
                int(row.get("_differential_source_index") or 0),
            ),
        )
        limit = top_per_group if top_per_group and top_per_group > 0 else len(ranked)
        kept = ranked[:limit]
        per_group_selected[group] = len(kept)
        selected.extend(kept)

    selected = sorted(
        selected,
        key=lambda row: (
            selection_sort_p(row),
            -abs(float(row.get("_differential_effect_value") or 0.0)),
            int(row.get("_differential_source_index") or 0),
        ),
    )
    if max_records and max_records > 0:
        selected = selected[:max_records]

    output: list[dict[str, Any]] = []
    for rank, row in enumerate(selected, start=1):
        clean = {key: value for key, value in row.items() if not str(key).startswith("_differential_")}
        clean["differential_selection_rank"] = rank
        clean["differential_selection_q_value"] = row.get("_differential_q_value")
        clean["differential_selection_q_field"] = row.get("_differential_q_field", "")
        clean["differential_selection_p_value"] = row.get("_differential_p_value")
        clean["differential_selection_p_field"] = row.get("_differential_p_field", "")
        clean["differential_selection_effect_value"] = row.get("_differential_effect_value")
        clean["differential_selection_effect_field"] = row.get("_differential_effect_field", "")
        output.append(clean)

    selected_groups = Counter(differential_group_key(row) for row in output)
    summary = {
        "contract_version": "differential_table_selection.v1",
        "input_count": len(records),
        "eligible_count": len(eligible),
        "selected_count": len(output),
        "q_threshold": q_threshold,
        "p_threshold": p_threshold,
        "min_abs_effect": min_abs_effect,
        "top_per_group": top_per_group,
        "max_records": max_records,
        "input_subtypes": dict(sorted(input_subtypes.items())),
        "selected_groups": dict(sorted(selected_groups.items())),
        "pre_cap_selected_by_group": dict(sorted(per_group_selected.items())),
        "rejected_counts": dict(sorted(rejected.items())),
        "selection_rule": "q or adjusted p ascending when present, then absolute effect descending; raw abundance columns are not required or read.",
        "interpretation_boundary": "Rows are treated as precomputed differential results. Effect fields drive seed direction and weight, not raw metabolite abundance.",
    }
    return output, summary


def trait_score_celltype_key(record: dict[str, Any]) -> str:
    normalized = normalized_record_lookup(record)
    return str(
        first_present(normalized, ("celltype_l", "cell_type_l", "celltype", "cell_type"))
        or "ALL_CELLTYPES"
    ).strip() or "ALL_CELLTYPES"


def select_trait_score_records(
    records: list[Any],
    *,
    q_threshold: float | None = 0.05,
    min_abs_effect: float = 0.0,
    top_per_celltype: int = 80,
    max_records: int = 120,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    eligible: list[dict[str, Any]] = []
    rejected = Counter()
    for index, record in enumerate(records):
        if not isinstance(record, dict) or not is_two_group_trait_comparison_record(record):
            rejected["not_trait_score_comparison"] += 1
            continue
        q_value, q_field = trait_score_numeric_value(record, TRAIT_SCORE_Q_FIELDS)
        effect_value, effect_field = trait_score_numeric_value(record, TRAIT_SCORE_EFFECT_FIELDS)
        if q_threshold is not None and q_value is not None and q_value > q_threshold:
            rejected["q_above_threshold"] += 1
            continue
        if effect_value is None:
            rejected["missing_effect"] += 1
            continue
        if abs(effect_value) < min_abs_effect:
            rejected["effect_below_threshold"] += 1
            continue
        row = dict(record)
        row["_trait_score_source_index"] = index
        row["_trait_score_q_value"] = q_value
        row["_trait_score_q_field"] = q_field
        row["_trait_score_effect_value"] = effect_value
        row["_trait_score_effect_field"] = effect_field
        eligible.append(row)

    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in eligible:
        grouped[trait_score_celltype_key(row)].append(row)

    def selection_sort_q(row: dict[str, Any]) -> float:
        value = row.get("_trait_score_q_value")
        return 1.0 if value is None else float(value)

    selected: list[dict[str, Any]] = []
    per_celltype_selected: dict[str, int] = {}
    for celltype, rows in sorted(grouped.items()):
        ranked = sorted(
            rows,
            key=lambda row: (
                selection_sort_q(row),
                -abs(float(row.get("_trait_score_effect_value") or 0.0)),
                int(row.get("_trait_score_source_index") or 0),
            ),
        )
        limit = top_per_celltype if top_per_celltype and top_per_celltype > 0 else len(ranked)
        kept = ranked[:limit]
        per_celltype_selected[celltype] = len(kept)
        selected.extend(kept)

    selected = sorted(
        selected,
        key=lambda row: (
            selection_sort_q(row),
            -abs(float(row.get("_trait_score_effect_value") or 0.0)),
            int(row.get("_trait_score_source_index") or 0),
        ),
    )
    if max_records and max_records > 0:
        selected = selected[:max_records]

    output: list[dict[str, Any]] = []
    for rank, row in enumerate(selected, start=1):
        clean = {key: value for key, value in row.items() if not str(key).startswith("_trait_score_")}
        normalized = normalized_record_lookup(clean)
        accession = normalize_gwas_accession(first_present(normalized, ("accession_id", "trait", "gcst", "gwas_trait")))
        celltype = trait_score_celltype_key(clean)
        clean.setdefault("input_id", f"{accession or 'trait'}|{celltype}|trait_score_rank_{rank}")
        clean["trait_score_selection_rank"] = rank
        clean["trait_score_selection_q_value"] = row.get("_trait_score_q_value")
        clean["trait_score_selection_q_field"] = row.get("_trait_score_q_field", "")
        clean["trait_score_selection_effect_value"] = row.get("_trait_score_effect_value")
        clean["trait_score_selection_effect_field"] = row.get("_trait_score_effect_field", "")
        output.append(clean)

    selected_celltypes = Counter(trait_score_celltype_key(row) for row in output)
    summary = {
        "contract_version": "trait_score_selection.v1",
        "input_count": len(records),
        "eligible_count": len(eligible),
        "selected_count": len(output),
        "q_threshold": q_threshold,
        "min_abs_effect": min_abs_effect,
        "top_per_celltype": top_per_celltype,
        "max_records": max_records,
        "selected_celltypes": dict(sorted(selected_celltypes.items())),
        "pre_cap_selected_by_celltype": dict(sorted(per_celltype_selected.items())),
        "rejected_counts": dict(sorted(rejected.items())),
        "selection_rule": "q ascending, then absolute effect descending; q defaults to 1 only for records with no q column.",
    }
    return output, summary


def infer_trait_score_context(records: list[dict[str, Any]]) -> dict[str, Any]:
    celltypes = [trait_score_celltype_key(row) for row in records]
    celltype_counts = Counter(celltypes)
    group_pairs = Counter(
        (
            str(normalized_record_lookup(row).get("group1") or "").strip(),
            str(normalized_record_lookup(row).get("group2") or "").strip(),
        )
        for row in records
        if isinstance(row, dict)
    )
    top_celltype = celltype_counts.most_common(1)[0][0] if celltype_counts else ""
    tokens = [token for token in re.split(r"[_\s]+", top_celltype) if token]
    terms = dedupe_preserve_order([top_celltype, *tokens])
    group1, group2 = group_pairs.most_common(1)[0][0] if group_pairs else ("", "")
    context = {
        "cell_type": top_celltype,
        "context_terms": terms,
        "comparison": f"{group1} vs {group2}".strip(),
        "analysis_mode": "two_group_trait_comparison_table",
    }
    if tokens:
        context["cancer_type"] = tokens[0]
    return {
        key: value
        for key, value in context.items()
        if value is not None and value != "" and value != []
    }


def infer_differential_table_context(records: list[dict[str, Any]]) -> dict[str, Any]:
    groups = Counter(differential_group_key(row) for row in records)
    group_pairs = Counter(
        (
            str(normalized_record_lookup(row).get("group1") or "").strip(),
            str(normalized_record_lookup(row).get("group2") or "").strip(),
        )
        for row in records
        if isinstance(row, dict)
    )
    top_group = groups.most_common(1)[0][0] if groups else ""
    tokens = [token for token in re.split(r"[_\s]+", top_group) if token]
    context = {
        "cell_type": top_group if top_group != "ALL_ROWS" else "",
        "context_terms": dedupe_preserve_order([top_group, *tokens]) if top_group != "ALL_ROWS" else [],
        "analysis_mode": "differential_table",
    }
    if tokens:
        context["cancer_type"] = tokens[0]
    if group_pairs:
        group1, group2 = group_pairs.most_common(1)[0][0]
        if group1 or group2:
            context["comparison"] = f"{group1} vs {group2}".strip()
    return {
        key: value
        for key, value in context.items()
        if value is not None and value != "" and value != []
    }


def split_outside_parentheses(text: str, separator: str = " to ") -> list[str]:
    depth = 0
    lowered = text.casefold()
    sep = separator.casefold()
    for index, char in enumerate(text):
        if char == "(":
            depth += 1
        elif char == ")" and depth:
            depth -= 1
        if depth == 0 and lowered.startswith(sep, index):
            return [text[:index], text[index + len(separator) :]]
    return [text]


def split_all_outside_parentheses(text: str, separator: str = " to ") -> list[str]:
    depth = 0
    lowered = text.casefold()
    sep = separator.casefold()
    parts: list[str] = []
    start = 0
    index = 0
    while index < len(text):
        char = text[index]
        if char == "(":
            depth += 1
        elif char == ")" and depth:
            depth -= 1
        if depth == 0 and lowered.startswith(sep, index):
            parts.append(text[start:index].strip())
            index += len(separator)
            start = index
            continue
        index += 1
    parts.append(text[start:].strip())
    return [part for part in parts if part]


def clean_reported_trait_name(value: Any) -> str:
    text = re.sub(r"\s+", " ", str(value or "").strip())
    text = re.sub(r"\s+levels$", "", text, flags=re.IGNORECASE).strip()
    text = re.sub(r"\s+\[\d+\]$", "", text).strip()
    return text


def dedupe_preserve_order(values: Iterable[Any]) -> list[str]:
    seen: set[str] = set()
    output: list[str] = []
    for value in values:
        text = str(value or "").strip()
        if not text:
            continue
        key = text.casefold()
        if key in seen:
            continue
        seen.add(key)
        output.append(text)
    return output


def ratio_trait_component_names(value: Any) -> list[str]:
    text = clean_reported_trait_name(value)
    ratio_match = re.search(r"\s+ratio$", text, flags=re.IGNORECASE)
    if not ratio_match:
        return []
    ratio_body = text[: ratio_match.start()].strip()
    return [clean_reported_trait_name(part) for part in split_all_outside_parentheses(ratio_body) if clean_reported_trait_name(part)]


def ratio_component_name_variants(value: Any) -> list[str]:
    text = clean_reported_trait_name(value)
    variants = [text, clean_reported_trait_name(re.sub(r"\([^)]*\)", "", text))]
    variants.extend(clean_reported_trait_name(match) for match in re.findall(r"\(([^)]*)\)", text))
    return dedupe_preserve_order(variant for variant in variants if variant)


def ratio_trait_descriptor(value: Any, effect_label: str = "", effect_value: Any = None, direction: str = "") -> dict[str, Any]:
    label = clean_reported_trait_name(value)
    components = ratio_trait_component_names(label)
    if len(components) < 2:
        return {}
    dual = len(components) == 2
    denominator_weight = 0.5 if dual else round(1.0 / max(1, len(components)), 6)
    ratio_components = []
    for index, component in enumerate(components):
        if dual and index == 0:
            role = "numerator"
            signed = "same_as_ratio"
        elif dual and index == 1:
            role = "denominator"
            signed = "opposite_to_ratio"
        else:
            role = f"component_{index + 1}"
            signed = "undirected_composite"
        ratio_components.append(
            {
                "name": component,
                "role": role,
                "signed_ratio_direction": signed,
                "component_weight": denominator_weight,
                "name_variants": ratio_component_name_variants(component),
            }
        )
    return {
        "contract_version": "ratio_trait.v1",
        "label": label,
        "component_count": len(ratio_components),
        "dual_component_ratio": dual,
        "components": ratio_components,
        "effect_label": effect_label,
        "effect_value": effect_value,
        "ratio_direction": direction,
        "interpretation_boundary": (
            "Ratio change constrains relative component balance only; it does not prove numerator abundance "
            "increased or denominator abundance decreased."
        ),
        "direction_policy": (
            "Dual-component ratios assign weak same/opposite ratio-consistent directions to numerator and "
            "denominator. Multi-component ratios stay undirected composite traits."
        ),
    }


def ratio_trait_for_record(record: Any) -> dict[str, Any]:
    if not isinstance(record, dict):
        return {}
    existing = record.get("ratio_trait")
    if isinstance(existing, dict) and existing.get("components"):
        return existing
    normalized = normalized_record_lookup(record)
    label = first_present(normalized, ("reported_trait", "reportedtrait", "name", "metabolite", "compound"))
    effect_value = first_present(normalized, ("effect_value", "log2fc", "log2_fc", "mean_diff", "cohen_d", "z_wilcoxon"))
    effect_label = first_present(normalized, ("effect_label", "comparison_effect_source")) or ""
    direction = first_present(normalized, ("direction", "regulation", "change_direction")) or ""
    return ratio_trait_descriptor(label, str(effect_label or ""), effect_value, str(direction or ""))


def ratio_component_match_score(candidate: dict[str, Any], component: dict[str, Any]) -> float:
    variants = [normalize_lookup_key(value) for value in component.get("name_variants", [])]
    variants = [variant for variant in variants if variant]
    if not variants:
        return 0.0
    candidate_texts = [candidate.get("display_name", ""), candidate.get("primary_external_id", "")]
    for match in candidate.get("matches", []) or []:
        candidate_texts.extend([match.get("raw_value", ""), match.get("field", ""), match.get("match_field", "")])
    candidate_keys = [normalize_lookup_key(value) for value in candidate_texts if value not in {None, ""}]
    for variant in variants:
        for key in candidate_keys:
            if not key:
                continue
            if key == variant:
                return 1.0
            if len(variant) >= 5 and (variant in key or key in variant):
                return max(0.5, min(len(variant), len(key)) / max(len(variant), len(key)))
    return 0.0


def ratio_component_for_candidate(record: Any, candidate: dict[str, Any]) -> dict[str, Any]:
    ratio = ratio_trait_for_record(record)
    components = ratio.get("components") or []
    if not components:
        return {}
    scored = [
        (ratio_component_match_score(candidate, component), index, component)
        for index, component in enumerate(components)
    ]
    scored.sort(key=lambda item: (-item[0], item[1]))
    best_score, index, component = scored[0]
    if best_score <= 0.0:
        if len(components) == 1:
            component = components[0]
        else:
            return {
                "ratio_trait": {key: value for key, value in ratio.items() if key != "components"},
                "component": {
                    "name": "",
                    "role": "unmatched_ratio_component",
                    "signed_ratio_direction": "unknown",
                    "component_weight": 0.0,
                    "match_score": 0.0,
                },
            }
    return {
        "ratio_trait": {key: value for key, value in ratio.items() if key != "components"},
        "component": {
            "name": component.get("name", ""),
            "role": component.get("role", ""),
            "signed_ratio_direction": component.get("signed_ratio_direction", "unknown"),
            "component_weight": component.get("component_weight", 0.0),
            "component_index": index,
            "match_score": round(best_score, 6),
        },
    }


def opposite_direction(value: Any) -> str:
    direction = str(value or "").strip().casefold()
    if direction == "up":
        return "down"
    if direction == "down":
        return "up"
    return "unknown"


LIPID_CLASS_PATTERNS = (
    ("phosphatidylcholine_gpc", re.compile(r"\b(?:gpc|phosphatidylcholine|phosphocholine|pc)\b", re.IGNORECASE)),
    ("sphingomyelin", re.compile(r"\b(?:sphingomyelin|sm)\b", re.IGNORECASE)),
    ("ceramide", re.compile(r"\b(?:ceramide|cer)\b", re.IGNORECASE)),
    ("sphingosine_sphingoid", re.compile(r"\b(?:sphingosine|sphinganine|sphingoid)\b", re.IGNORECASE)),
    ("fatty_acid", re.compile(r"\b(?:fatty acid|linoleoyl|palmitoyl|stearoyl|oleoyl|linoleate|palmitate|stearate)\b", re.IGNORECASE)),
    ("bile_acid", re.compile(r"\b(?:bile acid|cholate|deoxycholate|lithocholate)\b", re.IGNORECASE)),
    ("sterol_steroid", re.compile(r"\b(?:sterol|steroid|cholesterol)\b", re.IGNORECASE)),
)


def class_seed_descriptor(record: Any, candidate: dict[str, Any] | None = None) -> dict[str, Any]:
    if not isinstance(record, dict):
        return {}
    normalized = normalized_record_lookup(record)
    label = str(first_present(normalized, ("reported_trait", "reportedtrait", "name", "metabolite", "compound")) or "").strip()
    candidate = candidate or {}
    candidate_label = str(candidate.get("display_name") or "").strip()
    combined = " ".join(value for value in (label, candidate_label) if value)
    if not combined:
        return {}
    class_hits = [class_id for class_id, pattern in LIPID_CLASS_PATTERNS if pattern.search(combined)]
    chain_constraints = dedupe_preserve_order(re.findall(r"\b(?:d?\d{1,2}:\d(?:/\d{1,2}:\d)?(?:\([EZ0-9,]+\))?|[OP]-?\d{1,2}:\d)\b", combined, flags=re.IGNORECASE))
    stable_fields = [
        field
        for field in ("hmdb", "hmdb_id", "chebi", "chebi_id", "pubchem_cid", "cid", "kegg", "kegg_id", "inchikey", "inchi_key")
        if normalized.get(field) not in {None, ""}
    ]
    candidate_keys = candidate_identity_keys(candidate)
    structural_keys = [key for key in candidate_keys if key.startswith(("inchikey:", "inchikey_connectivity:", "hmdb:", "pubchem:"))]
    has_stable_identity_signal = bool(stable_fields or structural_keys)
    if has_stable_identity_signal and chain_constraints:
        allowed = "candidate_lipid_species_with_chain_constraints"
    elif has_stable_identity_signal:
        allowed = "candidate_identity_with_class_review"
    else:
        allowed = "class_or_pathway_only"
    return {
        "contract_version": "class_seed.v1",
        "label": label or candidate_label,
        "candidate_label": candidate_label,
        "class_ids": class_hits or ["class_or_pool"],
        "chain_constraints": chain_constraints,
        "stable_identifier_fields": stable_fields,
        "candidate_identity_keys": structural_keys[:8],
        "allowed_claim_level": allowed,
        "interpretation_boundary": (
            "Class or pool traits support class/pathway interpretation unless stable identifiers and chain constraints "
            "support a candidate lipid species; they do not by themselves prove a unique structural isomer."
        ),
    }


def candidate_identity_keys(candidate: dict[str, Any]) -> list[str]:
    keys: list[str] = []
    for xref in parse_list(candidate.get("external_xrefs")):
        text = str(xref or "").strip()
        if not text:
            continue
        upper = text.upper()
        if upper.startswith("INCHIKEY:"):
            value = upper.split(":", 1)[1]
            keys.append(f"inchikey:{value}")
            conn = inchikey_connectivity(value)
            if conn:
                keys.append(f"inchikey_connectivity:{conn}")
        elif upper.startswith("HMDB:"):
            keys.append(f"hmdb:{upper.split(':', 1)[1]}")
        elif upper.startswith(("PUBCHEM.COMPOUND:", "PUBCHEM:")):
            keys.append(f"pubchem:{upper.rsplit(':', 1)[1]}")
        elif upper.startswith("CHEBI:"):
            keys.append(f"chebi:{upper.split(':', 1)[1]}")
        elif upper.startswith("KEGG.COMPOUND:"):
            keys.append(f"kegg:{upper.rsplit(':', 1)[1]}")
    for match in candidate.get("matches", []) or []:
        namespace = str(match.get("namespace") or "").upper()
        raw = str(match.get("raw_value") or "").strip()
        match_field = str(match.get("match_field") or "").casefold()
        if not raw:
            continue
        raw_upper = raw.upper()
        if namespace == "HMDB" or raw_upper.startswith("HMDB"):
            keys.append(f"hmdb:{raw_upper.replace('HMDB:', '')}")
        elif namespace in {"CID", "PUBCHEM"}:
            keys.append(f"pubchem:{raw_upper.replace('PUBCHEM:', '').replace('PUBCHEM.COMPOUND:', '')}")
        elif namespace == "CHEBI" or raw_upper.startswith("CHEBI"):
            keys.append(f"chebi:{raw_upper.replace('CHEBI:', '')}")
        elif "inchikey" in match_field:
            keys.append(f"inchikey:{raw_upper}")
            conn = inchikey_connectivity(raw_upper)
            if conn:
                keys.append(f"inchikey_connectivity:{conn}")
    display_key = normalize_lookup_key(candidate.get("display_name", ""))
    if display_key:
        keys.append(f"name:{display_key}")
    return dedupe_preserve_order(keys)


def preferred_identity_cluster_key(candidate: dict[str, Any]) -> tuple[str, str]:
    keys = candidate_identity_keys(candidate)
    priority = (
        "inchikey_connectivity:",
        "inchikey:",
        "hmdb:",
        "pubchem:",
        "chebi:",
        "kegg:",
        "name:",
    )
    for prefix in priority:
        for key in keys:
            if key.startswith(prefix):
                return prefix.rstrip(":"), key
    return "entity", f"entity:{candidate.get('entity_uid', '')}"


def identity_clusters_for_candidates(candidates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    basis_by_key: dict[str, str] = {}
    for candidate in candidates:
        basis, key = preferred_identity_cluster_key(candidate)
        grouped[key].append(candidate)
        basis_by_key[key] = basis

    clusters: list[dict[str, Any]] = []
    for cluster_key, rows in grouped.items():
        ranked = sorted(
            rows,
            key=lambda row: (
                -float(row.get("score") or 0.0),
                normalize_lookup_key(row.get("display_name", "")),
                row.get("entity_uid", ""),
            ),
        )
        representative = ranked[0]
        clusters.append(
            {
                "contract_version": "identity_cluster.v1",
                "cluster_key": cluster_key,
                "cluster_basis": basis_by_key.get(cluster_key, "entity"),
                "representative_uid": representative.get("entity_uid", ""),
                "representative_name": representative.get("display_name", ""),
                "candidate_uids": [row.get("entity_uid", "") for row in ranked if row.get("entity_uid")],
                "candidate_names": [row.get("display_name", "") for row in ranked if row.get("display_name")],
                "candidate_count": len(ranked),
                "max_score": round(max((float(row.get("score") or 0.0) for row in ranked), default=0.0), 6),
                "allowed_claim_level": "identity_cluster_or_consensus_only",
                "interpretation_boundary": (
                    "Soft identity clusters preserve equivalent or near-equivalent candidates; claims should refer "
                    "to the cluster or shared downstream mechanisms unless manually reviewed."
                ),
                "representative": representative,
            }
        )
    clusters.sort(key=lambda row: (-float(row.get("max_score") or 0.0), row.get("cluster_key", "")))
    return clusters


def european_source_path(workspace: Path, filename: str) -> Path:
    for directory in ("European", "European_point"):
        path = workspace / "raw_lake" / directory / filename
        if path.exists():
            return path
    return workspace / "raw_lake" / "European" / filename


def reported_trait_candidate_names(value: Any) -> list[str]:
    text = clean_reported_trait_name(value)
    if not text:
        return []
    candidates: list[str] = []

    def add(name: str) -> None:
        cleaned = clean_reported_trait_name(name)
        if cleaned and cleaned not in candidates:
            candidates.append(cleaned)
        without_parentheses = clean_reported_trait_name(re.sub(r"\([^)]*\)", "", cleaned))
        if without_parentheses and without_parentheses not in candidates:
            candidates.append(without_parentheses)

    ratio_match = re.search(r"\s+ratio$", text, flags=re.IGNORECASE)
    if ratio_match:
        for part in ratio_trait_component_names(text):
            add(part)
    else:
        add(text)
    return candidates


def parse_chat_json_records(text: str) -> list[dict[str, Any]] | None:
    stripped = text.strip()
    if not stripped.startswith(("[", "{")):
        return None
    parsed = json.loads(stripped)
    records = parsed.get("records") if isinstance(parsed, dict) else parsed
    if not isinstance(records, list):
        raise ValueError("JSON input must be a list or an object with a records list.")
    return [record for index, row in enumerate(records) if (record := normalize_chat_record(row, index))]


def choose_chat_table_dialect(text: str) -> csv.Dialect:
    sample = "\n".join(line for line in text.splitlines()[:20])
    try:
        return csv.Sniffer().sniff(sample, delimiters=",\t;")
    except csv.Error:
        class FallbackDialect(csv.excel):
            delimiter = "\t" if "\t" in sample else "," if "," in sample else ";"

        return FallbackDialect


def parse_chat_table_text(table_text: Any) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    text = str(table_text or "").strip()
    metadata: dict[str, Any] = {
        "source": "table_text",
        "record_count": 0,
        "columns": [],
        "warnings": [],
    }
    if not text:
        metadata["warnings"].append({"code": "empty_table_text", "message": "No table text was provided."})
        return [], metadata
    try:
        json_records = parse_chat_json_records(text)
    except json.JSONDecodeError as exc:
        metadata["warnings"].append({"code": "invalid_json_records", "message": str(exc)})
        return [], metadata
    if json_records is not None:
        metadata["source"] = "json_records"
        metadata["record_count"] = len(json_records)
        metadata["columns"] = sorted({key for record in json_records for key in record})
        return json_records, metadata

    lines = [line for line in text.replace("\r\n", "\n").replace("\r", "\n").split("\n") if line.strip()]
    if not lines:
        metadata["warnings"].append({"code": "empty_table_text", "message": "No non-empty table rows were found."})
        return [], metadata
    dialect = choose_chat_table_dialect("\n".join(lines))
    rows = [row for row in csv.reader(io.StringIO("\n".join(lines)), dialect) if any(str(cell).strip() for cell in row)]
    if not rows:
        metadata["warnings"].append({"code": "empty_table_text", "message": "No non-empty table rows were found."})
        return [], metadata

    header = [canonical_chat_column(cell) for cell in rows[0]]
    header_hits = sum(1 for field in header if field in SUPPORTED_CHAT_RECORD_FIELDS)
    records = []
    if header_hits:
        metadata["columns"] = header
        repaired_rows = 0
        for index, row in enumerate(rows[1:]):
            if len(row) > len(header) and header and header[0] in {"name", "metabolite", "compound"}:
                extra_cells = len(row) - len(header)
                delimiter = getattr(dialect, "delimiter", ",") or ","
                row = [delimiter.join(row[: extra_cells + 1]), *row[extra_cells + 1 :]]
                repaired_rows += 1
            raw_record = {}
            for col_index, cell in enumerate(row):
                if col_index >= len(header):
                    continue
                field = header[col_index]
                if not field:
                    continue
                raw_record[field] = cell
            record = normalize_chat_record(raw_record, index)
            if record:
                records.append(record)
        if repaired_rows:
            metadata["warnings"].append(
                {
                    "code": "row_width_repaired",
                    "message": "Some rows had more cells than the header; extra leading cells were merged into the first metabolite/name column. Quote names containing delimiters for reproducibility.",
                    "count": repaired_rows,
                }
            )
    else:
        if len(rows) == 1 or any(re.search(r"(?:上调|下调|升高|降低|log2fc|p\s*[=:<：≤]|p值)", line, flags=re.IGNORECASE) for line in lines):
            free_records, free_meta = parse_chat_free_text(text)
            free_meta["source"] = "free_text_from_table_box"
            return free_records, free_meta
        metadata["columns"] = ["name", "log2FC", "pvalue"]
        metadata["warnings"].append(
            {
                "code": "header_not_detected",
                "message": "The table had no recognized header; rows were interpreted as name, optional log2FC, optional pvalue.",
            }
        )
        for index, row in enumerate(rows):
            cells = [str(cell).strip() for cell in row if str(cell).strip()]
            if not cells:
                continue
            raw_record = {"name": cells[0]}
            if len(cells) > 1:
                raw_record["log2FC"] = cells[1]
            if len(cells) > 2:
                raw_record["pvalue"] = cells[2]
            record = normalize_chat_record(raw_record, index)
            if record:
                records.append(record)

    metadata["record_count"] = len(records)
    if not records:
        metadata["warnings"].append({"code": "no_records_parsed", "message": "No usable records could be parsed."})
    return records, metadata


def parse_chat_metabolite_list(text: Any) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    raw = str(text or "").strip()
    metadata: dict[str, Any] = {
        "source": "metabolite_list",
        "record_count": 0,
        "columns": ["name"],
        "warnings": [],
    }
    if not raw:
        return [], metadata
    tokens = [token.strip() for token in re.split(r"[\n,;，；、]+", raw) if token.strip()]
    records = [normalize_chat_record({"name": token}, index) for index, token in enumerate(tokens)]
    records = [record for record in records if record]
    metadata["record_count"] = len(records)
    return records, metadata


FREE_TEXT_STOPWORDS = {
    "please",
    "explain",
    "interpret",
    "analysis",
    "result",
    "results",
    "pathway",
    "pathways",
    "target",
    "targets",
    "disease",
    "diseases",
    "evidence",
    "summary",
    "metabolite",
    "metabolites",
    "compound",
    "compounds",
    "help",
    "what",
    "why",
    "how",
    "show",
}


def compact_free_text_metadata(records: list[dict[str, Any]], warnings: list[dict[str, Any]] | None = None) -> dict[str, Any]:
    return {
        "source": "free_text",
        "record_count": len(records),
        "columns": sorted({key for record in records for key in record}),
        "warnings": warnings or [],
    }


def extract_first_number(pattern: str, text: str) -> float | None:
    match = re.search(pattern, text, flags=re.IGNORECASE)
    if not match:
        return None
    return maybe_number(match.group(1))


def strip_free_text_features(text: str) -> str:
    stripped = text
    patterns = [
        r"\b(?:log2fc|log2_fc|logfc|log_fc|log2\s*fold\s*change)\s*[=:：]?\s*[-+]?\d+(?:\.\d+)?(?:e[-+]?\d+)?",
        r"\b(?:fold_change|foldchange|fc)\s*[=:：]?\s*[-+]?\d+(?:\.\d+)?(?:e[-+]?\d+)?",
        r"\b(?:padj|fdr|qvalue|q_value|adjusted\s*p(?:\s*value)?)\s*[=:<：≤]?\s*[-+]?\d+(?:\.\d+)?(?:e[-+]?\d+)?",
        r"\b(?:pvalue|p_value|pval|p)\s*[=:<：≤]?\s*[-+]?\d+(?:\.\d+)?(?:e[-+]?\d+)?",
        r"p\s*值\s*[=:<：≤]?\s*[-+]?\d+(?:\.\d+)?(?:e[-+]?\d+)?",
        r"校正\s*p\s*值\s*[=:<：≤]?\s*[-+]?\d+(?:\.\d+)?(?:e[-+]?\d+)?",
        r"\b(?:mz|m/z)\s*[=:：]?\s*[-+]?\d+(?:\.\d+)?",
        r"\b(?:rt|retention\s*time)\s*[=:：]?\s*[-+]?\d+(?:\.\d+)?",
        r"\b(?:up|down|upregulated|downregulated|increased|decreased|increase|decrease|higher|lower)\b",
        r"(?:上调|下调|升高|降低|增加|减少|不变|无变化)",
    ]
    for pattern in patterns:
        stripped = re.sub(pattern, " ", stripped, flags=re.IGNORECASE)
    stripped = re.sub(r"\bHMDB\d+\b", " ", stripped, flags=re.IGNORECASE)
    stripped = re.sub(r"\bCHEBI[:：]?\s*\d+\b", " ", stripped, flags=re.IGNORECASE)
    stripped = re.sub(r"\b(?:KEGG[:：]?)?\s*[CDG]\d{5}\b", " ", stripped, flags=re.IGNORECASE)
    stripped = re.sub(r"\b(?:PubChem|CID)[:：]?\s*\d+\b", " ", stripped, flags=re.IGNORECASE)
    stripped = re.sub(r"[，,。.;；:：()\[\]{}<>|]+", " ", stripped)
    stripped = re.sub(r"\s+", " ", stripped).strip()
    return stripped


def free_text_name(fragment: str) -> str:
    cleaned = strip_free_text_features(fragment)
    cleaned = re.sub(r"(请|帮我|解释|分析|看看|这个|这些|结果|代谢组|代谢物|通路|靶点|疾病|证据|总结|说明|什么意思)", " ", cleaned)
    cleaned = re.sub(r"\s+", " ", cleaned).strip()
    if not cleaned:
        return ""
    # Prefer short English biochemical names from the fragment; Chinese disease/context words are usually question text.
    tokens = re.findall(r"[A-Za-z][A-Za-z0-9'α-ωΑ-Ω+\-]*", cleaned)
    tokens = [token for token in tokens if token.casefold() not in FREE_TEXT_STOPWORDS and len(token) > 1]
    if tokens:
        return " ".join(tokens[:4])
    if len(cleaned) <= 40 and not re.search(r"[?？]", cleaned):
        return cleaned
    return ""


def free_text_record_from_fragment(fragment: str, row_index: int) -> dict[str, Any]:
    text = fragment.strip()
    record: dict[str, Any] = {"input_id": f"row_{row_index}"}
    hmdb = re.search(r"\b(HMDB\d{5,9})\b", text, flags=re.IGNORECASE)
    if hmdb:
        record["hmdb"] = hmdb.group(1).upper()
    chebi = re.search(r"\bCHEBI[:：]?\s*(\d+)\b", text, flags=re.IGNORECASE)
    if chebi:
        record["chebi"] = f"CHEBI:{chebi.group(1)}"
    pubchem = re.search(r"\b(?:PubChem|CID)[:：]?\s*(\d+)\b", text, flags=re.IGNORECASE)
    if pubchem:
        record["pubchem_cid"] = pubchem.group(1)
    kegg = re.search(r"\b(?:KEGG[:：]?\s*)?([CDG]\d{5})\b", text, flags=re.IGNORECASE)
    if kegg:
        record["kegg"] = kegg.group(1).upper()
    inchikey = re.search(r"\b([A-Z]{14}-[A-Z]{10}-[A-Z])\b", text)
    if inchikey:
        record["inchikey"] = inchikey.group(1)

    name = free_text_name(text)
    if name and not any(key in record for key in ("hmdb", "chebi", "pubchem_cid", "kegg", "inchikey")):
        record["name"] = name

    for field, pattern in (
        ("log2FC", r"(?:log2fc|log2_fc|logfc|log_fc|log2\s*fold\s*change)\s*[=:：]?\s*([-+]?\d+(?:\.\d+)?(?:e[-+]?\d+)?)"),
        ("fold_change", r"(?:fold_change|foldchange|fc)\s*[=:：]?\s*([-+]?\d+(?:\.\d+)?(?:e[-+]?\d+)?)"),
        ("padj", r"(?:padj|fdr|qvalue|q_value|adjusted\s*p(?:\s*value)?|校正\s*p\s*值)\s*[=:<：≤]?\s*([-+]?\d+(?:\.\d+)?(?:e[-+]?\d+)?)"),
        ("pvalue", r"(?:pvalue|p_value|pval|p\s*值|\bp\b)\s*[=:<：≤]?\s*([-+]?\d+(?:\.\d+)?(?:e[-+]?\d+)?)"),
        ("mz", r"(?:mz|m/z)\s*[=:：]?\s*([-+]?\d+(?:\.\d+)?)"),
        ("rt", r"(?:rt|retention\s*time)\s*[=:：]?\s*([-+]?\d+(?:\.\d+)?)"),
    ):
        value = extract_first_number(pattern, text)
        if value is not None:
            record[field] = value

    direction = ""
    lowered = text.casefold()
    for token, normalized in DIRECTION_CHAT_ALIASES.items():
        if token and (token in lowered or token in text):
            direction = normalized
            break
    if direction:
        record["direction"] = direction

    if any(key in record for key in ("name", "hmdb", "chebi", "pubchem_cid", "kegg", "inchikey", "mz")):
        return record
    return {}


def parse_chat_free_text(text: Any) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    raw = str(text or "").strip()
    if not raw:
        return [], compact_free_text_metadata([])
    fragments = [part.strip() for part in re.split(r"[\n;；。]+", raw) if part.strip()]
    if len(fragments) == 1:
        # Commas often separate free-text metabolite clauses when no table header is present.
        fragments = [part.strip() for part in re.split(r"[,，、]+", raw) if part.strip()]
    records = []
    for fragment in fragments:
        record = free_text_record_from_fragment(fragment, len(records))
        if record:
            records.append(record)
    if not records:
        # Last fallback: pull obvious IDs and short English names from the whole text.
        records = [
            record
            for token in re.findall(r"\b(?:HMDB\d{5,9}|CHEBI[:：]?\s*\d+|[A-Za-z][A-Za-z0-9+\-]{2,})\b", raw, flags=re.IGNORECASE)
            if (record := free_text_record_from_fragment(token, len(records)))
        ]
    warnings = []
    if records:
        warnings.append(
            {
                "code": "free_text_extraction_used",
                "message": "Records were extracted from loose text; review matched, ambiguous, and unmatched rows before trusting rankings.",
            }
        )
    return records, compact_free_text_metadata(records, warnings)


def chat_records_from_body(body: dict[str, Any]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    raw_records = body.get("records")
    if isinstance(raw_records, list) and raw_records:
        records = [record for index, row in enumerate(raw_records) if (record := normalize_chat_record(row, index))]
        records, analysis_metadata = prepare_analysis_records(records, {"source": "records"})
        return records, {
            "source": "records",
            "record_count": len(records),
            "columns": sorted({key for record in records for key in record}),
            "warnings": [],
            **{key: value for key, value in analysis_metadata.items() if key not in {"source", "record_count", "columns", "warnings"}},
        }
    table_text = body.get("table_text") or body.get("table") or body.get("csv") or body.get("tsv")
    if table_text:
        records, metadata = parse_chat_table_text(table_text)
        return prepare_analysis_records(records, metadata)
    metabolite_list = body.get("metabolites") or body.get("metabolite_list") or body.get("compound_list")
    if metabolite_list:
        records, metadata = parse_chat_metabolite_list(metabolite_list)
        return prepare_analysis_records(records, metadata)
    free_text = body.get("free_text") or body.get("message") or body.get("question")
    if free_text:
        records, metadata = parse_chat_free_text(free_text)
        return prepare_analysis_records(records, metadata)
    return [], {
        "source": "none",
        "record_count": 0,
        "columns": [],
        "warnings": [{"code": "missing_records", "message": "Upload a table, paste a list, or provide records."}],
        "analysis_mode": "metabolite_table",
        "input_format_detected": "none",
        "converted_record_count": 0,
    }


def chat_question_from_body(body: dict[str, Any]) -> str:
    question = str(body.get("question") or body.get("message") or "").strip()
    if question:
        return question
    return "请解释这个代谢组结果，说明匹配情况、主要通路、靶点/疾病关联和可追溯证据。"


def neutral_mass_from_mz(mz: Any, adduct: Any = "", charge: Any = None) -> tuple[float | None, dict[str, Any]]:
    mz_value = maybe_number(mz)
    if mz_value is None or mz_value <= 0:
        return None, {"reason": "missing_mz"}
    key = normalize_adduct(adduct)
    spec = ADDUCTS.get(key)
    if spec:
        absolute_charge = abs(int(spec["charge"]))
        neutral = mz_value * absolute_charge - float(spec["delta"])
        return neutral, {"mz": mz_value, "adduct": key, "charge": spec["charge"], "neutral_mass": neutral}
    charge_value = maybe_number(charge)
    if charge_value:
        absolute_charge = abs(int(charge_value))
        if absolute_charge > 0:
            return mz_value * absolute_charge, {
                "mz": mz_value,
                "adduct": key,
                "charge": int(charge_value),
                "neutral_mass": mz_value * absolute_charge,
                "warning": "unknown_adduct_used_charge_only",
            }
    return mz_value, {"mz": mz_value, "adduct": key, "charge": charge, "neutral_mass": mz_value, "warning": "unknown_adduct_assumed_neutral"}


@dataclass(frozen=True)
class ServiceConfig:
    min_match_score: float = DEFAULT_MIN_MATCH_SCORE
    min_match_margin: float = DEFAULT_MIN_MATCH_MARGIN
    max_candidates: int = DEFAULT_MAX_CANDIDATES
    max_paths: int = DEFAULT_MAX_PATHS
    max_hops: int = DEFAULT_MAX_HOPS
    max_edges_per_node: int = DEFAULT_MAX_EDGES_PER_NODE
    propagation_max_adaptive_edges_per_node: int = DEFAULT_PROPAGATION_MAX_ADAPTIVE_EDGES_PER_NODE
    max_subgraph_edges: int = DEFAULT_MAX_SUBGRAPH_EDGES
    max_evidence_items: int = DEFAULT_MAX_EVIDENCE_ITEMS
    default_ppm_tolerance: float = DEFAULT_PPM_TOLERANCE
    default_rt_tolerance: float = DEFAULT_RT_TOLERANCE
    default_ms2_mz_tolerance: float = DEFAULT_MS2_MZ_TOLERANCE
    propagation_alpha: float = DEFAULT_PROPAGATION_ALPHA
    propagation_iterations: int = DEFAULT_PROPAGATION_ITERATIONS
    propagation_tolerance: float = DEFAULT_PROPAGATION_TOLERANCE
    max_propagation_nodes: int = DEFAULT_MAX_PROPAGATION_NODES
    max_propagation_edges: int = DEFAULT_MAX_PROPAGATION_EDGES
    propagation_beam_per_type: int = DEFAULT_PROPAGATION_BEAM_PER_TYPE
    propagation_retained_mass: float = DEFAULT_PROPAGATION_RETAINED_MASS
    propagation_degree_penalty: bool = True
    min_literature_overlay_prob: float = DEFAULT_MIN_LITERATURE_OVERLAY_PROB
    max_prediction_rows: int = DEFAULT_MAX_PREDICTION_ROWS

    def as_dict(self) -> dict[str, Any]:
        return {
            "min_match_score": self.min_match_score,
            "min_match_margin": self.min_match_margin,
            "max_candidates": self.max_candidates,
            "max_paths": self.max_paths,
            "max_hops": self.max_hops,
            "max_edges_per_node": self.max_edges_per_node,
            "propagation_max_adaptive_edges_per_node": self.propagation_max_adaptive_edges_per_node,
            "max_subgraph_edges": self.max_subgraph_edges,
            "max_evidence_items": self.max_evidence_items,
            "default_ppm_tolerance": self.default_ppm_tolerance,
            "default_rt_tolerance": self.default_rt_tolerance,
            "default_ms2_mz_tolerance": self.default_ms2_mz_tolerance,
            "propagation_alpha": self.propagation_alpha,
            "propagation_iterations": self.propagation_iterations,
            "propagation_tolerance": self.propagation_tolerance,
            "max_propagation_nodes": self.max_propagation_nodes,
            "max_propagation_edges": self.max_propagation_edges,
            "propagation_beam_per_type": self.propagation_beam_per_type,
            "propagation_retained_mass": self.propagation_retained_mass,
            "propagation_degree_penalty": self.propagation_degree_penalty,
            "min_literature_overlay_prob": self.min_literature_overlay_prob,
            "max_prediction_rows": self.max_prediction_rows,
        }


class MetaboService:
    def __init__(
        self,
        workspace: Path,
        normalized_root: Path | None = None,
        graph_root: Path | None = None,
        pubchem_root: Path | None = None,
        compound_root: Path | None = None,
        literature_root: Path | None = None,
        prediction_root: Path | None = None,
        release_id: str | None = None,
        config: ServiceConfig | None = None,
        llm_config: ExternalLLMConfig | None = None,
        llm_transport: Any | None = None,
    ):
        require_arrow()
        self.workspace = workspace.resolve()
        self.normalized_root = (normalized_root or self.workspace / DEFAULT_NORMALIZED_ROOT).resolve()
        self.graph_root = (graph_root or self.workspace / DEFAULT_GRAPH_ROOT).resolve()
        self.pubchem_root = (pubchem_root or self.workspace / DEFAULT_PUBCHEM_ROOT).resolve()
        self.compound_root = (compound_root or self.workspace / DEFAULT_COMPOUND_MATCH_ROOT).resolve()
        self.literature_root = (literature_root or self.workspace / DEFAULT_LITERATURE_ROOT).resolve()
        self.prediction_root = (prediction_root or self.workspace / DEFAULT_PREDICTION_OVERLAY_ROOT).resolve()
        self.database_accuracy_root = (self.workspace / DEFAULT_DATABASE_ACCURACY_ROOT).resolve()
        self.release_id = release_id or latest_release_id(self.graph_root)
        self.config = config or ServiceConfig()
        self.llm_config = llm_config or ExternalLLMConfig.from_env()
        self.llm_transport = llm_transport
        self.normalized_dir = self.normalized_root / self.release_id
        self.graph_dir = self.graph_root / self.release_id
        self.pubchem_dir = self.pubchem_root / self.release_id
        self.compound_dir = self.compound_root / self.release_id
        self.literature_dir = self.literature_root / self.release_id
        self.prediction_dir = self.prediction_root / self.release_id
        self.database_accuracy_dir = self.database_accuracy_root / self.release_id
        self.database_accuracy_manifest_path = self.database_accuracy_dir / "database_accuracy_store_manifest.json"
        self._node_cache: dict[str, dict[str, Any] | None] = {}
        self._node_idx_cache: dict[int, dict[str, Any] | None] = {}
        self._incident_edge_cache: dict[str, list[dict[str, Any]]] = {}
        self._edge_type_index_cache: dict[int, dict[str, Any]] | None = None
        self._sparse_incident_cache: dict[int, list[dict[str, Any]]] = {}
        self._pathway_index: dict[str, Any] | None = None
        self._literature_support_rows: list[dict[str, Any]] | None = None
        self._literature_support_by_edge_uid: dict[str, list[dict[str, Any]]] | None = None
        self._literature_support_by_triple: dict[tuple[str, str, str], list[dict[str, Any]]] | None = None
        self._literature_support_by_entity: dict[str, list[dict[str, Any]]] | None = None
        self._prediction_overlay_cache: dict[str, dict[str, Any]] = {}
        self._european_trait_index: dict[str, dict[str, Any]] | None = None
        self._european_trait_annotation_index: dict[str, dict[str, Any]] | None = None
        self._compound_identifier_query_cache: dict[tuple[str, str], list[dict[str, Any]]] = {}
        self._compound_name_query_cache: dict[str, list[dict[str, Any]]] = {}
        self._compound_formula_query_cache: dict[str, list[dict[str, Any]]] = {}
        self._compound_inchikey_full_cache: dict[str, list[dict[str, Any]]] = {}
        self._compound_inchikey_connectivity_cache: dict[str, list[dict[str, Any]]] = {}
        self._compound_mass_range_cache: dict[tuple[float, float], list[dict[str, Any]]] = {}
        self._compound_rt_candidate_cache: dict[tuple[str, ...], list[dict[str, Any]]] = {}
        self._compound_ms2_candidate_cache: dict[tuple[str, ...], list[dict[str, Any]]] = {}
        self._metabolite_graph_support_cache: dict[str, tuple[float, dict[str, Any]]] = {}
        self._biological_entity_pool_cache: dict[tuple[str, ...], list[dict[str, Any]]] = {}
        self._database_accuracy_tables: dict[str, list[dict[str, Any]]] = {}
        self._database_accuracy_indexes: dict[str, dict[str, Any]] = {}
        self._gold_standard_cache: dict[str, Any] | None = None

    def database_accuracy_store_available(self) -> bool:
        return self.database_accuracy_manifest_path.exists()

    def database_accuracy_manifest(self) -> dict[str, Any]:
        return read_json_file(self.database_accuracy_manifest_path) if self.database_accuracy_store_available() else {}

    def database_accuracy_table_path(self, section: str, table_name: str) -> Path:
        return self.database_accuracy_dir / section / f"{table_name}.parquet"

    def database_accuracy_rows(self, section: str, table_name: str, limit: int | None = None) -> list[dict[str, Any]]:
        key = f"{section}/{table_name}"
        if key not in self._database_accuracy_tables:
            path = self.database_accuracy_table_path(section, table_name)
            self._database_accuracy_tables[key] = table_rows(path) if path.exists() else []
        rows = self._database_accuracy_tables[key]
        return rows[:limit] if limit is not None else rows

    def database_accuracy_index(self, section: str, table_name: str, field: str) -> dict[str, Any]:
        key = f"{section}/{table_name}/{field}"
        if key not in self._database_accuracy_indexes:
            self._database_accuracy_indexes[key] = {
                str(row.get(field) or ""): row
                for row in self.database_accuracy_rows(section, table_name)
                if row.get(field) not in {None, ""}
            }
        return self._database_accuracy_indexes[key]

    def database_accuracy_multindex(self, section: str, table_name: str, field: str) -> dict[str, list[dict[str, Any]]]:
        key = f"{section}/{table_name}/{field}/multi"
        if key not in self._database_accuracy_indexes:
            grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
            for row in self.database_accuracy_rows(section, table_name):
                value = str(row.get(field) or "")
                if value:
                    grouped[value].append(row)
            self._database_accuracy_indexes[key] = dict(grouped)
        return self._database_accuracy_indexes[key]  # type: ignore[return-value]

    def database_accuracy_trait_row_for_record(self, record: Any) -> dict[str, Any]:
        if not isinstance(record, dict):
            return {}
        normalized = self.normalized_record_keys(record)
        accession = normalize_gwas_accession(
            first_present(normalized, ("accession_id", "accessionid", "trait", "gcst", "gwas_trait"))
        )
        if not accession:
            return {}
        return self.database_accuracy_index("entity_store", "trait_entities", "accession_id").get(accession, {})

    def database_accuracy_class_row_for_label(self, label: str) -> dict[str, Any]:
        key = normalize_lookup_key(label)
        if not key:
            return {}
        for row in self.database_accuracy_rows("entity_store", "metabolite_classes"):
            if normalize_lookup_key(row.get("class_name")) == key:
                return row
        return {}

    def database_accuracy_input_label(self, row_contract: dict[str, Any]) -> str:
        record = row_contract.get("record", {})
        normalized = self.normalized_record_keys(record)
        annotation = self.european_annotation_for_record(record)
        values = [
            first_present(normalized, ("reported_trait", "reportedtrait", "name", "metabolite", "compound", "trait", "accession_id", "gcst")),
            annotation.get("reported_trait", ""),
            row_contract.get("query", ""),
        ]
        for value in values:
            text = str(value or "").strip()
            if text:
                return text
        return str(row_contract.get("input_id") or "").strip()

    def database_accuracy_feature_type(self, row_contract: dict[str, Any], label: str) -> str:
        record = row_contract.get("record", {})
        normalized = self.normalized_record_keys(record)
        has_trait = any(
            normalized.get(field) not in {None, ""}
            for field in ("trait", "accession_id", "accessionid", "gcst", "gwas_trait", "reported_trait", "reportedtrait")
        )
        if has_trait:
            combined = " ".join(
                str(value or "")
                for value in (
                    label,
                    row_contract.get("query", ""),
                    self.european_annotation_for_record(record).get("mapped_names", ""),
                    self.european_annotation_for_record(record).get("candidate_names", ""),
                )
                if value
            )
            if trait_text_is_ratio_or_composite(combined):
                return "ratio"
            if trait_text_is_pool_or_class(combined):
                return "lipid_class"
            return "trait_score"
        if trait_text_is_pool_or_class(label):
            return "lipid_class"
        if re.match(r"^X[-_ ]?\d+$", str(label or "").strip(), flags=re.IGNORECASE):
            return "unknown_peak"
        return "direct_metabolite"

    def database_accuracy_input_feature(self, row_contract: dict[str, Any]) -> dict[str, Any]:
        record = row_contract.get("record", {})
        feature = analysis_feature_from_record(record, row_contract.get("input_id", ""))
        label = self.database_accuracy_input_label(row_contract)
        feature_type = self.database_accuracy_feature_type(row_contract, label)
        normalized = self.normalized_record_keys(record)
        return {
            "input_feature_uid": f"input_feature_{content_hash({'input_id': row_contract.get('input_id', ''), 'label': label, 'type': feature_type})[:20]}",
            "input_row_id": row_contract.get("input_id", ""),
            "feature_label": label,
            "feature_type": feature_type,
            "effect_value": feature.get("log2_fold_change"),
            "effect_label": "log2FC" if feature.get("log2_fold_change") is not None else "",
            "direction": feature.get("direction", "unknown"),
            "pvalue": feature.get("p_value"),
            "padj": feature.get("adjusted_p_value"),
            "sample_context_uid": str(normalized.get("sample_context_uid") or normalized.get("context_uid") or ""),
            "raw_input": record,
        }

    def database_accuracy_identity_candidates(
        self,
        row_contract: dict[str, Any],
        input_feature: dict[str, Any],
    ) -> list[dict[str, Any]]:
        feature_type = str(input_feature.get("feature_type") or "")
        label = str(input_feature.get("feature_label") or "")
        candidates: list[dict[str, Any]] = []
        trait_row = self.database_accuracy_trait_row_for_record(row_contract.get("record", {}))
        if trait_row:
            candidates.append(
                {
                    "candidate_uid": f"candidate_{content_hash({'feature': input_feature.get('input_feature_uid'), 'trait': trait_row.get('trait_uid')})[:20]}",
                    "input_feature_uid": input_feature.get("input_feature_uid", ""),
                    "candidate_entity_uid": trait_row.get("trait_uid", ""),
                    "candidate_entity_type": "trait",
                    "candidate_name": trait_row.get("reported_trait", "") or label,
                    "match_basis": "trait_annotation",
                    "match_score": 100.0,
                    "match_margin": 100.0,
                    "candidate_rank": 1,
                    "blocking_reason": "ratio_not_exact_abundance" if feature_type == "ratio" else "",
                }
            )
        if feature_type == "lipid_class":
            class_row = self.database_accuracy_class_row_for_label(label)
            class_uid = class_row.get("class_uid", "") or f"class_{content_hash(label)[:16]}"
            candidates.append(
                {
                    "candidate_uid": f"candidate_{content_hash({'feature': input_feature.get('input_feature_uid'), 'class': class_uid})[:20]}",
                    "input_feature_uid": input_feature.get("input_feature_uid", ""),
                    "candidate_entity_uid": class_uid,
                    "candidate_entity_type": "class",
                    "candidate_name": class_row.get("class_name", "") or label,
                    "match_basis": "class_name",
                    "match_score": 70.0 if class_row else 50.0,
                    "match_margin": 0.0,
                    "candidate_rank": len(candidates) + 1,
                    "blocking_reason": "class_name_not_exact_identity",
                }
            )
        for rank, candidate in enumerate((row_contract.get("resolution", {}) or {}).get("candidates") or [], start=1):
            score = float(candidate.get("score") or 0.0)
            score_components = candidate.get("score_components") or {}
            blocking_reason = ""
            if score_components.get("common_biochemical_zero_score_rescue"):
                blocking_reason = ""
            elif score <= 0.0:
                blocking_reason = "zero_score_candidate"
            elif feature_type == "ratio":
                blocking_reason = "ratio_component_not_abundance_seed"
            elif feature_type == "lipid_class":
                blocking_reason = "class_name_not_exact_identity"
            elif "high_risk_name" in (row_contract.get("resolution", {}) or {}).get("identity_review_reasons", []):
                blocking_reason = "high_risk_name"
            candidates.append(
                {
                    "candidate_uid": f"candidate_{content_hash({'feature': input_feature.get('input_feature_uid'), 'chemical': candidate.get('entity_uid', ''), 'rank': rank})[:20]}",
                    "input_feature_uid": input_feature.get("input_feature_uid", ""),
                    "candidate_entity_uid": candidate.get("entity_uid", ""),
                    "candidate_entity_type": "chemical",
                    "candidate_name": candidate.get("display_name", "") or candidate.get("entity_uid", ""),
                    "match_basis": candidate.get("match_field", "") or "compound_resolver",
                    "match_score": round(score, 6),
                    "match_margin": round(float(row_contract.get("resolution", {}).get("top_margin") or 0.0), 6),
                    "candidate_rank": rank + (1 if trait_row else 0),
                    "blocking_reason": blocking_reason,
                }
            )
        return candidates

    def database_accuracy_identity_decision(
        self,
        row_contract: dict[str, Any],
        input_feature: dict[str, Any],
        candidates: list[dict[str, Any]],
    ) -> dict[str, Any]:
        feature_type = str(input_feature.get("feature_type") or "")
        resolution = row_contract.get("resolution", {}) or {}
        status = str(resolution.get("status") or "unmatched")
        chemical_candidates = [row for row in candidates if row.get("candidate_entity_type") == "chemical"]
        trait_candidates = [row for row in candidates if row.get("candidate_entity_type") == "trait"]
        class_candidates = [row for row in candidates if row.get("candidate_entity_type") == "class"]
        rejected_reasons = sorted({row.get("blocking_reason", "") for row in candidates if row.get("blocking_reason")})
        top_chemical = chemical_candidates[0] if chemical_candidates else {}
        top_chemical_blocking_reason = str(top_chemical.get("blocking_reason") or "")
        decision_status = "unmatched"
        accepted_uid = ""
        accepted_type = ""
        decision_rule = "no_candidate"
        if feature_type == "ratio":
            if trait_candidates:
                decision_status = "accepted_trait"
                accepted_uid = str(trait_candidates[0].get("candidate_entity_uid") or "")
                accepted_type = "trait"
                decision_rule = "ratio_trait_preserved_not_abundance_seed"
            else:
                decision_status = "ambiguous"
                decision_rule = "ratio_trait_missing_trait_entity"
                rejected_reasons.append("ratio_not_exact_abundance")
        elif feature_type == "lipid_class":
            if class_candidates:
                decision_status = "accepted_class"
                accepted_uid = str(class_candidates[0].get("candidate_entity_uid") or "")
                accepted_type = "class"
                decision_rule = "class_identity_preserved_not_exact_compound"
            else:
                decision_status = "ambiguous"
                decision_rule = "class_name_not_exact_identity"
                rejected_reasons.append("class_name_not_exact_identity")
        elif chemical_candidates and float(top_chemical.get("match_score") or 0.0) <= 0.0:
            decision_status = "rejected"
            decision_rule = "zero_score_candidate"
            rejected_reasons.append("zero_score_candidate")
        elif status == "matched" and chemical_candidates and not top_chemical_blocking_reason:
            decision_status = "accepted_exact"
            accepted_uid = str(top_chemical.get("candidate_entity_uid") or "")
            accepted_type = "chemical"
            decision_rule = "stable_identifier_or_unique_high_quality_name"
            rejected_reasons = []
        elif status == "ambiguous":
            decision_status = "ambiguous"
            decision_rule = "resolver_ambiguous"
        elif status == "invalid":
            decision_status = "rejected"
            decision_rule = str(row_contract.get("invalid_reason") or "invalid_input")
        else:
            decision_status = "unmatched"
            decision_rule = "no_accepted_identity"
        return {
            "decision_uid": f"decision_{content_hash({'feature': input_feature.get('input_feature_uid'), 'status': decision_status, 'accepted': accepted_uid, 'rule': decision_rule})[:20]}",
            "input_feature_uid": input_feature.get("input_feature_uid", ""),
            "decision_status": decision_status,
            "accepted_entity_uid": accepted_uid,
            "accepted_entity_type": accepted_type,
            "decision_rule": decision_rule,
            "rejected_reasons": sorted(set(reason for reason in rejected_reasons if reason)),
            "requires_review": decision_status in {"ambiguous", "unmatched", "rejected"},
            "source_record_uid": "",
        }

    def identity_resolution_v2_for_row(self, row_contract: dict[str, Any]) -> dict[str, Any]:
        input_feature = self.database_accuracy_input_feature(row_contract)
        candidates = self.database_accuracy_identity_candidates(row_contract, input_feature)
        decision = self.database_accuracy_identity_decision(row_contract, input_feature, candidates)
        return {
            "contract_version": DATABASE_ACCURACY_CONTRACT_VERSION,
            "store_available": self.database_accuracy_store_available(),
            "input_feature": input_feature,
            "candidates": candidates,
            "decision": decision,
        }

    def identity_resolution_v2_summary(self, rows: list[dict[str, Any]]) -> dict[str, Any]:
        feature_types = Counter()
        decision_statuses = Counter()
        rejected_reasons = Counter()
        for row in rows:
            identity = row.get("identity_resolution_v2") or {}
            feature_types[str((identity.get("input_feature") or {}).get("feature_type") or "unknown")] += 1
            decision = identity.get("decision") or {}
            decision_statuses[str(decision.get("decision_status") or "unknown")] += 1
            for reason in decision.get("rejected_reasons") or []:
                rejected_reasons[str(reason)] += 1
        return {
            "contract_version": DATABASE_ACCURACY_CONTRACT_VERSION,
            "store_available": self.database_accuracy_store_available(),
            "feature_types": dict(sorted(feature_types.items())),
            "decision_statuses": dict(sorted(decision_statuses.items())),
            "rejected_reasons": dict(sorted(rejected_reasons.items())),
        }

    def cache_rows_by_field(
        self,
        cache: dict[Any, list[dict[str, Any]]],
        path: Path,
        field_name: str,
        keys: Iterable[Any],
        extra_filter: ds.Expression | None = None,
        chunk_size: int = 500,
    ) -> None:
        if not path.exists():
            return
        wanted = [key for key in dict.fromkeys(keys) if key not in {None, ""} and key not in cache]
        if not wanted:
            return
        for key in wanted:
            cache[key] = []
        for offset in range(0, len(wanted), chunk_size):
            chunk = wanted[offset : offset + chunk_size]
            filt = ds.field(field_name).isin(chunk)
            if extra_filter is not None:
                filt = extra_filter & filt
            for row in table_rows(path, None, filt):
                cache_key = row.get(field_name)
                if cache_key in cache:
                    cache[cache_key].append(row)

    def prime_compound_resolver_caches(self, records: list[Any]) -> None:
        if not self.compound_index_available() or not records:
            return
        identifier_keys: dict[str, set[str]] = defaultdict(set)
        name_keys: set[str] = set()
        formula_keys: set[str] = set()
        inchikey_keys: set[str] = set()
        connectivity_keys: set[str] = set()

        for record in records:
            parsed = self.compound_record_queries(record)
            for identifier in parsed.get("identifiers", []):
                namespace = str(identifier.get("namespace") or "").upper()
                value = str(identifier.get("value") or "").strip()
                if not namespace or not value:
                    continue
                if namespace == "CID":
                    value = normalize_id_token(value).removeprefix("pubchem:").removeprefix("cid:")
                elif ":" in value:
                    value = value.split(":", 1)[1]
                lookup_key = normalize_id_token(value)
                if lookup_key:
                    identifier_keys[namespace].add(lookup_key)
            for name, _mode in self.expanded_compound_name_queries(parsed.get("names", [])):
                key = normalize_lookup_key(name)
                if key:
                    name_keys.add(key)
            formula_key = normalize_formula(parsed.get("formula", ""))
            if formula_key:
                formula_keys.add(formula_key)
            inchikey = str(parsed.get("inchikey") or "").strip()
            if inchikey:
                key = normalize_id_token(inchikey).upper()
                if key:
                    inchikey_keys.add(key)
                    connectivity_keys.add(inchikey_connectivity(key))

        identifier_path = self.compound_dir / "compound_identifier_index.parquet"
        for namespace, keys in identifier_keys.items():
            missing = [
                key
                for key in keys
                if (namespace, key) not in self._compound_identifier_query_cache
            ]
            for key in missing:
                self._compound_identifier_query_cache[(namespace, key)] = []
            for offset in range(0, len(missing), 500):
                chunk = missing[offset : offset + 500]
                if not chunk or not identifier_path.exists():
                    continue
                filt = (ds.field("namespace") == namespace) & ds.field("lookup_key").isin(chunk)
                for row in table_rows(identifier_path, None, filt):
                    cache_key = (namespace, row.get("lookup_key"))
                    self._compound_identifier_query_cache.setdefault(cache_key, []).append(row)

        self.cache_rows_by_field(
            self._compound_name_query_cache,
            self.compound_dir / "compound_name_index.parquet",
            "name_key",
            name_keys,
        )
        self.cache_rows_by_field(
            self._compound_formula_query_cache,
            self.compound_dir / "compound_formula_mass_index.parquet",
            "formula_key",
            formula_keys,
        )
        self.cache_rows_by_field(
            self._compound_inchikey_full_cache,
            self.compound_dir / "compound_inchikey_index.parquet",
            "inchikey",
            inchikey_keys,
        )
        self.cache_rows_by_field(
            self._compound_inchikey_connectivity_cache,
            self.compound_dir / "compound_inchikey_index.parquet",
            "connectivity_key",
            connectivity_keys,
        )
        candidate_uids: set[str] = set()
        for key in name_keys:
            candidate_uids.update(str(row.get("metabolite_uid") or "") for row in self._compound_name_query_cache.get(key, []))
        for namespace, keys in identifier_keys.items():
            for key in keys:
                candidate_uids.update(
                    str(row.get("metabolite_uid") or "")
                    for row in self._compound_identifier_query_cache.get((namespace, key), [])
                )
        for key in formula_keys:
            candidate_uids.update(str(row.get("metabolite_uid") or "") for row in self._compound_formula_query_cache.get(key, []))
        for key in inchikey_keys:
            candidate_uids.update(str(row.get("metabolite_uid") or "") for row in self._compound_inchikey_full_cache.get(key, []))
        for key in connectivity_keys:
            candidate_uids.update(str(row.get("metabolite_uid") or "") for row in self._compound_inchikey_connectivity_cache.get(key, []))
        self.prime_metabolite_graph_support(uid for uid in candidate_uids if uid)

    def prime_metabolite_graph_support(self, metabolite_uids: Iterable[str]) -> None:
        missing = sorted({uid for uid in metabolite_uids if uid and uid not in self._metabolite_graph_support_cache})
        if not missing:
            return
        edge_counts: dict[str, int] = {uid: 0 for uid in missing}
        edge_path = self.graph_dir / "edges.parquet"
        for offset in range(0, len(missing), 500):
            chunk = missing[offset : offset + 500]
            if not chunk or not edge_path.exists():
                continue
            filt = ds.field("subject_uid").isin(chunk) | ds.field("object_uid").isin(chunk)
            for row in table_rows(edge_path, ["subject_uid", "object_uid"], filt):
                subject_uid = str(row.get("subject_uid") or "")
                object_uid = str(row.get("object_uid") or "")
                if subject_uid in edge_counts:
                    edge_counts[subject_uid] += 1
                if object_uid in edge_counts and object_uid != subject_uid:
                    edge_counts[object_uid] += 1

        self.load_literature_support()
        self.get_nodes(missing)
        for metabolite_uid in missing:
            support_rows = [
                row
                for row in (self._literature_support_by_entity or {}).get(metabolite_uid, [])
                if row.get("support_class") in {"confirm", "support_direction", "novel_candidate"}
                and float(row.get("p_literature") or 0.0) >= self.config.min_literature_overlay_prob
            ]
            max_p_literature = max((float(row.get("p_literature") or 0.0) for row in support_rows), default=0.0)
            score = 0.0
            edge_count = edge_counts.get(metabolite_uid, 0)
            if edge_count:
                score += 5.0 + min(13.0, math.log2(edge_count + 1.0) * 3.0)
            node = self.get_node(metabolite_uid) or {}
            xrefs = parse_list(node.get("external_xrefs"))
            if any(xref.startswith("HMDB:") for xref in xrefs):
                score += 2.0
            if any(xref.startswith("KEGG.COMPOUND:") for xref in xrefs):
                score += 1.0
            if support_rows:
                score += min(4.0, math.log2(len(support_rows) + 1.0))
                score += min(1.5, max_p_literature * 1.5)
            self._metabolite_graph_support_cache[metabolite_uid] = (
                round(min(18.0, score), 6),
                {
                    "graph_edge_count": edge_count,
                    "literature_support_count": len(support_rows),
                    "max_p_literature": round(max_p_literature, 6),
                },
            )

    def release_meta(self) -> dict[str, Any]:
        return {
            "release_id": self.release_id,
            "normalized_manifest": str(self.normalized_dir / "normalized_manifest.json"),
            "graph_manifest": str(self.graph_dir / "graph_manifest.json"),
            "pubchem_manifest": str(self.pubchem_dir / "pubchem_cid_cache_manifest.json"),
            "compound_match_manifest": str(self.compound_dir / "compound_match_index_manifest.json"),
            "literature_evidence_manifest": str(self.literature_dir / "literature_evidence_manifest.json"),
            "prediction_overlay_dir": str(self.prediction_dir),
            "config_hash": content_hash(self.config.as_dict())[:16],
        }

    def response_envelope(self, endpoint: str, request: dict[str, Any], payload: dict[str, Any]) -> dict[str, Any]:
        body = {
            "api_version": "mvp.v1",
            "endpoint": endpoint,
            "release": self.release_meta(),
            "determinism": {
                "request_hash": content_hash(request),
                "config": self.config.as_dict(),
            },
            **payload,
        }
        body["determinism"]["response_hash"] = content_hash({k: v for k, v in body.items() if k != "determinism"})
        return body

    def releases(self) -> dict[str, Any]:
        release_ids = sorted(
            {
                path.name
                for root in (self.normalized_root, self.graph_root, self.pubchem_root, self.compound_root)
                if root.exists()
                for path in root.iterdir()
                if path.is_dir() and path.name.startswith("mvp_")
            }
        )
        releases = []
        for release_id in release_ids:
            graph_manifest = read_json_file(self.graph_root / release_id / "graph_manifest.json")
            normalized_manifest = read_json_file(self.normalized_root / release_id / "normalized_manifest.json")
            pubchem_manifest = read_json_file(self.pubchem_root / release_id / "pubchem_cid_cache_manifest.json")
            compound_manifest = read_json_file(self.compound_root / release_id / "compound_match_index_manifest.json")
            literature_manifest = read_json_file(self.literature_root / release_id / "literature_evidence_manifest.json")
            releases.append(
                {
                    "release_id": release_id,
                    "active": release_id == self.release_id,
                    "has_normalized_store": (self.normalized_root / release_id).exists(),
                    "has_graph_projection": (self.graph_root / release_id).exists(),
                    "has_pubchem_cid_cache": (self.pubchem_root / release_id).exists(),
                    "has_compound_match_index": (self.compound_root / release_id).exists(),
                    "has_literature_evidence": (self.literature_root / release_id).exists(),
                    "node_count": graph_manifest.get("node_count"),
                    "edge_count": next(
                        (row.get("rows") for row in graph_manifest.get("tables", []) if row.get("table") == "edges"),
                        None,
                    ),
                    "normalized_tables": [
                        {"table": row.get("table"), "rows": row.get("rows")}
                        for row in normalized_manifest.get("tables", [])
                    ],
                    "pubchem_target_cid_count": pubchem_manifest.get("target_cid_count"),
                    "pubchem_matched_cid_count": pubchem_manifest.get("matched_cid_count"),
                    "compound_match_tables": [
                        {"table": row.get("table"), "rows": row.get("rows")}
                        for row in compound_manifest.get("tables", [])
                    ],
                    "literature_evidence_tables": [
                        {"table": row.get("table"), "rows": row.get("rows")}
                        for row in literature_manifest.get("tables", [])
                    ],
                    "literature_evidence_metrics": literature_manifest.get("metrics", {}),
                }
            )
        return self.response_envelope("/releases", {}, {"releases": releases})

    def prediction_overlay_rows(self, overlay_type: str) -> dict[str, Any]:
        if overlay_type in self._prediction_overlay_cache:
            return self._prediction_overlay_cache[overlay_type]
        names = PREDICTION_OVERLAY_FILES.get(overlay_type, (overlay_type,))
        rows: list[dict[str, Any]] = []
        files: list[str] = []
        for base_dir in (self.prediction_dir, self.prediction_root):
            if not base_dir.exists():
                continue
            for name in names:
                for suffix in (".json", ".jsonl", ".csv", ".tsv"):
                    path = base_dir / f"{name}{suffix}"
                    if not path.exists():
                        continue
                    loaded = load_records_file(path)
                    for index, row in enumerate(loaded):
                        enriched = dict(row)
                        enriched.setdefault("_overlay_type", overlay_type)
                        enriched.setdefault("_source_file", str(path))
                        enriched.setdefault("_source_row", index + 1)
                        rows.append(enriched)
                    files.append(str(path))
        rows.sort(key=lambda row: (str(row.get("_source_file", "")), int(row.get("_source_row", 0) or 0)))
        payload = {"rows": rows, "files": sorted(set(files)), "row_count": len(rows)}
        self._prediction_overlay_cache[overlay_type] = payload
        return payload

    def normalize_prediction_context(self, context: Any = None, records: list[Any] | None = None) -> dict[str, Any]:
        fields: dict[str, list[str]] = defaultdict(list)
        allowed_fields = {
            "cancer_type",
            "tissue",
            "cell_type",
            "cell_state",
            "cell_line",
            "context_terms",
            "experimental_context",
            "species",
            "platform",
            "sample_type",
        }

        def add_field(field: str, value: Any) -> None:
            canonical = CONTEXT_FIELD_ALIASES.get(str(field or "").strip().casefold().replace(" ", "_"), str(field or ""))
            if canonical not in allowed_fields:
                return
            for term in prediction_terms_from_value(value):
                if term not in fields[canonical]:
                    fields[canonical].append(term)

        if isinstance(context, str):
            stripped = context.strip()
            if stripped:
                try:
                    parsed = json.loads(stripped)
                except json.JSONDecodeError:
                    parsed = {"context_terms": stripped}
                context = parsed
        if isinstance(context, dict):
            for key, value in context.items():
                add_field(key, value)
        elif isinstance(context, (list, tuple, set)):
            add_field("context_terms", list(context))
        elif context:
            add_field("context_terms", context)

        for record in records or []:
            if not isinstance(record, dict):
                continue
            for key, value in record.items():
                canonical = CONTEXT_FIELD_ALIASES.get(str(key or "").strip().casefold().replace(" ", "_"))
                if canonical:
                    add_field(canonical, value)

        ordered_fields = {key: sorted(value) for key, value in sorted(fields.items())}
        all_terms = sorted({term for values in ordered_fields.values() for term in values})
        encoded = self.context_tags_for_terms(ordered_fields, all_terms)
        return {
            "fields": ordered_fields,
            "terms": all_terms,
            "has_context": bool(all_terms),
            "mode": "context_aware" if all_terms else "generalized",
            "context_tags": encoded["tags"],
            "context_vector": encoded["vector"],
        }

    def context_tags_for_terms(self, fields: dict[str, list[str]], terms: list[str]) -> dict[str, Any]:
        joined = " ".join(terms)
        tags: dict[str, float] = {}

        def add(tag: str, weight: float = 1.0) -> None:
            tags[tag] = round(max(tags.get(tag, 0.0), weight), 6)

        if any(term in joined for term in ("kidney", "renal", "nephron", "tubular", "proximal tubule", "proximal tubular")):
            add("renal_organic_anion_uremic_toxin_handling", 1.0)
            add("renal_osmolyte_homeostasis", 0.9)
            add("fatty_acid_oxidation_pressure", 0.75)
            add("proximal_tubule_solute_transport_pressure", 1.0 if "proximal" in joined or "tubular" in joined else 0.8)
        if any(term in joined for term in ("brain", "glioma", "glioblastoma", "astrocy", "glial", "neural", "neuronal")):
            add("glial_glutamate_gaba_metabolism", 1.0)
            add("idh_like_metabolic_pressure", 0.9 if any(term in joined for term in ("glioma", "glioblastoma", "idh")) else 0.65)
            add("hypoxia_glycolytic_shift", 0.7)
        if any(term in joined for term in EPITHELIAL_CONTEXT_TERMS):
            add("epithelial_barrier_membrane_repair", 1.0)
            add("nucleotide_sugar_glycosylation", 0.8)
            add("phospholipid_membrane_remodeling", 0.8)
        if any(term in joined for term in ("hypoxia", "low oxygen", "ischemia", "glycolytic")):
            add("hypoxia_glycolytic_shift", 1.0)
        vector = {theme_id: 0.0 for theme_id in sorted(BIOCHEMICAL_THEME_SPECS)}
        for tag, tag_weight in tags.items():
            for theme_id, boost in CONTEXT_THEME_PRIORS.get(tag, {}).items():
                vector[theme_id] = round(max(vector.get(theme_id, 0.0), float(boost) * float(tag_weight)), 6)
        vector = {theme_id: weight for theme_id, weight in vector.items() if weight > 0.0}
        return {
            "tags": [{"tag": tag, "weight": weight} for tag, weight in sorted(tags.items())],
            "vector": vector,
        }

    def context_match(self, row_terms: list[str], context_terms: list[str], context_mode: str = "soft") -> dict[str, Any]:
        mode = context_mode if context_mode in {"soft", "hard"} else "soft"
        requested = set(context_terms)
        available = set(row_terms)
        if not requested:
            return {"passed": True, "multiplier": 1.0, "matched_terms": [], "mode": mode}
        if not available:
            return {
                "passed": mode != "hard",
                "multiplier": 0.85 if mode != "hard" else 0.0,
                "matched_terms": [],
                "mode": mode,
                "note": "no_context_terms_on_overlay_row",
            }
        matched = sorted(requested & available)
        if not matched:
            return {
                "passed": mode != "hard",
                "multiplier": 0.75 if mode != "hard" else 0.0,
                "matched_terms": [],
                "mode": mode,
                "note": "context_terms_not_matched",
            }
        return {
            "passed": True,
            "multiplier": round(min(1.35, 1.0 + 0.15 * len(matched)), 6),
            "matched_terms": matched,
            "mode": mode,
        }

    def ranking_term_index(self, rows: list[dict[str, Any]], uid_key: str) -> tuple[dict[str, dict[str, Any]], float]:
        index: dict[str, dict[str, Any]] = {}
        max_score = max((float(row.get("score") or 0.0) for row in rows), default=0.0)
        for row in rows:
            uid = str(row.get(uid_key) or "")
            node = self.get_node(uid) if uid else None
            values = [
                uid,
                row.get("display_name", ""),
                row.get("name", ""),
                row.get("primary_external_id", ""),
            ]
            if node:
                values.extend(
                    [
                        node.get("canonical_name", ""),
                        node.get("display_name", ""),
                        node.get("primary_external_id", ""),
                        node.get("external_xrefs", []),
                    ]
                )
            for term in prediction_terms_from_value(values):
                index.setdefault(term, row)
        return index, max(max_score, 1e-12)

    def ranking_exact_term_index(self, rows: list[dict[str, Any]], uid_key: str) -> tuple[dict[str, dict[str, Any]], float]:
        index: dict[str, dict[str, Any]] = {}
        max_score = max((float(row.get("score") or 0.0) for row in rows), default=0.0)
        for row in rows:
            uid = str(row.get(uid_key) or "")
            node = self.get_node(uid) if uid else None
            values = [
                uid,
                row.get("display_name", ""),
                row.get("name", ""),
                row.get("primary_external_id", ""),
            ]
            if node:
                values.extend(
                    [
                        node.get("node_uid", ""),
                        node.get("canonical_name", ""),
                        node.get("display_name", ""),
                        node.get("primary_external_id", ""),
                        node.get("external_xrefs", []),
                    ]
                )
            for term in exact_prediction_terms_from_value(values):
                index.setdefault(term, row)
        return index, max(max_score, 1e-12)

    def evidence_ref_for_overlay_row(self, row: dict[str, Any], target_uid: str = "") -> dict[str, Any]:
        return {
            "source_type": "prediction_overlay",
            "overlay_type": row.get("_overlay_type", ""),
            "source_file": row.get("_source_file", ""),
            "source_row": row.get("_source_row", 0),
            "source_name": row_value(row, "source_name", "source"),
            "source_record_id": row_value(row, "source_record_id", "record_id", "id"),
            "evidence_level": row_value(row, "evidence_level", "evidence"),
            "license_id": row_value(row, "license_id", "license"),
            "target_uid": target_uid,
        }

    def confidence_label(self, score: float) -> str:
        if score >= 0.75:
            return "high"
        if score >= 0.4:
            return "medium"
        return "low"

    def build_drug_predictions(
        self,
        targets: list[dict[str, Any]],
        context: dict[str, Any],
        context_mode: str,
    ) -> list[dict[str, Any]]:
        overlay = self.prediction_overlay_rows("drug_targets")
        if not overlay["rows"] or not targets:
            return []
        target_index, max_target_score = self.ranking_exact_term_index(targets, "target_uid")
        drug_rows: dict[str, dict[str, Any]] = {}
        tier_rank = {"low": 0, "exploratory": 1, "medium": 2, "high": 3}
        for row in overlay["rows"]:
            target_terms = split_overlay_exact_terms(row, "target_uid", "target_id", "target_symbol", "target_name", "gene_symbol", "gene")
            matched_target = next((target_index[term] for term in target_terms if term in target_index), None)
            if not matched_target:
                continue
            row_context_terms = split_overlay_terms(row, "context_terms", "cancer_type", "tissue", "cell_type", "cell_state", "cell_line")
            context_eval = self.context_match(row_context_terms, context.get("terms", []), context_mode)
            if not context_eval["passed"]:
                continue
            target_uid = str(matched_target.get("target_uid") or "")
            target_score_norm = clamp_unit(float(matched_target.get("score") or 0.0) / max_target_score)
            edge_confidence = clamp_unit(row_value(row, "confidence", "p_evidence", "activity_confidence"), 0.5)
            trace_multiplier = self.prediction_target_trace_multiplier(matched_target)
            component_score = clamp_unit(
                target_score_norm * edge_confidence * float(context_eval["multiplier"]) * trace_multiplier
            )
            drug_id = str(row_value(row, "drug_id", "chembl_id", "drug_uid", "id") or row_value(row, "drug_name", "name"))
            drug_name = str(row_value(row, "drug_name", "name", "pref_name") or drug_id)
            if not drug_id or not drug_name:
                continue
            normalized_drug = normalize_lookup_key(f"{drug_id} {drug_name} {row_value(row, 'mechanism', 'action_type', 'activity_type')}")
            downgrade_terms = [term for term in DRUG_OVERLAY_DOWNGRADE_TERMS if term in normalized_drug]
            candidate = drug_rows.setdefault(
                drug_id,
                {
                    "drug_id": drug_id,
                    "drug_name": drug_name,
                    "prediction_type": "drug_target_prioritization",
                    "_remaining": 1.0,
                    "matched_targets": [],
                    "mechanisms": set(),
                    "approval_statuses": set(),
                    "upstream_confidence_tiers": set(),
                    "downgrade_reasons": set(),
                    "evidence_refs": [],
                    "score_components": {
                        "formula": "1-prod(1-target_score_norm*edge_confidence*context_multiplier*trace_multiplier)",
                        "target_components": [],
                    },
                },
            )
            candidate["_remaining"] *= 1.0 - component_score
            mechanism = str(row_value(row, "mechanism", "action_type", "activity_type") or "unknown")
            approval_status = str(row_value(row, "approval_status", "max_phase", "development_phase") or "")
            candidate["mechanisms"].add(mechanism)
            if approval_status:
                candidate["approval_statuses"].add(approval_status)
            target_tier = str(matched_target.get("confidence_tier") or "low")
            candidate["upstream_confidence_tiers"].add(target_tier)
            if downgrade_terms:
                candidate["downgrade_reasons"].add("non_therapeutic_or_generic_overlay_entity")
            candidate["matched_targets"].append(
                {
                    "target_uid": target_uid,
                    "display_name": matched_target.get("display_name", ""),
                    "primary_external_id": matched_target.get("primary_external_id", ""),
                    "target_score": matched_target.get("score", 0.0),
                    "target_confidence_tier": target_tier,
                    "target_calibrated_confidence": matched_target.get("calibrated_confidence", 0.0),
                    "component_score": round(component_score, 6),
                    "context_match": context_eval,
                    "trace_multiplier": trace_multiplier,
                }
            )
            candidate["evidence_refs"].append(self.evidence_ref_for_overlay_row(row, target_uid=target_uid))
            candidate["score_components"]["target_components"].append(
                {
                    "target_uid": target_uid,
                    "target_score_norm": round(target_score_norm, 6),
                    "edge_confidence": round(edge_confidence, 6),
                    "context_multiplier": context_eval["multiplier"],
                    "trace_multiplier": trace_multiplier,
                    "component_score": round(component_score, 6),
                }
            )

        predictions = []
        for candidate in drug_rows.values():
            raw_overlay_score = clamp_unit(1.0 - float(candidate.pop("_remaining")))
            display_score = raw_overlay_score
            downgrade_reasons = candidate.pop("downgrade_reasons")
            if downgrade_reasons:
                display_score *= 0.35
            calibrated = self.bounded_confidence(display_score * PREDICTION_TYPE_WEIGHTS["drug"], precision=6)
            upstream_tiers = sorted(candidate.pop("upstream_confidence_tiers"), key=lambda tier: tier_rank.get(tier, 0), reverse=True)
            upstream_cap = upstream_tiers[0] if upstream_tiers else "low"
            initial_tier = self.confidence_tier(calibrated)
            final_tier = initial_tier
            if tier_rank.get(final_tier, 0) > tier_rank.get(upstream_cap, 0):
                final_tier = upstream_cap
                downgrade_reasons.add("drug_confidence_inherits_target_confidence")
            candidate["raw_overlay_score"] = round(raw_overlay_score, 6)
            candidate["score"] = round(display_score, 6)
            candidate["confidence"] = final_tier
            candidate["display_confidence"] = final_tier
            candidate["prediction"] = f"The model predicts {candidate.get('drug_name') or candidate.get('drug_id')} as a research-prioritization drug overlay lead."
            candidate["prediction_task"] = "drug_hypothesis_prediction"
            candidate["result_type"] = "drug"
            candidate["confidence_tier"] = final_tier
            candidate["calibrated_confidence"] = calibrated
            candidate["calibration_status"] = "overlay_informed"
            candidate["score_components"]["raw_overlay_score"] = round(raw_overlay_score, 6)
            candidate["score_components"]["display_score"] = round(display_score, 6)
            candidate["score_components"]["display_confidence"] = final_tier
            candidate["input_support_count"] = len({row.get("target_uid", "") for row in candidate.get("matched_targets", []) if row.get("target_uid")})
            candidate["graph_distance"] = 3
            candidate["is_directly_supported"] = False
            candidate["is_extrapolated"] = True
            candidate["needs_validation"] = True
            candidate["boundary"] = "Overlay-derived drug lead; not efficacy, causality, or clinical advice."
            candidate["upstream_confidence_tier"] = upstream_cap
            candidate["downgrade_reason"] = ";".join(sorted(downgrade_reasons))
            candidate["appendix"] = True
            candidate["appendix_reason"] = candidate["downgrade_reason"] or "overlay_informed_drug_hypothesis"
            candidate["mechanisms"] = sorted(candidate["mechanisms"])
            candidate["approval_statuses"] = sorted(candidate["approval_statuses"])
            candidate["matched_targets"].sort(key=lambda row: (-float(row.get("component_score") or 0.0), row.get("target_uid", "")))
            candidate["evidence_refs"] = candidate["evidence_refs"][:10]
            candidate["claim_refs"] = {
                "traceability_passed": bool(candidate["evidence_refs"]),
                "edge_uids": [],
                "evidence_ref_uids": [
                    str(ref.get("source_record_id") or ref.get("source_file") or "")
                    for ref in candidate["evidence_refs"]
                    if ref.get("source_record_id") or ref.get("source_file")
                ],
            }
            predictions.append(candidate)
        predictions.sort(key=lambda row: (-float(row.get("score") or 0.0), row.get("drug_name", ""), row.get("drug_id", "")))
        return predictions[: self.config.max_prediction_rows]

    def prediction_target_trace_multiplier(self, target: dict[str, Any]) -> float:
        components = target.get("score_components", {})
        propagation_score = float(components.get("propagation_score") or 0.0)
        path_confidence = float(components.get("best_path_confidence") or target.get("path_confidence") or 0.0)
        literature_count = int(components.get("literature_support_count") or target.get("literature_support_count") or 0)
        if path_confidence > 0.0 or propagation_score >= 0.00025:
            return 1.0
        if literature_count > 0:
            return 0.55
        return 0.35

    def build_cell_type_predictions(
        self,
        analysis_features: list[dict[str, Any]],
        pathways: list[dict[str, Any]],
        targets: list[dict[str, Any]],
        context: dict[str, Any],
        context_mode: str,
    ) -> list[dict[str, Any]]:
        overlay = self.prediction_overlay_rows("cell_type_signatures")
        if not overlay["rows"]:
            return []
        target_index, max_target_score = self.ranking_term_index(targets, "target_uid")
        pathway_index, max_pathway_score = self.ranking_term_index(pathways, "pathway_uid")
        metabolite_terms = {
            term
            for feature in analysis_features
            for term in prediction_terms_from_value(
                [feature.get("display_name", ""), feature.get("input_name", ""), feature.get("primary_external_id", "")]
            )
        }
        predictions = []
        for row in overlay["rows"]:
            row_context_terms = split_overlay_terms(row, "context_terms", "cancer_type", "tissue", "cell_type", "cell_state", "cell_line")
            context_eval = self.context_match(row_context_terms, context.get("terms", []), context_mode)
            if not context_eval["passed"]:
                continue
            target_matches = []
            target_support = 0.0
            for term in split_overlay_terms(row, "target_uids", "target_symbols", "gene_symbols", "marker_genes", "markers"):
                matched = target_index.get(term)
                if not matched:
                    continue
                support = clamp_unit(float(matched.get("score") or 0.0) / max_target_score)
                target_support = max(target_support, support)
                target_matches.append({"target_uid": matched.get("target_uid", ""), "display_name": matched.get("display_name", ""), "support": round(support, 6)})
            pathway_matches = []
            pathway_support = 0.0
            for term in split_overlay_terms(row, "pathway_uids", "pathway_terms", "pathways"):
                matched = pathway_index.get(term)
                if not matched:
                    for pathway in pathways:
                        name = normalize_prediction_text(pathway.get("name") or pathway.get("display_name") or "")
                        if term and term in name:
                            matched = pathway
                            break
                if not matched:
                    continue
                support = clamp_unit(float(matched.get("score") or 0.0) / max_pathway_score)
                pathway_support = max(pathway_support, support)
                pathway_matches.append({"pathway_uid": matched.get("pathway_uid", ""), "display_name": matched.get("name") or matched.get("display_name", ""), "support": round(support, 6)})
            requested_metabolite_terms = split_overlay_terms(row, "metabolite_terms", "metabolites")
            metabolite_matches = sorted(set(requested_metabolite_terms) & metabolite_terms)
            metabolite_support = 1.0 if metabolite_matches else 0.0
            if not target_matches and not pathway_matches and not metabolite_matches and not context_eval["matched_terms"]:
                continue
            row_confidence = clamp_unit(row_value(row, "confidence", "signature_confidence"), 0.5)
            score = clamp_unit(
                (0.6 * target_support + 0.3 * pathway_support + 0.1 * metabolite_support)
                * row_confidence
                * float(context_eval["multiplier"])
            )
            if score <= 0.0 and context_eval["matched_terms"]:
                score = clamp_unit(0.15 * row_confidence * float(context_eval["multiplier"]))
            cell_type_id = str(row_value(row, "cell_type_id", "cell_id", "id") or row_value(row, "cell_type", "cell_type_name", "name"))
            cell_type_name = str(row_value(row, "cell_type_name", "cell_type", "name") or cell_type_id)
            if not cell_type_id or not cell_type_name:
                continue
            calibrated = self.bounded_confidence(score * PREDICTION_TYPE_WEIGHTS["cell_context_overlay"], precision=6)
            predictions.append(
                {
                    "cell_type_id": cell_type_id,
                    "cell_type_name": cell_type_name,
                    "prediction_type": "cell_context_compatibility",
                    "prediction": f"The model predicts {cell_type_name} as a research-prioritization context overlay lead.",
                    "prediction_task": "cell_context_prediction",
                    "result_type": "cell_context_overlay",
                    "confidence_tier": self.confidence_tier(calibrated),
                    "calibrated_confidence": calibrated,
                    "calibration_status": "overlay_informed",
                    "input_support_count": len(metabolite_matches),
                    "graph_distance": 3,
                    "is_directly_supported": False,
                    "is_extrapolated": True,
                    "needs_validation": True,
                    "boundary": "Overlay-derived context lead; validate independently before interpretation.",
                    "score": round(score, 6),
                    "confidence": self.confidence_label(score),
                    "context_match": context_eval,
                    "matched_targets": target_matches[:10],
                    "matched_pathways": pathway_matches[:10],
                    "matched_metabolites": metabolite_matches[:10],
                    "score_components": {
                        "formula": "(0.6*target_support+0.3*pathway_support+0.1*metabolite_support)*signature_confidence*context_multiplier",
                        "target_support": round(target_support, 6),
                        "pathway_support": round(pathway_support, 6),
                        "metabolite_support": round(metabolite_support, 6),
                        "signature_confidence": round(row_confidence, 6),
                        "context_multiplier": context_eval["multiplier"],
                    },
                    "evidence_refs": [self.evidence_ref_for_overlay_row(row)],
                    "claim_refs": {
                        "traceability_passed": True,
                        "edge_uids": [],
                        "evidence_ref_uids": [str(row_value(row, "source_record_id", "cell_type_id", "cell_type_name") or "")],
                    },
                }
            )
        predictions.sort(key=lambda row: (-float(row.get("score") or 0.0), row.get("cell_type_name", ""), row.get("cell_type_id", "")))
        return predictions[: self.config.max_prediction_rows]

    def build_prediction_pack(
        self,
        records: list[Any],
        analysis_features: list[dict[str, Any]],
        pathways: list[dict[str, Any]],
        targets: list[dict[str, Any]],
        diseases: list[dict[str, Any]],
        context: Any = None,
        context_mode: str = "soft",
    ) -> dict[str, Any]:
        mode = context_mode if context_mode in {"soft", "hard"} else "soft"
        normalized_context = (
            context
            if isinstance(context, dict) and "terms" in context and "fields" in context
            else self.normalize_prediction_context(context, records)
        )
        drug_overlay = self.prediction_overlay_rows("drug_targets")
        cell_overlay = self.prediction_overlay_rows("cell_type_signatures")
        drug_rankings = self.build_drug_predictions(targets, normalized_context, mode)
        cell_type_rankings = self.build_cell_type_predictions(analysis_features, pathways, targets, normalized_context, mode)
        blocked_reasons = []
        if not drug_overlay["rows"]:
            blocked_reasons.append({"code": "drug_overlay_missing", "message": "No drug_targets overlay rows were found for this release."})
        if not cell_overlay["rows"]:
            blocked_reasons.append(
                {"code": "cell_type_overlay_missing", "message": "No cell_type_signatures overlay rows were found for this release."}
            )
        if drug_overlay["rows"] and not drug_rankings:
            blocked_reasons.append({"code": "no_drug_target_overlap", "message": "Drug overlay rows did not overlap ranked targets or context filters."})
        if cell_overlay["rows"] and not cell_type_rankings:
            blocked_reasons.append(
                {"code": "no_cell_signature_overlap", "message": "Cell-type overlay rows did not overlap ranked targets/pathways/metabolites or context filters."}
            )
        pack = {
            "contract_version": PREDICTION_PACK_CONTRACT_VERSION,
            "prediction_scope": "research_prioritization_not_clinical_decision_support",
            "context_mode": mode,
            "context": normalized_context,
            "overlay_status": {
                "drug_targets": {
                    "row_count": drug_overlay["row_count"],
                    "files": drug_overlay["files"],
                },
                "cell_type_signatures": {
                    "row_count": cell_overlay["row_count"],
                    "files": cell_overlay["files"],
                },
            },
            "drug_rankings": drug_rankings,
            "cell_type_rankings": cell_type_rankings,
            "blocked_reasons": blocked_reasons,
            "limitations": [
                "Predictions are generated from ranked analysis targets/pathways plus local overlay rows; they are not canonical graph facts.",
                "Use hard context mode only when overlay rows have explicit tissue, cancer, cell-type, or cell-line annotations.",
                "Candidate drugs require independent review of mechanism, assay context, dose, and safety before any experimental use.",
            ],
            "determinism": {
                "release_id": self.release_id,
                "context_hash": content_hash(normalized_context),
                "overlay_hash": content_hash(
                    {
                        "drug_targets": drug_overlay["rows"],
                        "cell_type_signatures": cell_overlay["rows"],
                    }
                ),
            },
        }
        pack["determinism"]["prediction_pack_hash"] = content_hash({key: value for key, value in pack.items() if key != "determinism"})
        return pack

    def _direct_uid_candidate(self, query: str, entity_type: str | None) -> list[dict[str, Any]]:
        raw = str(query or "").strip()
        if not raw:
            return []
        node = self.get_node(raw)
        if not node:
            return []
        if entity_type and node.get("node_type") != entity_type:
            return []
        return [
            {
                "entity_uid": node["node_uid"],
                "entity_type": node["node_type"],
                "namespace": "UID",
                "lookup_key": node["node_uid"],
                "raw_value": node["node_uid"],
                "match_field": "node_uid",
                "rank": 100.0,
                "source_table": "nodes",
                "source_release": node.get("source_release", ""),
                "license_id": node.get("license_id", ""),
                "_query_priority": 120,
            }
        ]

    def _resolver_rows(
        self,
        query: str,
        entity_type: str | None = None,
        namespace: str | None = None,
    ) -> list[dict[str, Any]]:
        rows = self._direct_uid_candidate(query, entity_type)
        resolver_path = self.graph_dir / "resolver_index.parquet"
        if namespace:
            key = normalize_lookup_key(query) if namespace.upper() == "TEXT" else normalize_id_token(query)
            candidates = [(namespace.upper(), key, 100)]
        else:
            candidates = candidate_queries(query)

        if resolver_path.exists():
            dataset = table_dataset(resolver_path)
            for ns, key, priority in candidates:
                filt = (ds.field("namespace") == ns) & (ds.field("lookup_key") == key)
                if entity_type:
                    filt = filt & (ds.field("entity_type") == entity_type)
                for row in dataset.to_table(filter=filt).to_pylist():
                    row["_query_priority"] = priority
                    rows.append(row)

        pubchem_lookup = self.pubchem_dir / "cid_lookup_index.parquet"
        if pubchem_lookup.exists() and (entity_type in {None, "metabolite"}):
            dataset = table_dataset(pubchem_lookup)
            for ns, key, priority in candidates:
                filt = (ds.field("namespace") == ns) & (ds.field("lookup_key") == key)
                for row in dataset.to_table(filter=filt).to_pylist():
                    for metabolite_uid in parse_list(row.get("metabolite_uids")):
                        rows.append(
                            {
                                "entity_uid": metabolite_uid,
                                "entity_type": "metabolite",
                                "namespace": row.get("namespace", ""),
                                "lookup_key": row.get("lookup_key", ""),
                                "raw_value": row.get("raw_value", ""),
                                "match_field": f"pubchem_cache:{row.get('match_field', '')}",
                                "rank": row.get("rank", 0.0),
                                "source_table": "cid_lookup_index",
                                "source_release": row.get("source_release", ""),
                                "license_id": row.get("license_id", ""),
                                "_query_priority": priority,
                                "pubchem_cid": row.get("pubchem_cid", ""),
                            }
                        )
        return rows

    def resolve(
        self,
        query: str,
        entity_type: str | None = None,
        namespace: str | None = None,
        limit: int | None = None,
    ) -> dict[str, Any]:
        request = {
            "query": query,
            "entity_type": entity_type,
            "namespace": namespace,
            "limit": limit or self.config.max_candidates,
        }
        rows = self._resolver_rows(query, entity_type=entity_type, namespace=namespace)
        by_entity: dict[str, dict[str, Any]] = {}
        for row in rows:
            entity_uid = row.get("entity_uid", "")
            if not entity_uid:
                continue
            score, components = score_resolver_row(row)
            candidate = by_entity.get(entity_uid)
            match = {
                "namespace": row.get("namespace", ""),
                "lookup_key": row.get("lookup_key", ""),
                "raw_value": row.get("raw_value", ""),
                "match_field": row.get("match_field", ""),
                "rank": float(row.get("rank") or 0.0),
                "source_table": row.get("source_table", ""),
                "license_id": row.get("license_id", ""),
                "source_release": row.get("source_release", ""),
            }
            if row.get("pubchem_cid"):
                match["pubchem_cid"] = row.get("pubchem_cid")
            if candidate is None:
                by_entity[entity_uid] = {
                    "entity_uid": entity_uid,
                    "entity_type": row.get("entity_type", ""),
                    "score": score,
                    "score_components": components,
                    "best_match": match,
                    "matches": [match],
                }
            else:
                candidate["matches"].append(match)
                if score > candidate["score"]:
                    candidate["score"] = score
                    candidate["score_components"] = components
                    candidate["best_match"] = match
        candidates = sorted(
            by_entity.values(),
            key=lambda item: (-item["score"], item["entity_type"], item["entity_uid"]),
        )
        for candidate in candidates:
            node = self.get_node(candidate["entity_uid"])
            if node:
                candidate["display_name"] = node.get("display_name", "")
                candidate["primary_external_id"] = node.get("primary_external_id", "")
                candidate["external_xrefs"] = node.get("external_xrefs", [])
            candidate["matches"] = sorted(
                candidate["matches"],
                key=lambda item: (
                    -float(item.get("rank") or 0.0),
                    item.get("namespace", ""),
                    item.get("match_field", ""),
                    item.get("raw_value", ""),
                ),
            )[:5]
        candidates = candidates[: limit or self.config.max_candidates]
        top1 = candidates[0]["score"] if candidates else 0.0
        top2 = candidates[1]["score"] if len(candidates) > 1 else 0.0
        margin = round(top1 - top2, 6)
        if not candidates:
            status = "unmatched"
        elif top1 >= self.config.min_match_score and margin >= self.config.min_match_margin:
            status = "matched"
        else:
            status = "ambiguous"
        payload = {
            "query": query,
            "entity_type_filter": entity_type,
            "status": status,
            "top_score": top1,
            "top_margin": margin,
            "thresholds": {
                "top1_min": self.config.min_match_score,
                "top1_minus_top2_min": self.config.min_match_margin,
            },
            "candidates": candidates,
        }
        return self.response_envelope("/resolve", request, payload)

    def compound_index_available(self) -> bool:
        return all(
            (self.compound_dir / f"{name}.parquet").exists()
            for name in (
                "compound_identifier_index",
                "compound_name_index",
                "compound_inchikey_index",
                "compound_formula_mass_index",
            )
        )

    def normalized_record_keys(self, record: Any) -> dict[str, Any]:
        if not isinstance(record, dict):
            return {}
        return {str(key).strip().casefold().replace(" ", "_"): value for key, value in record.items()}

    def european_trait_index(self) -> dict[str, dict[str, Any]]:
        if self._european_trait_index is not None:
            return self._european_trait_index
        path = european_source_path(self.workspace, "European.csv")
        index: dict[str, dict[str, Any]] = {}
        if not path.exists():
            self._european_trait_index = index
            return index
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.reader(handle)
            try:
                header = next(reader)
            except StopIteration:
                self._european_trait_index = index
                return index
            canonical_header = [canonical_chat_column(cell) for cell in header]
            raw_header = [str(cell or "").strip() for cell in header]
            for row in reader:
                fields: dict[str, Any] = {}
                for col_index, cell in enumerate(row):
                    if col_index >= len(canonical_header):
                        continue
                    canonical = canonical_header[col_index]
                    raw = raw_header[col_index]
                    value = str(cell or "").strip()
                    if not value:
                        continue
                    if canonical and canonical not in fields:
                        fields[canonical] = value
                    if raw and raw not in fields:
                        fields[raw] = value
                accession = fields.get("accession_id") or fields.get("accessionId") or fields.get("trait")
                accession_key = normalize_gwas_accession(accession)
                if accession_key:
                    fields.setdefault("accession_id", accession_key)
                    index[accession_key] = fields
        self._european_trait_index = index
        return index

    def european_trait_annotation_index(self) -> dict[str, dict[str, Any]]:
        if self._european_trait_annotation_index is not None:
            return self._european_trait_annotation_index
        path = european_source_path(self.workspace, "European_trait_annotations.csv")
        index: dict[str, dict[str, Any]] = {}
        if not path.exists():
            self._european_trait_annotation_index = index
            return index
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle):
                accession = normalize_gwas_accession(row.get("accession_id"))
                if accession:
                    index[accession] = dict(row)
        self._european_trait_annotation_index = index
        return index

    def split_annotation_values(self, value: Any) -> list[str]:
        text = str(value or "").strip()
        if not text:
            return []
        delimiter = "|" if "|" in text else ";"
        return [part.strip() for part in text.split(delimiter) if part.strip()]

    def enrich_compound_record_with_european_trait(self, record: Any) -> Any:
        if not isinstance(record, dict):
            return record
        normalized = self.normalized_record_keys(record)
        if any(normalized.get(key) not in {None, ""} for key in ("name", "metabolite", "compound", "hmdb", "hmdb_id", "chebi", "chebi_id", "pubchem_cid", "cid", "kegg", "kegg_id", "inchikey", "inchi_key")):
            return record
        trait_value = normalized.get("trait") or normalized.get("accession_id") or normalized.get("accessionid") or normalized.get("gcst")
        reported_trait = normalized.get("reported_trait") or normalized.get("reportedtrait")
        trait_row = None
        accession_key = normalize_gwas_accession(trait_value)
        annotation_row = self.european_trait_annotation_index().get(accession_key) if accession_key else None
        annotation_names = self.split_annotation_values(annotation_row.get("mapped_names") if annotation_row else "")
        local_annotation_cids = self.split_annotation_values(annotation_row.get("local_pubchem_cid_hits") if annotation_row else "")
        strict_manual_scope = any(
            scope.strip().casefold() == "strict_identity"
            for scope in self.split_annotation_values(annotation_row.get("manual_identity_scope") if annotation_row else "")
        )
        if annotation_row and not reported_trait:
            reported_trait = annotation_row.get("reported_trait")
        if accession_key:
            trait_row = self.european_trait_index().get(accession_key)
            if trait_row and not reported_trait:
                reported_trait = trait_row.get("reported_trait") or trait_row.get("reportedTrait")
        candidate_names = dedupe_preserve_order([*annotation_names, *reported_trait_candidate_names(reported_trait)])
        if not candidate_names:
            return record
        enriched = dict(record)
        enriched.setdefault("reported_trait", str(reported_trait or ""))
        if accession_key:
            enriched.setdefault("accession_id", accession_key)
        ratio_trait = ratio_trait_descriptor(
            reported_trait,
            str(enriched.get("effect_label") or enriched.get("comparison_effect_source") or ""),
            enriched.get("effect_value", enriched.get("log2FC")),
            str(enriched.get("direction") or ""),
        )
        if ratio_trait:
            enriched["ratio_trait"] = ratio_trait
        if local_annotation_cids:
            enriched.setdefault("pubchem_cid", local_annotation_cids[0])
        if strict_manual_scope and annotation_row:
            manual_cids = self.split_annotation_values(annotation_row.get("reviewed_pubchem_cid"))
            manual_hmdb_ids = self.split_annotation_values(annotation_row.get("reviewed_hmdb_id"))
            manual_chebi_ids = self.split_annotation_values(annotation_row.get("reviewed_chebi_id"))
            manual_inchikeys = self.split_annotation_values(annotation_row.get("reviewed_inchikey"))
            if manual_cids:
                enriched.setdefault("pubchem_cid", manual_cids[0])
            if manual_hmdb_ids:
                enriched.setdefault("hmdb", manual_hmdb_ids[0])
            if manual_chebi_ids:
                enriched.setdefault("chebi", manual_chebi_ids[0])
            if manual_inchikeys:
                enriched.setdefault("inchikey", manual_inchikeys[0])
        if trait_row:
            for target_key, source_keys in {
                "summary_statistics_url": ("summary_statistics_url", "summaryStatistics"),
                "pubmed_id": ("pubmed_id", "pubmedId"),
                "paper_title": ("paper_title", "title"),
                "journal": ("journal",),
                "publication_date": ("publication_date", "publicationDate"),
                "efo_traits": ("efo_traits", "efoTraits"),
                "bg_traits": ("bg_traits", "bgTraits"),
                "initial_sample_description": ("initial_sample_description", "initialSampleDescription"),
                "discovery_sample_ancestry": ("discovery_sample_ancestry", "discoverySampleAncestry"),
            }.items():
                for source_key in source_keys:
                    value = trait_row.get(source_key)
                    if value not in {None, ""}:
                        enriched.setdefault(target_key, value)
                        break
        enriched["european_trait_candidate_names"] = candidate_names
        enriched["name"] = candidate_names[0]
        return enriched

    def add_compound_match(
        self,
        candidates: dict[str, dict[str, Any]],
        metabolite_uid: str,
        component: str,
        component_score: float,
        match: dict[str, Any],
    ) -> None:
        if not metabolite_uid:
            return
        candidate = candidates.setdefault(
            metabolite_uid,
            {
                "entity_uid": metabolite_uid,
                "entity_type": "metabolite",
                "score": 0.0,
                "score_components": {
                    "identifier": 0.0,
                    "structure": 0.0,
                    "mass": 0.0,
                    "formula": 0.0,
                    "name": 0.0,
                    "rt": 0.0,
                    "ms2": 0.0,
                    "context": 0.0,
                    "graph": 0.0,
                    "ambiguity_penalty": 0.0,
                },
                "matches": [],
            },
        )
        candidate["score_components"][component] = max(
            float(candidate["score_components"].get(component, 0.0)),
            float(component_score),
        )
        candidate["matches"].append(match)

    def expanded_compound_name_queries(self, names: list[str]) -> list[tuple[str, str]]:
        queries: list[tuple[str, str]] = []
        seen: set[str] = set()
        for name in names:
            raw = str(name or "").strip()
            if not raw:
                continue
            aliases = COMMON_METABOLITE_NAME_ALIASES.get(normalize_lookup_key(raw), [])
            spelling_variants = [
                raw,
                re.sub(r"(?<=[A-Za-z])-(?=\d)", " ", raw),
                raw.replace("-", " "),
            ]
            for value, mode in (
                [(variant, "input_name_variant" if variant != raw else "input_name") for variant in spelling_variants]
                + [(alias, f"common_alias:{raw}") for alias in aliases]
            ):
                key = normalize_lookup_key(value)
                if not key or key in seen:
                    continue
                seen.add(key)
                queries.append((value, mode))
        return queries

    def compound_record_has_orthogonal_support(self, parsed: dict[str, Any]) -> bool:
        return bool(
            parsed.get("identifiers")
            or parsed.get("inchikey")
            or parsed.get("formula")
            or parsed.get("mz") not in {None, ""}
            or parsed.get("rt") not in {None, ""}
            or parsed.get("ms2_peaks")
        )

    def compound_name_triage_policy(
        self,
        name_key: str,
        rows: list[dict[str, Any]],
        parsed: dict[str, Any],
    ) -> dict[str, Any]:
        uids = sorted({str(row.get("metabolite_uid") or "") for row in rows if row.get("metabolite_uid")})
        reasons = []
        has_orthogonal_support = self.compound_record_has_orthogonal_support(parsed)
        lipid_like = bool(
            re.search(r"\b[A-Z]{1,4}\(?\d{1,2}:\d", name_key, flags=re.IGNORECASE)
            and re.search(r"[/_]\d{1,2}:\d", name_key)
        )
        delimiter_sensitive = any(token in name_key for token in (",", ";"))
        max_rank = max((float(row.get("rank") or 0.0) for row in rows), default=0.0)
        if len(uids) > 1:
            reasons.append("name_maps_to_multiple_metabolites")
        if lipid_like:
            reasons.append("lipid_shorthand_or_chain_notation")
        if delimiter_sensitive:
            reasons.append("delimiter_sensitive_name")
        if rows and max_rank < 50.0:
            reasons.append("synonym_only_low_rank_surface")

        policy = "monitor"
        skip_name_component = False
        score_cap = 65.0
        if "name_maps_to_multiple_metabolites" in reasons and not has_orthogonal_support:
            policy = "auto_abstain_name_only"
            skip_name_component = True
            score_cap = 0.0
        elif "lipid_shorthand_or_chain_notation" in reasons and not has_orthogonal_support:
            policy = "requires_identifier_or_orthogonal_feature"
            skip_name_component = True
            score_cap = 0.0
        elif "name_maps_to_multiple_metabolites" in reasons or "lipid_shorthand_or_chain_notation" in reasons:
            policy = "monitor_or_downweight"
            score_cap = 25.0
        elif "delimiter_sensitive_name" in reasons:
            policy = "format_warning"
        elif "synonym_only_low_rank_surface" in reasons:
            policy = "monitor_or_downweight"
            score_cap = 35.0

        return {
            "policy": policy,
            "reasons": reasons,
            "skip_name_component": skip_name_component,
            "score_cap": score_cap,
            "unique_metabolite_uid_count": len(uids),
            "has_orthogonal_support": has_orthogonal_support,
            "max_rank": round(max_rank, 6),
        }

    def metabolite_graph_support_score(self, metabolite_uid: str) -> tuple[float, dict[str, Any]]:
        if not metabolite_uid:
            return 0.0, {"graph_edge_count": 0, "literature_support_count": 0, "max_p_literature": 0.0}
        if metabolite_uid in self._metabolite_graph_support_cache:
            score, meta = self._metabolite_graph_support_cache[metabolite_uid]
            return score, dict(meta)
        edge_count = 0
        edge_path = self.graph_dir / "edges.parquet"
        if edge_path.exists():
            edge_count = len(
                table_rows(
                    edge_path,
                    ["edge_uid"],
                    (ds.field("subject_uid") == metabolite_uid) | (ds.field("object_uid") == metabolite_uid),
                )
            )
        self.load_literature_support()
        support_rows = [
            row
            for row in (self._literature_support_by_entity or {}).get(metabolite_uid, [])
            if row.get("support_class") in {"confirm", "support_direction", "novel_candidate"}
            and float(row.get("p_literature") or 0.0) >= self.config.min_literature_overlay_prob
        ]
        max_p_literature = max((float(row.get("p_literature") or 0.0) for row in support_rows), default=0.0)
        score = 0.0
        if edge_count:
            score += 5.0 + min(13.0, math.log2(edge_count + 1.0) * 3.0)
        node = self.get_node(metabolite_uid) or {}
        xrefs = parse_list(node.get("external_xrefs"))
        if any(xref.startswith("HMDB:") for xref in xrefs):
            score += 2.0
        if any(xref.startswith("KEGG.COMPOUND:") for xref in xrefs):
            score += 1.0
        if support_rows:
            score += min(4.0, math.log2(len(support_rows) + 1.0))
            score += min(1.5, max_p_literature * 1.5)
        result = round(min(18.0, score), 6), {
            "graph_edge_count": edge_count,
            "literature_support_count": len(support_rows),
            "max_p_literature": round(max_p_literature, 6),
        }
        self._metabolite_graph_support_cache[metabolite_uid] = result
        return result[0], dict(result[1])

    def query_compound_identifier(self, namespace: str, value: str) -> list[dict[str, Any]]:
        path = self.compound_dir / "compound_identifier_index.parquet"
        if not path.exists():
            return []
        namespace = namespace.upper()
        lookup_key = normalize_id_token(value)
        if not lookup_key:
            return []
        cache_key = (namespace, lookup_key)
        if cache_key not in self._compound_identifier_query_cache:
            self._compound_identifier_query_cache[cache_key] = table_rows(
                path,
                None,
                (ds.field("namespace") == namespace) & (ds.field("lookup_key") == lookup_key),
            )
        return [dict(row) for row in self._compound_identifier_query_cache[cache_key]]

    def query_compound_name(self, name_key: str) -> list[dict[str, Any]]:
        path = self.compound_dir / "compound_name_index.parquet"
        if not path.exists() or not name_key:
            return []
        if name_key not in self._compound_name_query_cache:
            self._compound_name_query_cache[name_key] = table_rows(path, None, ds.field("name_key") == name_key)
        return [dict(row) for row in self._compound_name_query_cache[name_key]]

    def query_compound_formula(self, formula_key: str) -> list[dict[str, Any]]:
        path = self.compound_dir / "compound_formula_mass_index.parquet"
        if not path.exists() or not formula_key:
            return []
        if formula_key not in self._compound_formula_query_cache:
            self._compound_formula_query_cache[formula_key] = table_rows(path, None, ds.field("formula_key") == formula_key)
        return [dict(row) for row in self._compound_formula_query_cache[formula_key]]

    def query_compound_mass_range(self, neutral_mass: float, ppm_tolerance: float) -> list[dict[str, Any]]:
        path = self.compound_dir / "compound_formula_mass_index.parquet"
        if not path.exists() or not neutral_mass:
            return []
        cache_key = (round(float(neutral_mass), 8), round(float(ppm_tolerance), 6))
        if cache_key not in self._compound_mass_range_cache:
            delta = neutral_mass * ppm_tolerance / 1_000_000.0
            filt = (ds.field("exact_mass") >= neutral_mass - delta) & (ds.field("exact_mass") <= neutral_mass + delta)
            self._compound_mass_range_cache[cache_key] = table_rows(path, None, filt)
        return [dict(row) for row in self._compound_mass_range_cache[cache_key]]

    def query_compound_inchikey(self, inchikey: str) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
        path = self.compound_dir / "compound_inchikey_index.parquet"
        if not path.exists() or not inchikey:
            return [], []
        key = normalize_id_token(inchikey).upper()
        if key not in self._compound_inchikey_full_cache:
            self._compound_inchikey_full_cache[key] = table_rows(path, None, ds.field("inchikey") == key)
        conn_key = inchikey_connectivity(key)
        if conn_key not in self._compound_inchikey_connectivity_cache:
            self._compound_inchikey_connectivity_cache[conn_key] = table_rows(path, None, ds.field("connectivity_key") == conn_key)
        return (
            [dict(row) for row in self._compound_inchikey_full_cache[key]],
            [dict(row) for row in self._compound_inchikey_connectivity_cache[conn_key]],
        )

    def query_compound_rt_candidates(self, candidate_uids: list[str]) -> list[dict[str, Any]]:
        path = self.compound_dir / "compound_rt_index.parquet"
        cache_key = tuple(sorted(candidate_uids))
        if not path.exists() or not cache_key:
            return []
        if cache_key not in self._compound_rt_candidate_cache:
            self._compound_rt_candidate_cache[cache_key] = table_rows(path, None, ds.field("metabolite_uid").isin(list(cache_key)))
        return [dict(row) for row in self._compound_rt_candidate_cache[cache_key]]

    def query_compound_ms2_candidates(self, candidate_uids: list[str]) -> list[dict[str, Any]]:
        path = self.compound_dir / "compound_ms2_index.parquet"
        cache_key = tuple(sorted(candidate_uids))
        if not path.exists() or not cache_key:
            return []
        if cache_key not in self._compound_ms2_candidate_cache:
            self._compound_ms2_candidate_cache[cache_key] = table_rows(path, None, ds.field("metabolite_uid").isin(list(cache_key)))
        return [dict(row) for row in self._compound_ms2_candidate_cache[cache_key]]

    def compound_record_queries(self, record: Any) -> dict[str, Any]:
        if isinstance(record, str):
            return {"input_id": "", "record": {"name": record}, "names": [record], "identifiers": []}
        record = self.enrich_compound_record_with_european_trait(record)
        normalized = self.normalized_record_keys(record)
        identifiers = []
        for key, namespace in [
            ("hmdb_id", "HMDB"),
            ("hmdb", "HMDB"),
            ("chebi_id", "CHEBI"),
            ("chebi", "CHEBI"),
            ("pubchem_cid", "CID"),
            ("cid", "CID"),
            ("pubchem", "CID"),
            ("kegg_id", "KEGG"),
            ("kegg", "KEGG"),
        ]:
            value = normalized.get(key)
            if value not in {None, ""}:
                identifiers.append({"field": key, "namespace": namespace, "value": str(value).strip()})
        names = [
            str(normalized[key]).strip()
            for key in ("name", "metabolite", "compound")
            if normalized.get(key) not in {None, ""}
        ]
        for name in normalized.get("european_trait_candidate_names") or []:
            text = str(name or "").strip()
            if text and text not in names:
                names.append(text)
        return {
            "input_id": str(normalized.get("input_id") or normalized.get("id") or ""),
            "record": record if isinstance(record, dict) else {"name": record},
            "identifiers": identifiers,
            "names": names,
            "inchikey": str(normalized.get("inchikey") or normalized.get("inchi_key") or "").strip(),
            "formula": str(normalized.get("formula") or "").strip(),
            "mz": normalized.get("mz", normalized.get("m/z")),
            "adduct": normalized.get("adduct", ""),
            "charge": normalized.get("charge"),
            "ppm_tolerance": normalized.get("ppm_tolerance", normalized.get("ppm")),
            "rt": normalized.get("rt", normalized.get("retention_time")),
            "rt_tolerance": normalized.get("rt_tolerance"),
            "rt_method": normalized.get("rt_method", normalized.get("method", normalized.get("lc_method", "default"))),
            "ms2_peaks": parse_ms2_peaks(
                normalized.get("ms2_peaks")
                or normalized.get("msms_peaks")
                or normalized.get("fragments")
                or normalized.get("fragment_peaks")
            ),
            "ms2_mz_tolerance": normalized.get("ms2_mz_tolerance", normalized.get("fragment_mz_tolerance")),
            "ion_mode": str(normalized.get("ion_mode") or normalized.get("polarity") or "").strip(),
        }

    def finalize_compound_candidates(self, candidates: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
        rows = list(candidates.values())
        name_candidate_count = sum(1 for row in rows if row["score_components"].get("name", 0.0) > 0)
        mass_candidate_count = sum(1 for row in rows if row["score_components"].get("mass", 0.0) > 0)
        for row in rows:
            components = row["score_components"]
            if components.get("name", 0.0) > 0:
                graph_score, graph_meta = self.metabolite_graph_support_score(row["entity_uid"])
                components["graph"] = graph_score
                row["_graph_support_meta"] = graph_meta
            strong = max(components["identifier"], components["structure"], components["name"], components["mass"])
            support = 0.0
            if components["formula"] and strong < 95.0:
                support += components["formula"]
            if components["rt"] and strong >= 40.0:
                support += components["rt"]
            if components["ms2"] and strong >= 40.0:
                support += components["ms2"]
            if components["context"] and strong >= 40.0:
                support += components["context"]
            if components.get("graph") and components["name"] > 0 and strong >= 45.0:
                support += components["graph"]
            if components["identifier"] and components["structure"]:
                support += 5.0
            if components["name"] and strong >= 70.0:
                support += min(5.0, components["name"] / 13.0)
            if components.get("common_biochemical_zero_score_rescue"):
                strong = max(strong, self.config.min_match_score + min(10.0, float(components["common_biochemical_zero_score_rescue"]) / 4.0))
            penalty = 0.0
            if components["identifier"] == 0 and components["structure"] < 90:
                if components["name"] > 0 and name_candidate_count > 1:
                    penalty += min(12.0, math.log2(name_candidate_count) * 4.0)
                if components["mass"] > 0 and mass_candidate_count > 10:
                    penalty += min(18.0, math.log10(mass_candidate_count) * 8.0)
            components["ambiguity_penalty"] = round(-penalty, 6)
            row["score"] = round(max(0.0, min(100.0, strong + support - penalty)), 6)
            node = self.get_node(row["entity_uid"])
            if node:
                row["display_name"] = node.get("display_name", "")
                row["primary_external_id"] = node.get("primary_external_id", "")
                row["external_xrefs"] = node.get("external_xrefs", [])
            else:
                name_match = next((item for item in row["matches"] if item.get("component") == "name"), None)
                identifier_match = next((item for item in row["matches"] if item.get("component") == "identifier"), None)
                structure_match = next((item for item in row["matches"] if item.get("component") == "structure"), None)
                if name_match:
                    row["display_name"] = str(name_match.get("raw_value") or "")
                elif identifier_match:
                    namespace = str(identifier_match.get("namespace") or "").strip()
                    raw_value = str(identifier_match.get("raw_value") or "").strip()
                    row["display_name"] = f"{namespace}:{raw_value}" if namespace and raw_value else row["entity_uid"]
                else:
                    row["display_name"] = row["entity_uid"]
                external_xrefs = []
                if identifier_match:
                    namespace = str(identifier_match.get("namespace") or "").strip()
                    raw_value = str(identifier_match.get("raw_value") or "").strip()
                    if namespace and raw_value:
                        external_xrefs.append(f"{namespace}:{raw_value}")
                if structure_match and structure_match.get("raw_value"):
                    external_xrefs.append(f"INCHIKEY:{structure_match.get('raw_value')}")
                row["external_xrefs"] = list(dict.fromkeys(external_xrefs))
                row["primary_external_id"] = row["external_xrefs"][0] if row["external_xrefs"] else ""
            row["matches"] = sorted(
                row["matches"],
                key=lambda item: (
                    item.get("component", ""),
                    -float(item.get("component_score", 0.0)),
                    item.get("source_table", ""),
                    item.get("raw_value", ""),
                ),
            )[:10]
        return sorted(rows, key=lambda item: (-item["score"], item["entity_uid"]))

    def classify_candidate_list(self, candidates_list: list[dict[str, Any]]) -> tuple[str, float, float]:
        top1 = candidates_list[0]["score"] if candidates_list else 0.0
        top2 = candidates_list[1]["score"] if len(candidates_list) > 1 else 0.0
        margin = round(top1 - top2, 6)
        if not candidates_list:
            return "unmatched", top1, margin
        top_components = candidates_list[0].get("score_components", {})
        if top1 >= self.config.min_match_score and margin >= self.config.min_match_margin:
            return "matched", top1, margin
        second_components = candidates_list[1].get("score_components", {}) if len(candidates_list) > 1 else {}
        if (
            top1 >= 95.0
            and float(top_components.get("identifier") or 0.0) >= 88.0
            and float(top_components.get("structure") or 0.0) >= 95.0
            and float(second_components.get("identifier") or 0.0) < 80.0
        ):
            return "matched", top1, margin
        if (
            top1 >= self.config.min_match_score + 10.0
            and float(top_components.get("name") or 0.0) > 0.0
            and float(top_components.get("graph") or 0.0) >= 8.0
            and (
                margin >= 0.5
                or (
                    top1 >= self.config.min_match_score + 15.0
                    and float(top_components.get("graph") or 0.0) >= 16.0
                )
            )
        ):
            return "matched", top1, margin
        return "ambiguous", top1, margin

    def maybe_rescue_common_biochemical_match(
        self,
        status: str,
        candidates_list: list[dict[str, Any]],
        parsed: dict[str, Any],
        query_notes: list[dict[str, Any]],
    ) -> str:
        if status == "matched" or not candidates_list:
            return status
        names = [str(name or "").strip() for name in parsed.get("names", []) if str(name or "").strip()]
        if not names:
            return status
        input_name = names[0]
        if not is_common_biochemical_rescue_query(input_name):
            return status
        zero_score_rescue = [
            (common_biochemical_zero_score_rescue_rank(input_name, candidate), index, candidate)
            for index, candidate in enumerate(candidates_list)
            if float(candidate.get("score") or 0.0) <= 0.0
        ]
        zero_score_rescue = [item for item in zero_score_rescue if item[0] >= 10.0]
        if zero_score_rescue:
            zero_score_rescue.sort(key=lambda item: (-item[0], item[1]))
            selected_rank, _selected_index, selected = zero_score_rescue[0]
            selected["score"] = max(float(selected.get("score") or 0.0), self.config.min_match_score + min(10.0, selected_rank / 4.0))
            selected.setdefault("score_components", {})["common_biochemical_zero_score_rescue"] = round(selected_rank, 6)
            selected.setdefault("matches", []).append(
                {
                    "component": "common_biochemical_zero_score_rescue",
                    "component_score": round(selected_rank, 6),
                    "field": "name",
                    "raw_value": input_name,
                    "match_field": "common_biochemical_exact_name",
                    "reason": "common biochemical exact-name rescue after conservative name-only abstention",
                    "source_table": "service_common_biochemical_rescue",
                }
            )
            candidates_list.sort(key=lambda item: (-float(item.get("score") or 0.0), item.get("entity_uid", "")))
            query_notes.append(
                {
                    "mode": "common_biochemical_zero_score_rescue",
                    "active": True,
                    "input_name": input_name,
                    "selected_entity_uid": selected.get("entity_uid", ""),
                    "selected_display_name": selected.get("display_name", ""),
                    "previous_status": status,
                    "rescue_rank": round(selected_rank, 6),
                }
            )
            return "matched"
        top = candidates_list[0]
        top_score = float(top.get("score") or 0.0)
        if top_score < 45.0:
            return status
        if not biochemical_rescue_compatible(input_name, top.get("display_name", "")):
            return status
        top.setdefault("score_components", {})["biological_exact_rescue"] = round(max(0.0, self.config.min_match_score - top_score), 6)
        top.setdefault("matches", []).append(
            {
                "component": "biological_exact_rescue",
                "component_score": 1.0,
                "field": "name",
                "raw_value": input_name,
                "match_field": "common_biochemical_name",
                "reason": "common biochemical exact-name rescue for salt, charge, stereochemistry, or naming variant",
                "source_table": "service_common_biochemical_rescue",
            }
        )
        query_notes.append(
            {
                "mode": "biological_exact_rescue",
                "active": True,
                "input_name": input_name,
                "selected_entity_uid": top.get("entity_uid", ""),
                "selected_display_name": top.get("display_name", ""),
                "previous_status": status,
            }
        )
        return "matched"

    def apply_pathway_context_rerank(self, row_contracts: list[dict[str, Any]]) -> None:
        seed_uids = [
            row["resolution"]["candidates"][0]["entity_uid"]
            for row in row_contracts
            if row["resolution"].get("resolver_version") == "compound_resolver_v2"
            and row["resolution"]["status"] == "matched"
            and row["resolution"].get("candidates")
        ]
        if len(seed_uids) < 1 or len(row_contracts) < 2:
            return
        index = self.build_pathway_index()
        seed_counts: dict[str, int] = defaultdict(int)
        for seed_uid in seed_uids:
            for pathway_uid in index["metabolite_to_pathways"].get(seed_uid, set()):
                seed_counts[pathway_uid] += 1
        if not seed_counts:
            return

        for row in row_contracts:
            if row["resolution"].get("resolver_version") != "compound_resolver_v2":
                continue
            candidates_list = row["resolution"].get("candidates") or []
            if not candidates_list:
                continue
            candidate_map = {candidate["entity_uid"]: candidate for candidate in candidates_list}
            changed = False
            for candidate in candidate_map.values():
                candidate_pathways = index["metabolite_to_pathways"].get(candidate["entity_uid"], set())
                shared_weight = sum(seed_counts.get(pathway_uid, 0) for pathway_uid in candidate_pathways)
                if shared_weight <= 0:
                    continue
                support = min(8.0, 4.0 * shared_weight / max(1, len(seed_uids)))
                components = candidate.setdefault("score_components", {})
                if support <= float(components.get("context", 0.0)):
                    continue
                components["context"] = round(support, 6)
                shared_pathways = sorted(pathway_uid for pathway_uid in candidate_pathways if seed_counts.get(pathway_uid, 0))
                candidate.setdefault("matches", []).append(
                    {
                        "component": "context",
                        "component_score": round(support, 6),
                        "field": "batch_pathway_context",
                        "raw_value": f"{len(shared_pathways)} shared pathway(s)",
                        "match_field": "pathway_membership",
                        "shared_pathway_uids": shared_pathways[:10],
                        "seed_metabolite_count": len(seed_uids),
                        "source_table": "metabolite_pathway_edges/reaction_participants",
                    }
                )
                changed = True
            if not changed:
                continue
            refreshed = self.finalize_compound_candidates(candidate_map)[: self.config.max_candidates]
            status, top_score, top_margin = self.classify_candidate_list(refreshed)
            row["resolution"]["query_notes"] = normalize_query_notes(row["resolution"].get("query_notes", []))
            if any(
                note.get("mode") in {"biological_exact_rescue", "common_biochemical_zero_score_rescue"} and note.get("active")
                for note in row["resolution"].get("query_notes", [])
            ):
                status = "matched"
            row["resolution"]["candidates"] = refreshed
            row["resolution"]["status"] = status
            row["resolution"]["top_score"] = top_score
            row["resolution"]["top_margin"] = top_margin
            row["resolution"].setdefault("query_notes", []).append(
                {
                    "mode": "pathway_context",
                    "active": True,
                    "seed_metabolite_count": len(seed_uids),
                    "seed_pathway_count": len(seed_counts),
                }
            )

    def add_rt_support(self, candidates: dict[str, dict[str, Any]], parsed: dict[str, Any], query_notes: list[dict[str, Any]]) -> None:
        rt_value = maybe_number(parsed.get("rt"))
        if rt_value is None:
            return
        path = self.compound_dir / "compound_rt_index.parquet"
        if not path.exists():
            query_notes.append({"mode": "rt", "active": False, "reason": "compound_rt_index_missing"})
            return
        candidate_uids = sorted(candidates)
        if not candidate_uids:
            query_notes.append({"mode": "rt", "active": True, "hits": 0, "reason": "no_prior_candidates"})
            return
        tolerance = maybe_number(parsed.get("rt_tolerance")) or self.config.default_rt_tolerance
        method_key = normalize_method_key(parsed.get("rt_method"))
        rows = self.query_compound_rt_candidates(candidate_uids)
        hits = 0
        for row in rows:
            row_method = str(row.get("method_key") or "default")
            if method_key != "default" and row_method not in {method_key, "default"}:
                continue
            candidate_rt = maybe_number(row.get("rt"))
            if candidate_rt is None:
                continue
            delta = abs(candidate_rt - rt_value)
            if delta > tolerance:
                continue
            method_factor = 1.0 if row_method == method_key else 0.85
            component_score = max(0.0, 12.0 * (1.0 - delta / tolerance) * method_factor)
            hits += 1
            self.add_compound_match(
                candidates,
                row.get("metabolite_uid", ""),
                "rt",
                component_score,
                {
                    "component": "rt",
                    "component_score": round(component_score, 6),
                    "field": "rt",
                    "raw_value": row.get("rt", ""),
                    "match_field": "retention_time",
                    "rt_delta": round(delta, 6),
                    "rt_tolerance": tolerance,
                    "method_key": row_method,
                    "source_table": "compound_rt_index",
                    "source_name": row.get("source_name", ""),
                    "source_record_id": row.get("source_record_id", ""),
                    "source_release": row.get("source_release", ""),
                    "license_id": row.get("license_id", ""),
                },
            )
        query_notes.append(
            {
                "mode": "rt",
                "active": True,
                "rt": rt_value,
                "rt_tolerance": tolerance,
                "method_key": method_key,
                "candidate_rows": len(rows),
                "hits": hits,
            }
        )

    def add_ms2_support(self, candidates: dict[str, dict[str, Any]], parsed: dict[str, Any], query_notes: list[dict[str, Any]]) -> None:
        query_peaks = parsed.get("ms2_peaks") or []
        if not query_peaks:
            return
        path = self.compound_dir / "compound_ms2_index.parquet"
        if not path.exists():
            query_notes.append({"mode": "ms2", "active": False, "reason": "compound_ms2_index_missing"})
            return
        candidate_uids = sorted(candidates)
        if not candidate_uids:
            query_notes.append({"mode": "ms2", "active": True, "hits": 0, "reason": "no_prior_candidates"})
            return
        tolerance = maybe_number(parsed.get("ms2_mz_tolerance")) or self.config.default_ms2_mz_tolerance
        input_adduct = normalize_adduct(parsed.get("adduct"))
        input_ion_mode = str(parsed.get("ion_mode") or "").strip().casefold()
        rows = self.query_compound_ms2_candidates(candidate_uids)
        by_uid: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for row in rows:
            row_adduct = normalize_adduct(row.get("adduct"))
            row_ion_mode = str(row.get("ion_mode") or "").strip().casefold()
            if input_adduct and row_adduct and row_adduct != input_adduct:
                continue
            if input_ion_mode and row_ion_mode and row_ion_mode != input_ion_mode:
                continue
            by_uid[row.get("metabolite_uid", "")].append(row)

        total_intensity = sum(max(0.0, float(intensity or 0.0)) for _mz, intensity in query_peaks) or float(len(query_peaks))
        hits = 0
        for metabolite_uid, lib_rows in by_uid.items():
            used: set[int] = set()
            matched = []
            matched_intensity = 0.0
            for query_mz, query_intensity in sorted(query_peaks, key=lambda item: (-item[1], item[0])):
                best_idx = None
                best_delta = tolerance + 1.0
                best_row = None
                for idx, row in enumerate(lib_rows):
                    if idx in used:
                        continue
                    lib_mz = maybe_number(row.get("fragment_mz"))
                    if lib_mz is None:
                        continue
                    delta = abs(lib_mz - query_mz)
                    if delta <= tolerance and delta < best_delta:
                        best_idx = idx
                        best_delta = delta
                        best_row = row
                if best_idx is None or best_row is None:
                    continue
                used.add(best_idx)
                matched_intensity += max(0.0, float(query_intensity or 0.0))
                matched.append(
                    {
                        "query_mz": round(query_mz, 6),
                        "library_mz": round(float(best_row.get("fragment_mz") or 0.0), 6),
                        "delta": round(best_delta, 6),
                        "source_record_id": best_row.get("source_record_id", ""),
                    }
                )
            if not matched:
                continue
            intensity_fraction = matched_intensity / total_intensity
            count_fraction = len(matched) / max(1, min(len(query_peaks), 6))
            component_score = min(20.0, 20.0 * (0.7 * intensity_fraction + 0.3 * min(1.0, count_fraction)))
            hits += 1
            self.add_compound_match(
                candidates,
                metabolite_uid,
                "ms2",
                component_score,
                {
                    "component": "ms2",
                    "component_score": round(component_score, 6),
                    "field": "ms2_peaks",
                    "raw_value": f"{len(query_peaks)} query peaks",
                    "match_field": "fragment_mz",
                    "matched_fragments": matched[:10],
                    "matched_peak_count": len(matched),
                    "query_peak_count": len(query_peaks),
                    "fragment_mz_tolerance": tolerance,
                    "source_table": "compound_ms2_index",
                    "source_record_ids": sorted({item["source_record_id"] for item in matched if item.get("source_record_id")})[:5],
                },
            )
        query_notes.append(
            {
                "mode": "ms2",
                "active": True,
                "query_peak_count": len(query_peaks),
                "fragment_mz_tolerance": tolerance,
                "candidate_rows": len(rows),
                "hits": hits,
            }
        )

    def compound_resolve_record(self, record: Any, row_index: int = 0) -> dict[str, Any]:
        parsed = self.compound_record_queries(record)
        input_id = parsed.get("input_id") or f"row_{row_index}"
        if not self.compound_index_available():
            fallback = self.detect_record_query(record, row_index)
            if fallback.get("invalid_reason"):
                return {
                    "input_id": input_id,
                    "record": record,
                    "status": "invalid",
                    "invalid_reason": fallback["invalid_reason"],
                    "top_score": 0.0,
                    "top_margin": 0.0,
                    "candidates": [],
                }
            result = self.resolve(fallback["query"], entity_type="metabolite")
            return {
                "input_id": input_id,
                "query": fallback["query"],
                "query_field": fallback["query_field"],
                "record": fallback["record"],
                "resolver_version": "resolver_index_fallback",
                "status": result["status"],
                "top_score": result["top_score"],
                "top_margin": result["top_margin"],
                "candidates": result["candidates"],
            }

        candidates: dict[str, dict[str, Any]] = {}
        query_notes = []
        for identifier in parsed["identifiers"]:
            namespace = identifier["namespace"]
            value = identifier["value"]
            clean_value = value
            if namespace == "CID":
                clean_value = normalize_id_token(clean_value).removeprefix("pubchem:").removeprefix("cid:")
            elif ":" in clean_value:
                clean_value = clean_value.split(":", 1)[1]
            rows = self.query_compound_identifier(namespace, clean_value)
            query_notes.append({"mode": "identifier", "field": identifier["field"], "namespace": namespace, "hits": len(rows)})
            for row in rows:
                component_score = 92.0 if namespace != "CID" else 88.0
                self.add_compound_match(
                    candidates,
                    row.get("metabolite_uid", ""),
                    "identifier",
                    component_score,
                    {
                        "component": "identifier",
                        "component_score": component_score,
                        "field": identifier["field"],
                        "namespace": namespace,
                        "raw_value": row.get("raw_value", ""),
                        "match_field": row.get("match_field", ""),
                        "rank": row.get("rank", 0.0),
                        "source_table": row.get("source_table", ""),
                        "source_release": row.get("source_release", ""),
                        "license_id": row.get("license_id", ""),
                    },
                )

        inchikey = parsed.get("inchikey", "")
        if inchikey:
            key = normalize_id_token(inchikey).upper()
            full_rows, conn_rows = self.query_compound_inchikey(key)
            query_notes.append({"mode": "inchikey", "hits": len(full_rows), "connectivity_hits": len(conn_rows)})
            full_uids = {row.get("metabolite_uid") for row in full_rows}
            for row in full_rows:
                self.add_compound_match(
                    candidates,
                    row.get("metabolite_uid", ""),
                    "structure",
                    95.0,
                    {
                        "component": "structure",
                        "component_score": 95.0,
                        "field": "inchikey",
                        "raw_value": row.get("raw_value", ""),
                        "match_field": "inchikey_full",
                        "source_table": row.get("source_table", ""),
                        "source_release": row.get("source_release", ""),
                        "license_id": row.get("license_id", ""),
                    },
                )
            for row in conn_rows:
                if row.get("metabolite_uid") in full_uids:
                    continue
                self.add_compound_match(
                    candidates,
                    row.get("metabolite_uid", ""),
                    "structure",
                    70.0,
                    {
                        "component": "structure",
                        "component_score": 70.0,
                        "field": "inchikey",
                        "raw_value": row.get("raw_value", ""),
                        "match_field": "inchikey_connectivity",
                        "source_table": row.get("source_table", ""),
                        "source_release": row.get("source_release", ""),
                        "license_id": row.get("license_id", ""),
                    },
                )

        formula_key = normalize_formula(parsed.get("formula", ""))
        mass_path = self.compound_dir / "compound_formula_mass_index.parquet"
        if formula_key and mass_path.exists():
            rows = self.query_compound_formula(formula_key)
            query_notes.append({"mode": "formula", "formula_key": formula_key, "hits": len(rows)})
            for row in rows:
                self.add_compound_match(
                    candidates,
                    row.get("metabolite_uid", ""),
                    "formula",
                    20.0,
                    {
                        "component": "formula",
                        "component_score": 20.0,
                        "field": "formula",
                        "raw_value": row.get("formula", ""),
                        "match_field": row.get("match_field", ""),
                        "source_table": row.get("source_table", ""),
                        "source_release": row.get("source_release", ""),
                        "license_id": row.get("license_id", ""),
                    },
                )

        neutral_mass, mass_meta = neutral_mass_from_mz(parsed.get("mz"), parsed.get("adduct"), parsed.get("charge"))
        if neutral_mass and mass_path.exists():
            ppm_tolerance = maybe_number(parsed.get("ppm_tolerance")) or self.config.default_ppm_tolerance
            rows = self.query_compound_mass_range(neutral_mass, ppm_tolerance)
            query_notes.append(
                {
                    "mode": "mass",
                    "neutral_mass": round(neutral_mass, 8),
                    "ppm_tolerance": ppm_tolerance,
                    "hits": len(rows),
                    **mass_meta,
                }
            )
            for row in rows:
                ppm_error = abs(float(row.get("exact_mass", 0.0)) - neutral_mass) / neutral_mass * 1_000_000.0
                component_score = max(0.0, 50.0 * (1.0 - min(ppm_error, ppm_tolerance) / ppm_tolerance))
                self.add_compound_match(
                    candidates,
                    row.get("metabolite_uid", ""),
                    "mass",
                    component_score,
                    {
                        "component": "mass",
                        "component_score": round(component_score, 6),
                        "field": "mz",
                        "raw_value": row.get("exact_mass", ""),
                        "match_field": row.get("match_field", ""),
                        "mass_error_ppm": round(ppm_error, 6),
                        "neutral_mass": round(neutral_mass, 8),
                        "ppm_tolerance": ppm_tolerance,
                        "source_table": row.get("source_table", ""),
                        "source_release": row.get("source_release", ""),
                        "license_id": row.get("license_id", ""),
                    },
                )

        for name, name_mode in self.expanded_compound_name_queries(parsed.get("names", [])):
            key = normalize_lookup_key(name)
            rows = self.query_compound_name(key)
            triage = self.compound_name_triage_policy(key, rows, parsed)
            query_notes.append(
                {
                    "mode": "name",
                    "name_key": key,
                    "name_mode": name_mode,
                    "hits": len(rows),
                    "triage_policy": triage["policy"],
                    "triage_reasons": triage["reasons"],
                    "unique_metabolite_uid_count": triage["unique_metabolite_uid_count"],
                    "has_orthogonal_support": triage["has_orthogonal_support"],
                }
            )
            if triage["skip_name_component"]:
                for row in rows:
                    self.add_compound_match(
                        candidates,
                        row.get("metabolite_uid", ""),
                        "name",
                        0.0,
                        {
                            "component": "name",
                            "component_score": 0.0,
                            "field": "name",
                            "raw_value": row.get("raw_value", ""),
                            "match_field": row.get("match_field", ""),
                            "name_mode": name_mode,
                            "rank": row.get("rank", 0.0),
                            "source_table": row.get("source_table", ""),
                            "source_release": row.get("source_release", ""),
                            "license_id": row.get("license_id", ""),
                            "triage_policy": triage["policy"],
                            "triage_blocked_scoring": True,
                        },
                    )
                continue
            for row in rows:
                component_score = min(float(triage["score_cap"]), float(row.get("rank") or 0.0) * 0.65)
                if name_mode.startswith("common_alias:") and row.get("match_field") in {"canonical_name", "pubchem_title", "synonym"}:
                    component_score = max(component_score, min(65.0, float(triage["score_cap"])))
                self.add_compound_match(
                    candidates,
                    row.get("metabolite_uid", ""),
                    "name",
                    component_score,
                    {
                        "component": "name",
                        "component_score": round(component_score, 6),
                        "field": "name",
                        "raw_value": row.get("raw_value", ""),
                        "match_field": row.get("match_field", ""),
                        "name_mode": name_mode,
                        "rank": row.get("rank", 0.0),
                        "source_table": row.get("source_table", ""),
                        "source_release": row.get("source_release", ""),
                        "license_id": row.get("license_id", ""),
                    },
                )

        self.add_rt_support(candidates, parsed, query_notes)
        self.add_ms2_support(candidates, parsed, query_notes)

        candidates_list = self.finalize_compound_candidates(candidates)[: self.config.max_candidates]
        status, top1, margin = self.classify_candidate_list(candidates_list)
        status = self.maybe_rescue_common_biochemical_match(status, candidates_list, parsed, query_notes)
        if status == "matched" and candidates_list:
            top1 = float(candidates_list[0].get("score") or top1)
            top2 = float(candidates_list[1].get("score") or 0.0) if len(candidates_list) > 1 else 0.0
            margin = round(top1 - top2, 6)
        has_supported_input = bool(
            parsed["identifiers"]
            or parsed.get("names")
            or parsed.get("inchikey")
            or parsed.get("formula")
            or parsed.get("mz") not in {None, ""}
            or parsed.get("rt") not in {None, ""}
            or parsed.get("ms2_peaks")
        )
        if not has_supported_input:
            status = "invalid"
        return {
            "input_id": input_id,
            "record": parsed["record"],
            "resolver_version": "compound_resolver_v2",
            "status": status,
            "top_score": top1,
            "top_margin": margin,
            "query_notes": query_notes,
            "thresholds": {
                "top1_min": self.config.min_match_score,
                "top1_minus_top2_min": self.config.min_match_margin,
            },
            "candidates": candidates_list,
            "invalid_reason": "" if has_supported_input else "no_supported_compound_identifier_or_feature",
        }

    def european_annotation_for_record(self, record: Any) -> dict[str, Any]:
        if not isinstance(record, dict):
            return {}
        normalized = self.normalized_record_keys(record)
        accession = normalize_gwas_accession(
            first_present(normalized, ("accession_id", "accessionid", "trait", "gcst", "gwas_trait"))
        )
        return self.european_trait_annotation_index().get(accession, {}) if accession else {}

    def european_trait_source_evidence_for_record(self, record: Any) -> dict[str, Any]:
        if not isinstance(record, dict):
            return {}
        normalized = self.normalized_record_keys(record)
        accession = normalize_gwas_accession(
            first_present(normalized, ("accession_id", "accessionid", "trait", "gcst", "gwas_trait"))
        )
        annotation = self.european_trait_annotation_index().get(accession, {}) if accession else {}
        trait_row = self.european_trait_index().get(accession, {}) if accession else {}

        def pick(*keys: str) -> str:
            for source in (normalized, annotation, trait_row):
                for key in keys:
                    value = source.get(key)
                    if value is not None and str(value).strip():
                        return str(value).strip()
            return ""

        summary_url = pick("summary_statistics_url", "summarystatistics", "summaryStatistics")
        reported_trait = pick("reported_trait", "reportedtrait", "reportedTrait")
        if not accession and not summary_url and not reported_trait:
            return {}
        evidence = {
            "accession_id": accession,
            "reported_trait": reported_trait,
            "summary_statistics_url": summary_url,
            "pubmed_id": pick("pubmed_id", "pubmedid", "pubmedId"),
            "paper_title": pick("paper_title", "title"),
            "journal": pick("journal"),
            "publication_date": pick("publication_date", "publicationdate", "publicationDate"),
            "efo_traits": pick("efo_traits", "efotraits", "efoTraits"),
            "bg_traits": pick("bg_traits", "bgtraits", "bgTraits"),
            "initial_sample_description": pick(
                "initial_sample_description", "initialsampledescription", "initialSampleDescription"
            ),
            "discovery_sample_ancestry": pick(
                "discovery_sample_ancestry", "discoverysampleancestry", "discoverySampleAncestry"
            ),
            "annotation_status": str(annotation.get("resolution_status") or ""),
            "annotation_source_name": str(annotation.get("source_name") or ""),
            "annotation_pubchem_cids": str(annotation.get("pubchem_cids") or ""),
            "annotation_pubchem_titles": str(annotation.get("pubchem_titles") or ""),
            "manual_review_status": str(annotation.get("manual_review_status") or ""),
            "manual_identity_scope": str(annotation.get("manual_identity_scope") or ""),
            "manual_evidence_source": str(annotation.get("evidence_source") or ""),
            "manual_evidence_url": str(annotation.get("evidence_url") or ""),
        }
        if summary_url:
            evidence["summary_statistics_url_type"] = "gwas_summary_statistics"
        return {key: value for key, value in evidence.items() if value not in {None, ""}}

    def trait_identity_review_reasons(self, row_contract: dict[str, Any]) -> list[str]:
        record = row_contract.get("record", {})
        annotation = self.european_annotation_for_record(record)
        normalized = self.normalized_record_keys(record)
        label = " ".join(
            str(value or "")
            for value in (
                first_present(normalized, ("reported_trait", "reportedtrait", "name", "metabolite", "compound")),
                annotation.get("reported_trait", ""),
                annotation.get("mapped_names", ""),
                annotation.get("candidate_names", ""),
            )
            if value
        )
        reasons: list[str] = []
        if trait_text_is_ratio_or_composite(label):
            reasons.append("ratio_or_composite_trait_not_strict_identity")
        if trait_text_is_pool_or_class(label):
            reasons.append("pool_or_class_trait_not_strict_identity")
        return reasons

    def apply_trait_identity_resolution_policy(self, row_contract: dict[str, Any]) -> dict[str, Any]:
        resolution = row_contract.get("resolution", {})
        reasons = self.trait_identity_review_reasons(row_contract)
        if reasons:
            resolution["identity_review_reasons"] = reasons
            notes = list(resolution.get("query_notes") or [])
            for reason in reasons:
                if reason not in notes:
                    notes.append(reason)
            resolution["query_notes"] = notes
            if resolution.get("status") == "matched":
                resolution["status"] = "ambiguous"
        return row_contract

    def expanded_candidate_class_for_row(self, row_contract: dict[str, Any]) -> tuple[str, str]:
        record = row_contract.get("record", {})
        annotation = self.european_annotation_for_record(record)
        ann_status = str(annotation.get("resolution_status") or "")
        candidates = row_contract.get("resolution", {}).get("candidates") or []
        if not candidates:
            return "unresolved", "no_candidate_entities"
        normalized = self.normalized_record_keys(record)
        label = " ".join(
            str(value or "")
            for value in (
                first_present(normalized, ("reported_trait", "name", "metabolite", "compound")),
                annotation.get("reported_trait", ""),
                annotation.get("mapped_names", ""),
                annotation.get("candidate_names", ""),
            )
            if value
        )
        if trait_text_is_ratio_or_composite(label):
            return "ratio_component", "european_ratio_or_composite_trait"
        if trait_text_is_pool_or_class(label):
            return "class_or_pool", "lipid_or_pool_like_trait"
        if ann_status == "pubchem_cid_in_local_index":
            return "strict_identity", "european_pubchem_cid_in_local_index"
        if ann_status == "local_name":
            return "soft_identity", "european_local_name_candidate"
        if ann_status == "pubchem_only":
            return "analog_candidate", "european_pubchem_only_not_in_strict_local_index"
        top_score = float(row_contract.get("resolution", {}).get("top_score") or 0.0)
        if top_score >= 45.0:
            return "soft_identity", "resolver_score_above_soft_threshold"
        if top_score > 0.0:
            return "analog_candidate", "low_score_name_or_graph_candidate"
        return "unresolved", ann_status or "insufficient_candidate_evidence"

    def expanded_candidate_limit(self, seed_class: str) -> int:
        if seed_class == "strict_identity":
            return 1
        if seed_class == "soft_identity":
            return 2
        if seed_class in {"ratio_component", "class_or_pool"}:
            return 3
        if seed_class == "analog_candidate":
            return 2
        return 0

    def expanded_rows_for_resolution_row(self, row_contract: dict[str, Any]) -> list[dict[str, Any]]:
        seed_class, reason = self.expanded_candidate_class_for_row(row_contract)
        total_weight = float(EXPANDED_SEED_WEIGHTS.get(seed_class, 0.0))
        limit = self.expanded_candidate_limit(seed_class)

        def candidate_score(candidate: dict[str, Any]) -> float:
            try:
                return float(candidate.get("score") or 0.0)
            except (TypeError, ValueError):
                return 0.0

        candidates = [
            candidate
            for candidate in (row_contract.get("resolution", {}).get("candidates") or [])
            if candidate.get("entity_uid") and candidate_score(candidate) > 0.0
        ][:limit]
        identity_cluster_by_uid: dict[str, dict[str, Any]] = {}
        if seed_class == "soft_identity":
            clusters = identity_clusters_for_candidates(candidates)
            candidates = [cluster["representative"] for cluster in clusters if cluster.get("representative", {}).get("entity_uid")]
            for cluster in clusters:
                representative_uid = str(cluster.get("representative_uid") or "")
                if representative_uid:
                    identity_cluster_by_uid[representative_uid] = {key: value for key, value in cluster.items() if key != "representative"}
            candidates = candidates[:limit]
        if not candidates or total_weight <= 0.0:
            return []
        per_candidate_weight = total_weight / max(1, len(candidates))
        rows = []
        for index, candidate in enumerate(candidates):
            expanded = copy.deepcopy(row_contract)
            expanded["metabolite_uid"] = candidate.get("entity_uid", "")
            expanded["input_id"] = f"{row_contract.get('input_id', '')}::expanded_{index + 1}"
            expanded["resolution"] = {
                **dict(row_contract.get("resolution", {})),
                "status": "expanded_candidate",
                "candidates": [candidate],
                "top_score": candidate.get("score", row_contract.get("resolution", {}).get("top_score", 0.0)),
                "top_margin": row_contract.get("resolution", {}).get("top_margin", 0.0),
            }
            expanded["expanded_candidate"] = {
                "track": "expanded",
                "seed_class": seed_class,
                "weight_multiplier": round(per_candidate_weight, 6),
                "total_class_weight": total_weight,
                "candidate_count_for_input": len(candidates),
                "reason": reason,
                "source_status": row_contract.get("resolution", {}).get("status", ""),
                "original_input_id": row_contract.get("input_id", ""),
            }
            if seed_class == "ratio_component":
                ratio_component = ratio_component_for_candidate(row_contract.get("record", {}), candidate)
                if ratio_component:
                    expanded["expanded_candidate"]["ratio_trait"] = ratio_component.get("ratio_trait", {})
                    expanded["expanded_candidate"]["ratio_component"] = ratio_component.get("component", {})
            elif seed_class == "class_or_pool":
                class_seed = class_seed_descriptor(row_contract.get("record", {}), candidate)
                if class_seed:
                    expanded["expanded_candidate"]["class_seed"] = class_seed
            elif seed_class == "soft_identity":
                cluster = identity_cluster_by_uid.get(str(candidate.get("entity_uid") or ""))
                if cluster:
                    expanded["expanded_candidate"]["identity_cluster"] = cluster
                    expanded["expanded_candidate"]["reason"] = "identity_cluster_consensus_candidate"
            rows.append(expanded)
        return rows

    def expanded_candidate_rows(self, ambiguous_rows: list[dict[str, Any]], unmatched_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
        rows = []
        for row in [*ambiguous_rows, *unmatched_rows]:
            rows.extend(self.expanded_rows_for_resolution_row(row))
        return rows

    def expanded_candidate_summary(
        self,
        matched_rows: list[dict[str, Any]],
        expanded_rows: list[dict[str, Any]],
        ambiguous_rows: list[dict[str, Any]],
        unmatched_rows: list[dict[str, Any]],
        invalid_rows: list[dict[str, Any]],
    ) -> dict[str, Any]:
        by_class = Counter(row.get("expanded_candidate", {}).get("seed_class", "unknown") for row in expanded_rows)
        expanded_input_ids = {
            row.get("expanded_candidate", {}).get("original_input_id") or row.get("input_id", "")
            for row in expanded_rows
        }
        strict_input_ids = {row.get("input_id", "") for row in matched_rows}
        unresolved_input_ids = {
            row.get("input_id", "")
            for row in [*ambiguous_rows, *unmatched_rows, *invalid_rows]
            if row.get("input_id", "") not in expanded_input_ids
        }
        weighted_seed_mass = len(strict_input_ids) + sum(
            float(row.get("expanded_candidate", {}).get("weight_multiplier") or 0.0) for row in expanded_rows
        )
        return {
            "strict_matched_count": len(matched_rows),
            "strict_input_count": len(strict_input_ids),
            "expanded_candidate_count": len(expanded_rows),
            "expanded_candidate_input_count": len(expanded_input_ids),
            "expanded_input_count": len(expanded_input_ids),
            "expanded_by_class": dict(sorted(by_class.items())),
            "unresolved_input_count": len(unresolved_input_ids),
            "analysis_seed_input_count": len(strict_input_ids | expanded_input_ids),
            "analysis_seed_feature_count": len(matched_rows) + len(expanded_rows),
            "weighted_seed_mass": round(weighted_seed_mass, 6),
            "weight_policy": EXPANDED_SEED_WEIGHTS,
            "strict_semantics": "matched rows remain high-confidence near-unique chemical identities; expanded rows are low-weight research signals.",
        }

    def detect_record_query(self, record: Any, row_index: int) -> dict[str, Any]:
        if isinstance(record, str):
            return {"input_id": f"row_{row_index}", "query_field": "name", "query": record, "record": {"name": record}}
        if not isinstance(record, dict):
            return {"input_id": f"row_{row_index}", "query_field": "", "query": "", "record": record, "invalid_reason": "record_is_not_object_or_string"}
        record = self.enrich_compound_record_with_european_trait(record)
        normalized = {str(key).strip().casefold().replace(" ", "_"): value for key, value in record.items()}
        field_order = [
            ("metabolite_uid", "metabolite_uid", None),
            ("hmdb_id", "HMDB", "HMDB"),
            ("hmdb", "HMDB", "HMDB"),
            ("chebi_id", "CHEBI", "CHEBI"),
            ("chebi", "CHEBI", "CHEBI"),
            ("pubchem_cid", "CID", "CID"),
            ("cid", "CID", "CID"),
            ("pubchem", "CID", "CID"),
            ("kegg_id", "KEGG", "KEGG"),
            ("kegg", "KEGG", "KEGG"),
            ("inchikey", "INCHIKEY", "INCHIKEY"),
            ("inchi_key", "INCHIKEY", "INCHIKEY"),
            ("name", "TEXT", None),
            ("metabolite", "TEXT", None),
            ("compound", "TEXT", None),
        ]
        for key, label, namespace in field_order:
            value = normalized.get(key)
            if value is None or str(value).strip() == "":
                continue
            query = str(value).strip()
            if label in {"HMDB", "CHEBI", "KEGG", "CID"} and ":" not in query:
                query = f"{label}:{query}" if label != "CID" else f"PubChem:{query}"
            return {
                "input_id": str(record.get("input_id") or record.get("id") or f"row_{row_index}"),
                "query_field": key,
                "query": query,
                "namespace": namespace,
                "record": record,
            }
        return {
            "input_id": str(record.get("input_id") or record.get("id") or f"row_{row_index}"),
            "query_field": "",
            "query": "",
            "record": record,
            "invalid_reason": "no_supported_identifier_column",
        }

    def precheck_metabolites(self, records: list[Any], input_normalization: dict[str, Any] | None = None) -> dict[str, Any]:
        if input_normalization is None:
            records, input_normalization = prepare_analysis_records(records)
        request = {"records": records}
        self.prime_compound_resolver_caches(records)
        row_contracts = []
        for idx, record in enumerate(records):
            result = self.compound_resolve_record(record, idx)
            row_contract = {
                "input_id": result["input_id"],
                "query": result.get("query", ""),
                "query_field": result.get("query_field", ""),
                "record": result["record"],
                "resolution": {
                    "status": result["status"],
                    "top_score": result["top_score"],
                    "top_margin": result["top_margin"],
                    "resolver_version": result.get("resolver_version", ""),
                    "query_notes": result.get("query_notes", []),
                    "candidates": result["candidates"],
                },
            }
            if result["status"] == "invalid":
                row_contract["invalid_reason"] = result.get("invalid_reason", "")
            self.apply_trait_identity_resolution_policy(row_contract)
            row_contracts.append(row_contract)

        self.apply_pathway_context_rerank(row_contracts)
        for row_contract in row_contracts:
            self.apply_trait_identity_resolution_policy(row_contract)
            row_contract["identity_resolution_v2"] = self.identity_resolution_v2_for_row(row_contract)

        matched = []
        ambiguous = []
        unmatched = []
        invalid = []
        for row_contract in row_contracts:
            status = row_contract["resolution"]["status"]
            candidates = row_contract["resolution"].get("candidates") or []
            if status == "matched":
                row_contract["metabolite_uid"] = candidates[0]["entity_uid"]
                matched.append(row_contract)
            elif status == "ambiguous":
                ambiguous.append(row_contract)
            elif status == "invalid":
                invalid.append(row_contract)
            else:
                unmatched.append(row_contract)
        expanded_candidates = self.expanded_candidate_rows(ambiguous, unmatched)
        expanded_summary = self.expanded_candidate_summary(matched, expanded_candidates, ambiguous, unmatched, invalid)
        payload = {
            "analysis_mode": input_normalization.get("analysis_mode", "metabolite_table"),
            "input_normalization": input_normalization,
            "resolver_version": "compound_resolver_v2" if self.compound_index_available() else "resolver_index_fallback",
            "input_count": len(records),
            "summary": {
                "matched": len(matched),
                "ambiguous": len(ambiguous),
                "unmatched": len(unmatched),
                "invalid": len(invalid),
            },
            "matched": matched,
            "ambiguous": ambiguous,
            "unmatched": unmatched,
            "invalid": invalid,
            "expanded_candidates": expanded_candidates,
            "expanded_summary": expanded_summary,
            "identity_resolution_v2_summary": self.identity_resolution_v2_summary(row_contracts),
            "supported_identifier_columns": [
                "HMDB",
                "ChEBI",
                "PubChem CID",
                "KEGG",
                "InChIKey",
                "name",
                "formula",
                "m/z + adduct + ppm_tolerance",
                "rt + rt_method + rt_tolerance",
                "ms2_peaks + ms2_mz_tolerance",
            ],
        }
        return self.response_envelope("/precheck/metabolites", request, payload)

    def get_node(self, node_uid: str) -> dict[str, Any] | None:
        if node_uid in self._node_cache:
            return self._node_cache[node_uid]
        path = self.graph_dir / "nodes.parquet"
        if not path.exists():
            self._node_cache[node_uid] = None
            return None
        rows = table_rows(path, NODE_COLUMNS, ds.field("node_uid") == node_uid)
        node = rows[0] if rows else None
        self._node_cache[node_uid] = node
        return node

    def get_nodes(self, node_uids: list[str]) -> list[dict[str, Any]]:
        unique = sorted(set(node_uids))
        missing = [uid for uid in unique if uid not in self._node_cache]
        if missing:
            path = self.graph_dir / "nodes.parquet"
            if path.exists():
                rows = table_rows(path, NODE_COLUMNS, ds.field("node_uid").isin(missing))
                for row in rows:
                    self._node_cache[row["node_uid"]] = row
            for uid in missing:
                self._node_cache.setdefault(uid, None)
        return [self._node_cache[uid] for uid in unique if self._node_cache.get(uid)]

    def get_nodes_by_idx(self, node_idxs: list[int]) -> list[dict[str, Any]]:
        unique = sorted({int(idx) for idx in node_idxs})
        missing = [idx for idx in unique if idx not in self._node_idx_cache]
        if missing:
            path = self.graph_dir / "nodes.parquet"
            if path.exists():
                rows = table_rows(path, NODE_COLUMNS, ds.field("node_idx").isin(missing))
                for row in rows:
                    idx = int(row["node_idx"])
                    self._node_idx_cache[idx] = row
                    self._node_cache[row["node_uid"]] = row
            for idx in missing:
                self._node_idx_cache.setdefault(idx, None)
        return [self._node_idx_cache[idx] for idx in unique if self._node_idx_cache.get(idx)]

    def edge_type_index(self) -> dict[int, dict[str, Any]]:
        if self._edge_type_index_cache is not None:
            return self._edge_type_index_cache
        path = self.graph_dir / "edge_type_index.parquet"
        rows = table_rows(path) if path.exists() else []
        self._edge_type_index_cache = {int(row["edge_type_id"]): row for row in rows}
        return self._edge_type_index_cache

    def sparse_incident_edges(self, node_idx: int) -> list[dict[str, Any]]:
        node_idx = int(node_idx)
        return self.sparse_incident_edges_many([node_idx]).get(node_idx, [])

    def sparse_incident_edges_many(self, node_idxs: list[int]) -> dict[int, list[dict[str, Any]]]:
        requested = sorted({int(idx) for idx in node_idxs})
        missing = [idx for idx in requested if idx not in self._sparse_incident_cache]
        if not missing:
            return {idx: self._sparse_incident_cache.get(idx, []) for idx in requested}
        missing_set = set(missing)
        buckets: dict[int, list[dict[str, Any]]] = {idx: [] for idx in missing}
        path = self.graph_dir / "sparse_edges.parquet"
        type_index = self.edge_type_index()
        self.load_literature_support()
        by_edge = self._literature_support_by_edge_uid or {}
        if path.exists():
            rows = table_rows(
                path,
                None,
                ds.field("source_idx").isin(missing) | ds.field("target_idx").isin(missing),
            )
            for row in rows:
                raw_type_id = row.get("edge_type_id")
                type_row = type_index.get(int(raw_type_id) if raw_type_id is not None else -1, {})
                edge_type = type_row.get("edge_type", "")
                if edge_type not in ALLOWED_EXPLANATION_EDGE_TYPES:
                    continue
                source_idx = int(row.get("source_idx"))
                target_idx = int(row.get("target_idx"))
                incident_idxs = [idx for idx in (source_idx, target_idx) if idx in missing_set]
                if not incident_idxs:
                    continue
                weight = clamp_probability(row.get("weight"), 0.5)
                support_rows = [
                    support
                    for support in by_edge.get(str(row.get("edge_uid") or ""), [])
                    if support.get("support_class") in {"confirm", "support_direction"}
                ]
                if support_rows:
                    p_literature = independent_probability_union([float(support.get("p_literature") or 0.0) for support in support_rows])
                    weight = clamp_probability(1.0 - ((1.0 - weight) * (1.0 - p_literature)), weight)
                edge = {
                    "edge_uid": row.get("edge_uid", ""),
                    "edge_type": edge_type,
                    "subject_type": type_row.get("subject_type", ""),
                    "object_type": type_row.get("object_type", ""),
                    "source_idx": source_idx,
                    "target_idx": target_idx,
                    "weight": weight,
                    "p_final": weight,
                }
                for incident_idx in incident_idxs:
                    buckets[incident_idx].append(edge)
        node_rows = self.get_nodes_by_idx(missing)
        for node in node_rows:
            node_idx = int(node["node_idx"])
            current_uid = node.get("node_uid", "")
            for edge in self.literature_overlay_edges_for_entity(current_uid):
                subject = self.get_node(edge.get("subject_uid", ""))
                obj = self.get_node(edge.get("object_uid", ""))
                if not subject or not obj or subject.get("node_idx") is None or obj.get("node_idx") is None:
                    continue
                buckets.setdefault(node_idx, []).append(
                    {
                        "edge_uid": edge.get("edge_uid", ""),
                        "edge_type": edge.get("edge_type", ""),
                        "subject_type": edge.get("subject_type", ""),
                        "object_type": edge.get("object_type", ""),
                        "source_idx": int(subject["node_idx"]),
                        "target_idx": int(obj["node_idx"]),
                        "weight": clamp_probability(edge.get("p_final"), edge.get("weight", 0.5)),
                        "p_final": clamp_probability(edge.get("p_final"), edge.get("weight", 0.5)),
                        "literature_support": edge.get("literature_support", []),
                    }
                )
        for idx in missing:
            edges = buckets.get(idx, [])
            edges.sort(key=lambda row: (-row["weight"], row["edge_uid"]))
            self._sparse_incident_cache[idx] = edges
        return {idx: self._sparse_incident_cache.get(idx, []) for idx in requested}

    def entity_detail(self, entity_id: str, edge_limit: int = 25) -> dict[str, Any]:
        request = {"entity_id": entity_id, "edge_limit": edge_limit}
        node = self.get_node(entity_id)
        if not node:
            resolved = self.resolve(entity_id, limit=1)
            if resolved["status"] == "matched":
                node = self.get_node(resolved["candidates"][0]["entity_uid"])
        if not node:
            return self.response_envelope("/entity/{id}", request, {"status": "not_found", "entity": None, "edges": []})
        edges = self.incident_edges(node["node_uid"], limit=edge_limit)
        payload = {
            "status": "found",
            "entity": node,
            "edges": [self.edge_contract(edge) for edge in edges],
        }
        return self.response_envelope("/entity/{id}", request, payload)

    def pubchem(self, cid: str) -> dict[str, Any]:
        request = {"cid": str(cid)}
        cid = normalize_id_token(str(cid)).removeprefix("pubchem:").removeprefix("cid:")
        properties = None
        prop_path = self.pubchem_dir / "cid_properties.parquet"
        if prop_path.exists():
            rows = table_rows(prop_path, None, ds.field("pubchem_cid") == cid)
            properties = rows[0] if rows else None
        missing = None
        missing_path = self.pubchem_dir / "cid_missing.parquet"
        if properties is None and missing_path.exists():
            rows = table_rows(missing_path, None, ds.field("pubchem_cid") == cid)
            missing = rows[0] if rows else None
        payload = {
            "status": "found" if properties else "missing" if missing else "not_found",
            "cid": cid,
            "properties": properties,
            "missing": missing,
        }
        return self.response_envelope("/pubchem/{cid}", request, payload)

    def literature_manifest(self) -> dict[str, Any]:
        return read_json_file(self.literature_dir / "literature_evidence_manifest.json")

    def literature_support_contract(self, row: dict[str, Any]) -> dict[str, Any]:
        score_components = safe_metadata(row.get("score_components_json"))
        return {
            "support_uid": row.get("support_uid", ""),
            "subject_uid": row.get("subject_uid", ""),
            "subject_type": row.get("subject_type", ""),
            "predicate": row.get("predicate", ""),
            "object_uid": row.get("object_uid", ""),
            "object_type": row.get("object_type", ""),
            "polarity_set": parse_list(row.get("polarity_set")),
            "support_class": row.get("support_class", ""),
            "supported_existing_edge_uids": parse_list(row.get("supported_existing_edge_uids")),
            "evidence_relation_uids": parse_list(row.get("evidence_relation_uids")),
            "sentence_uids": parse_list(row.get("sentence_uids")),
            "pmids": parse_list(row.get("pmids")),
            "pmcids": parse_list(row.get("pmcids")),
            "evidence_sentence_count": int(row.get("evidence_sentence_count") or 0),
            "distinct_article_count": int(row.get("distinct_article_count") or 0),
            "p_literature": round(float(row.get("p_literature") or 0.0), 6),
            "raw_score_max": round(float(row.get("raw_score_max") or 0.0), 6),
            "calibrated_prob_max": round(float(row.get("calibrated_prob_max") or 0.0), 6),
            "score_components": score_components,
            "license_id": row.get("license_id", ""),
            "source_release": row.get("source_release", ""),
            "parser_hash": row.get("parser_hash", ""),
            "config_hash": row.get("config_hash", ""),
        }

    def load_literature_support(self) -> None:
        if self._literature_support_rows is not None:
            return
        support_path = self.literature_dir / "literature_edge_support.parquet"
        rows: list[dict[str, Any]] = []
        by_edge: dict[str, list[dict[str, Any]]] = defaultdict(list)
        by_triple: dict[tuple[str, str, str], list[dict[str, Any]]] = defaultdict(list)
        by_entity: dict[str, list[dict[str, Any]]] = defaultdict(list)
        if support_path.exists():
            for row in table_rows(support_path):
                contract = self.literature_support_contract(row)
                rows.append(contract)
                by_triple[(contract["subject_uid"], contract["predicate"], contract["object_uid"])].append(contract)
                by_entity[contract["subject_uid"]].append(contract)
                by_entity[contract["object_uid"]].append(contract)
                for edge_uid in contract.get("supported_existing_edge_uids", []):
                    by_edge[edge_uid].append(contract)
        rows.sort(key=lambda row: (-float(row.get("p_literature") or 0.0), row.get("support_uid", "")))
        for bucket in by_edge.values():
            bucket.sort(key=lambda row: (-float(row.get("p_literature") or 0.0), row.get("support_uid", "")))
        for bucket in by_triple.values():
            bucket.sort(key=lambda row: (-float(row.get("p_literature") or 0.0), row.get("support_uid", "")))
        for bucket in by_entity.values():
            bucket.sort(key=lambda row: (-float(row.get("p_literature") or 0.0), row.get("support_uid", "")))
        self._literature_support_rows = rows
        self._literature_support_by_edge_uid = dict(by_edge)
        self._literature_support_by_triple = dict(by_triple)
        self._literature_support_by_entity = dict(by_entity)

    def literature_support_for_edge(self, edge: dict[str, Any], supportive_only: bool = False) -> list[dict[str, Any]]:
        self.load_literature_support()
        by_edge = self._literature_support_by_edge_uid or {}
        by_triple = self._literature_support_by_triple or {}
        rows = list(by_edge.get(str(edge.get("edge_uid") or ""), []))
        literature_predicate = GRAPH_EDGE_TO_LITERATURE_PREDICATE.get(str(edge.get("edge_type") or ""))
        if literature_predicate:
            triple_key = (str(edge.get("subject_uid") or ""), literature_predicate, str(edge.get("object_uid") or ""))
            known = {row.get("support_uid") for row in rows}
            rows.extend(row for row in by_triple.get(triple_key, []) if row.get("support_uid") not in known)
        if supportive_only:
            rows = [row for row in rows if row.get("support_class") in {"confirm", "support_direction"}]
        return rows

    def literature_overlay_support_rows_for_entity(self, entity_uid: str) -> list[dict[str, Any]]:
        self.load_literature_support()
        rows = []
        for row in (self._literature_support_by_entity or {}).get(entity_uid, []):
            if row.get("subject_uid") != entity_uid:
                continue
            if row.get("support_class") != "novel_candidate":
                continue
            if row.get("predicate") not in ALLOWED_EXPLANATION_EDGE_TYPES:
                continue
            if float(row.get("p_literature") or 0.0) < self.config.min_literature_overlay_prob:
                continue
            if not row.get("subject_uid") or not row.get("object_uid"):
                continue
            if row.get("supported_existing_edge_uids"):
                continue
            rows.append(row)
        rows.sort(key=lambda row: (-float(row.get("p_literature") or 0.0), row.get("support_uid", "")))
        return rows[: min(self.config.max_edges_per_node, 16)]

    def literature_overlay_edge_from_support(self, support: dict[str, Any]) -> dict[str, Any] | None:
        subject_uid = str(support.get("subject_uid") or "")
        object_uid = str(support.get("object_uid") or "")
        if not subject_uid or not object_uid:
            return None
        subject = self.get_node(subject_uid)
        obj = self.get_node(object_uid)
        if not subject or not obj:
            return None
        support_uid = support.get("support_uid", "") or f"litsup_{content_hash([subject_uid, support.get('predicate', ''), object_uid])[:20]}"
        metadata = {
            "support_uid": support_uid,
            "support_class": support.get("support_class", ""),
            "polarity_set": parse_list(support.get("polarity_set")),
            "evidence_sentence_count": support.get("evidence_sentence_count", 0),
            "distinct_article_count": support.get("distinct_article_count", 0),
            "literature_overlay": True,
        }
        edge = {
            "edge_uid": f"literature_edge_{content_hash(support_uid)[:20]}",
            "edge_type": support.get("predicate", ""),
            "subject_uid": subject_uid,
            "subject_type": support.get("subject_type", "") or subject.get("node_type", ""),
            "predicate": support.get("predicate", ""),
            "object_uid": object_uid,
            "object_type": support.get("object_type", "") or obj.get("node_type", ""),
            "weight": round(float(support.get("p_literature") or 0.0), 6),
            "source_table": "literature_edge_support",
            "source_name": "literature_evidence",
            "source_record_id": support_uid,
            "evidence_level": "literature_overlay",
            "source_release": support.get("source_release", ""),
            "license_id": support.get("license_id", ""),
            "parser_hash": support.get("parser_hash", ""),
            "metadata_json": stable_json(metadata),
            "metadata": metadata,
            "literature_support": [support],
        }
        p_final, components = edge_probability(edge)
        edge["score_components"] = components
        edge["p_final"] = p_final
        edge["cost"] = round(-math.log(p_final), 6)
        return edge

    def literature_overlay_edges_for_entity(self, entity_uid: str) -> list[dict[str, Any]]:
        edges = []
        seen: set[str] = set()
        support_rows = self.literature_overlay_support_rows_for_entity(entity_uid)
        preload_uids = sorted(
            {
                uid
                for support in support_rows
                for uid in (support.get("subject_uid"), support.get("object_uid"))
                if uid
            }
        )
        if preload_uids:
            self.get_nodes(preload_uids)
        for support in support_rows:
            edge = self.literature_overlay_edge_from_support(support)
            if not edge or edge["edge_uid"] in seen:
                continue
            seen.add(edge["edge_uid"])
            edges.append(edge)
        edges.sort(key=lambda row: (-float(row.get("p_final") or 0.0), row.get("edge_uid", "")))
        return edges

    def incident_edges(self, node_uid: str, limit: int | None = None) -> list[dict[str, Any]]:
        if node_uid not in self._incident_edge_cache:
            edge_path = self.graph_dir / "edges.parquet"
            if not edge_path.exists():
                rows = []
            else:
                filt = (ds.field("subject_uid") == node_uid) | (ds.field("object_uid") == node_uid)
                rows = table_rows(edge_path, EDGE_COLUMNS, filt)
                rows = [self.annotate_edge(row) for row in rows]
            rows.sort(key=lambda row: (-row["score_components"]["p_final"], row["edge_uid"]))
            self._incident_edge_cache[node_uid] = rows
        rows = self._incident_edge_cache[node_uid]
        return rows[:limit] if limit else rows

    def annotate_edge(self, edge: dict[str, Any]) -> dict[str, Any]:
        if "score_components" in edge:
            return edge
        edge = dict(edge)
        edge["metadata"] = safe_metadata(edge.get("metadata_json"))
        edge["literature_support"] = self.literature_support_for_edge(edge, supportive_only=True)
        p_final, components = edge_probability(edge)
        edge["score_components"] = components
        edge["p_final"] = p_final
        edge["cost"] = round(-math.log(p_final), 6)
        return edge

    def edge_contract(self, edge: dict[str, Any]) -> dict[str, Any]:
        edge = self.annotate_edge(edge)
        return {
            "edge_uid": edge.get("edge_uid", ""),
            "edge_type": edge.get("edge_type", ""),
            "subject_uid": edge.get("subject_uid", ""),
            "subject_type": edge.get("subject_type", ""),
            "predicate": edge.get("predicate", ""),
            "object_uid": edge.get("object_uid", ""),
            "object_type": edge.get("object_type", ""),
            "source_table": edge.get("source_table", ""),
            "source_name": edge.get("source_name", ""),
            "source_record_id": edge.get("source_record_id", ""),
            "evidence_level": edge.get("evidence_level", ""),
            "source_release": edge.get("source_release", ""),
            "license_id": edge.get("license_id", ""),
            "parser_hash": edge.get("parser_hash", ""),
            "weight": edge.get("weight", 0.0),
            "p_final": edge.get("p_final"),
            "cost": edge.get("cost"),
            "score_components": edge.get("score_components"),
            "literature_support": edge.get("literature_support", []),
            "metadata": edge.get("metadata", {}),
        }

    def subgraph(
        self,
        seed_ids: list[str],
        max_hops: int = 1,
        edge_types: list[str] | None = None,
        max_edges: int | None = None,
    ) -> dict[str, Any]:
        max_edges = max_edges or self.config.max_subgraph_edges
        request = {"seed_ids": seed_ids, "max_hops": max_hops, "edge_types": edge_types, "max_edges": max_edges}
        seen_nodes = set(seed_ids)
        seen_edges: dict[str, dict[str, Any]] = {}
        frontier = set(seed_ids)
        allowed = set(edge_types or [])
        for _hop in range(max_hops):
            if not frontier or len(seen_edges) >= max_edges:
                break
            next_frontier = set()
            for node_uid in sorted(frontier):
                for edge in self.incident_edges(node_uid, limit=self.config.max_edges_per_node):
                    if allowed and edge.get("edge_type") not in allowed:
                        continue
                    seen_edges.setdefault(edge["edge_uid"], edge)
                    for endpoint in (edge["subject_uid"], edge["object_uid"]):
                        if endpoint not in seen_nodes:
                            seen_nodes.add(endpoint)
                            next_frontier.add(endpoint)
                    if len(seen_edges) >= max_edges:
                        break
                if len(seen_edges) >= max_edges:
                    break
            frontier = next_frontier
        nodes = self.get_nodes(list(seen_nodes))
        payload = {
            "nodes": sorted(nodes, key=lambda row: (row["node_type"], row["node_uid"])),
            "edges": [self.edge_contract(edge) for edge in sorted(seen_edges.values(), key=lambda row: row["edge_uid"])],
        }
        return self.response_envelope("/subgraph", request, payload)

    def relation_candidate_contract(self, row: dict[str, Any]) -> dict[str, Any]:
        return {
            "relation_uid": row.get("relation_uid", ""),
            "subject_uid": row.get("subject_uid", ""),
            "subject_type": row.get("subject_type", ""),
            "predicate": row.get("predicate", ""),
            "object_uid": row.get("object_uid", ""),
            "object_type": row.get("object_type", ""),
            "polarity": row.get("polarity", ""),
            "cancer_context": row.get("cancer_context", ""),
            "sentence_uid": row.get("sentence_uid", ""),
            "pmid": row.get("pmid", ""),
            "pmcid": row.get("pmcid", ""),
            "section": row.get("section", ""),
            "trigger_phrase": row.get("trigger_phrase", ""),
            "extraction_rule_id": row.get("extraction_rule_id", ""),
            "raw_score": round(float(row.get("raw_score") or 0.0), 6),
            "calibrated_prob": round(float(row.get("calibrated_prob") or 0.0), 6),
            "score_components": safe_metadata(row.get("score_components_json")),
            "license_id": row.get("license_id", ""),
            "source_release": row.get("source_release", ""),
            "parser_hash": row.get("parser_hash", ""),
            "config_hash": row.get("config_hash", ""),
        }

    def query_relation_candidates(
        self,
        relation_uids: list[str] | None = None,
        subject_uid: str = "",
        predicate: str = "",
        object_uid: str = "",
        entity_uid: str = "",
        sentence_uid: str = "",
        limit: int | None = None,
    ) -> list[dict[str, Any]]:
        path = self.literature_dir / "relation_candidates.parquet"
        if not path.exists():
            return []
        filt: ds.Expression | None = None

        def add_filter(new_filter: ds.Expression) -> None:
            nonlocal filt
            filt = new_filter if filt is None else filt & new_filter

        if relation_uids:
            add_filter(ds.field("relation_uid").isin(sorted(set(relation_uids))))
        if subject_uid:
            add_filter(ds.field("subject_uid") == subject_uid)
        if predicate:
            add_filter(ds.field("predicate") == predicate)
        if object_uid:
            add_filter(ds.field("object_uid") == object_uid)
        if entity_uid:
            add_filter((ds.field("subject_uid") == entity_uid) | (ds.field("object_uid") == entity_uid))
        if sentence_uid:
            add_filter(ds.field("sentence_uid") == sentence_uid)
        if filt is None:
            return []
        rows = [self.relation_candidate_contract(row) for row in table_rows(path, None, filt)]
        rows.sort(key=lambda row: (-float(row.get("calibrated_prob") or 0.0), row.get("relation_uid", "")))
        return rows[:limit] if limit else rows

    def query_sentence_mentions(self, sentence_uids: list[str]) -> list[dict[str, Any]]:
        path = self.literature_dir / "sentence_mentions.parquet"
        wanted = sorted(set(sentence_uids))
        if not path.exists() or not wanted:
            return []
        rows = table_rows(path, None, ds.field("sentence_uid").isin(wanted))
        rows.sort(key=lambda row: (row.get("sentence_uid", ""), int(row.get("start_offset") or 0), row.get("entity_type", ""), row.get("entity_uid", "")))
        return rows

    def query_sentences(self, sentence_uids: list[str]) -> list[dict[str, Any]]:
        path = self.normalized_dir / "sentences.parquet"
        wanted = sorted(set(sentence_uids))
        if not path.exists() or not wanted:
            return []
        rows = table_rows(
            path,
            ["sentence_uid", "article_uid", "pmid", "pmcid", "section", "sentence_text", "text_hash", "start_offset", "end_offset", "source_release", "parser_hash"],
            ds.field("sentence_uid").isin(wanted),
        )
        rows.sort(key=lambda row: row.get("sentence_uid", ""))
        return rows

    def evidence(
        self,
        edge_uid: str = "",
        subject_uid: str = "",
        predicate: str = "",
        object_uid: str = "",
        entity_uid: str = "",
        sentence_uid: str = "",
        relation_uid: str = "",
        limit: int | None = None,
    ) -> dict[str, Any]:
        limit = min(100, max(1, int(limit or self.config.max_evidence_items)))
        request = {
            "edge_uid": edge_uid,
            "subject_uid": subject_uid,
            "predicate": predicate,
            "object_uid": object_uid,
            "entity_uid": entity_uid,
            "sentence_uid": sentence_uid,
            "relation_uid": relation_uid,
            "limit": limit,
        }
        manifest = self.literature_manifest()
        if not self.literature_dir.exists():
            return self.response_envelope(
                "/evidence",
                request,
                {
                    "status": "missing_literature_evidence",
                    "manifest": manifest,
                    "support": [],
                    "relation_candidates": [],
                    "sentences": [],
                    "mentions": [],
                },
            )

        self.load_literature_support()
        support_rows: list[dict[str, Any]] = []
        all_support = self._literature_support_rows or []
        by_edge = self._literature_support_by_edge_uid or {}
        by_triple = self._literature_support_by_triple or {}
        if edge_uid:
            support_rows = list(by_edge.get(edge_uid, []))
        elif subject_uid and object_uid and predicate:
            support_rows = list(by_triple.get((subject_uid, predicate, object_uid), []))
        elif subject_uid and object_uid:
            support_rows = [row for row in all_support if row.get("subject_uid") == subject_uid and row.get("object_uid") == object_uid]
        elif entity_uid:
            support_rows = [row for row in all_support if row.get("subject_uid") == entity_uid or row.get("object_uid") == entity_uid]
        support_rows.sort(key=lambda row: (-float(row.get("p_literature") or 0.0), row.get("support_uid", "")))
        support_rows = support_rows[:limit]

        relation_uids = [relation_uid] if relation_uid else []
        for row in support_rows:
            relation_uids.extend(parse_list(row.get("evidence_relation_uids")))
        relation_uids = sorted(set(uid for uid in relation_uids if uid))
        relation_limit = max(limit, min(250, limit * 5))
        if relation_uids:
            relation_rows = self.query_relation_candidates(relation_uids=relation_uids, limit=relation_limit)
        else:
            relation_rows = self.query_relation_candidates(
                subject_uid=subject_uid,
                predicate=predicate,
                object_uid=object_uid,
                entity_uid=entity_uid,
                sentence_uid=sentence_uid,
                limit=relation_limit,
            )
        relation_rows = relation_rows[:limit]

        sentence_uids = [sentence_uid] if sentence_uid else []
        for row in support_rows:
            sentence_uids.extend(parse_list(row.get("sentence_uids")))
        sentence_uids.extend(row.get("sentence_uid", "") for row in relation_rows if row.get("sentence_uid"))
        sentence_uids = sorted(set(uid for uid in sentence_uids if uid))
        sentences = self.query_sentences(sentence_uids)
        mentions = self.query_sentence_mentions(sentence_uids)
        payload = {
            "status": "found" if support_rows or relation_rows or sentences or mentions else "not_found",
            "manifest": {
                "path": str(self.literature_dir / "literature_evidence_manifest.json"),
                "manifest_hash": manifest.get("manifest_hash", ""),
                "builder_version": manifest.get("builder_version", ""),
                "metrics": manifest.get("metrics", {}),
            },
            "support": support_rows,
            "relation_candidates": relation_rows,
            "sentences": sentences,
            "mentions": mentions,
            "score_notes": {
                "p_literature": "p_literature(edge)=1-prod(1-max_sentence_probability_per_article)",
                "p_final": "p_final(e)=1-(1-p_curated)(1-p_literature)(1-p_topology)(1-p_user)",
                "evidence_role": "Literature evidence is a validation and calibration overlay; it does not overwrite curated graph facts.",
            },
        }
        return self.response_envelope("/evidence", request, payload)

    def build_pathway_index(self) -> dict[str, Any]:
        if self._pathway_index is not None:
            return self._pathway_index
        metabolite_to_pathways: dict[str, set[str]] = defaultdict(set)
        pathway_to_metabolites: dict[str, set[str]] = defaultdict(set)
        evidence: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)

        direct_path = self.normalized_dir / "metabolite_pathway_edges.parquet"
        if direct_path.exists():
            rows = table_rows(
                direct_path,
                [
                    "edge_uid",
                    "metabolite_uid",
                    "pathway_uid",
                    "source_name",
                    "source_record_id",
                    "evidence_level",
                    "license_id",
                    "source_release",
                ],
            )
            for row in rows:
                metabolite_uid = row.get("metabolite_uid")
                pathway_uid = row.get("pathway_uid")
                if not metabolite_uid or not pathway_uid:
                    continue
                metabolite_to_pathways[metabolite_uid].add(pathway_uid)
                pathway_to_metabolites[pathway_uid].add(metabolite_uid)
                evidence[(metabolite_uid, pathway_uid)].append(row)

        reaction_to_pathway: dict[str, str] = {}
        reactions_path = self.normalized_dir / "reactions.parquet"
        if reactions_path.exists():
            for row in table_rows(reactions_path, ["reaction_uid", "pathway_uid"]):
                if row.get("reaction_uid") and row.get("pathway_uid"):
                    reaction_to_pathway[row["reaction_uid"]] = row["pathway_uid"]
        participants_path = self.normalized_dir / "reaction_participants.parquet"
        if participants_path.exists() and reaction_to_pathway:
            rows = table_rows(
                participants_path,
                [
                    "edge_uid",
                    "reaction_uid",
                    "participant_uid",
                    "participant_type",
                    "source_name",
                    "source_record_id",
                    "license_id",
                    "source_release",
                ],
            )
            for row in rows:
                if str(row.get("participant_type") or "").lower() != "metabolite":
                    continue
                metabolite_uid = row.get("participant_uid")
                pathway_uid = reaction_to_pathway.get(row.get("reaction_uid"))
                if not metabolite_uid or not pathway_uid:
                    continue
                metabolite_to_pathways[metabolite_uid].add(pathway_uid)
                pathway_to_metabolites[pathway_uid].add(metabolite_uid)
                evidence[(metabolite_uid, pathway_uid)].append({**row, "evidence_level": "curated_reaction_participant"})

        pathway_names: dict[str, dict[str, Any]] = {}
        pathways_path = self.normalized_dir / "pathways.parquet"
        if pathways_path.exists():
            for row in table_rows(pathways_path, ["pathway_uid", "name", "primary_external_id", "source_name", "species"]):
                pathway_names[row["pathway_uid"]] = row

        self._pathway_index = {
            "metabolite_to_pathways": metabolite_to_pathways,
            "pathway_to_metabolites": pathway_to_metabolites,
            "evidence": evidence,
            "pathway_names": pathway_names,
            "universe_metabolites": set(metabolite_to_pathways),
        }
        return self._pathway_index

    def pathway_enrichment(self, metabolite_uids: list[str], seed_presence_weights: dict[str, float] | None = None) -> list[dict[str, Any]]:
        index = self.build_pathway_index()
        seeds = sorted(set(metabolite_uids) & index["universe_metabolites"])
        if not seeds:
            return []
        seed_presence_weights = seed_presence_weights or {}
        universe_size = len(index["universe_metabolites"])
        draws = len(seeds)
        weighted_draws = sum(float(seed_presence_weights.get(uid, 1.0) or 0.0) for uid in seeds) or float(draws)
        pathway_hits: dict[str, set[str]] = defaultdict(set)
        for uid in seeds:
            for pathway_uid in index["metabolite_to_pathways"].get(uid, set()):
                pathway_hits[pathway_uid].add(uid)
        rows = []
        for pathway_uid, hit_uids in pathway_hits.items():
            pathway_metabolites = index["pathway_to_metabolites"][pathway_uid]
            overlap = len(hit_uids)
            weighted_overlap = sum(float(seed_presence_weights.get(uid, 1.0) or 0.0) for uid in hit_uids)
            pathway_size = len(pathway_metabolites)
            p_value = hypergeom_sf(overlap, universe_size, pathway_size, draws)
            pathway = index["pathway_names"].get(pathway_uid, {})
            evidence_refs = []
            for metabolite_uid in sorted(hit_uids):
                for edge_ref in index["evidence"].get((metabolite_uid, pathway_uid), [])[:3]:
                    evidence_refs.append(
                        {
                            "metabolite_uid": metabolite_uid,
                            "edge_uid": edge_ref.get("edge_uid", ""),
                            "source_name": edge_ref.get("source_name", ""),
                            "source_record_id": edge_ref.get("source_record_id", ""),
                            "evidence_level": edge_ref.get("evidence_level", ""),
                            "license_id": edge_ref.get("license_id", ""),
                            "source_release": edge_ref.get("source_release", ""),
                        }
                    )
            rows.append(
                {
                    "pathway_uid": pathway_uid,
                    "name": pathway.get("name", ""),
                    "primary_external_id": pathway.get("primary_external_id", ""),
                    "source_name": pathway.get("source_name", ""),
                    "species": pathway.get("species", ""),
                    "overlap_count": overlap,
                    "weighted_overlap": round(weighted_overlap, 6),
                    "input_count": draws,
                    "weighted_input_count": round(weighted_draws, 6),
                    "pathway_size": pathway_size,
                    "universe_size": universe_size,
                    "coverage": round(overlap / max(1, draws), 6),
                    "weighted_coverage": round(weighted_overlap / max(1e-12, weighted_draws), 6),
                    "p_value": round(p_value, 12),
                    "matched_metabolite_uids": sorted(hit_uids),
                    "evidence_refs": evidence_refs[:10],
                }
            )
        benjamini_hochberg(rows)
        for row in rows:
            fdr = max(float(row.get("fdr", 1.0)), 1e-300)
            row["score"] = round(-math.log10(fdr) * row["coverage"], 6)
            weighted_coverage = float(row.get("weighted_coverage", row["coverage"]) or 0.0)
            row["weighted_score"] = round(-math.log10(fdr) * weighted_coverage, 6)
            row["score_components"] = {
                "overlap_count": row["overlap_count"],
                "weighted_overlap": row.get("weighted_overlap", row["overlap_count"]),
                "coverage": row["coverage"],
                "weighted_coverage": row.get("weighted_coverage", row["coverage"]),
                "p_value": row["p_value"],
                "fdr": row["fdr"],
                "seed_weight_policy": "strict_identity=1.0; expanded seed classes use configured fractional weights split across candidates",
            }
            row["score"] = row["weighted_score"]
        return sorted(rows, key=lambda row: (-row["score"], row["fdr"], row["pathway_uid"]))

    def matched_analysis_features(self, matched_rows: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], dict[str, list[dict[str, Any]]]]:
        features = []
        by_uid: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for row in matched_rows:
            metabolite_uid = row.get("metabolite_uid", "")
            feature = analysis_feature_from_record(row.get("record", {}), row.get("input_id", ""))
            feature["metabolite_uid"] = metabolite_uid
            candidate = (row.get("resolution", {}).get("candidates") or [{}])[0]
            feature["display_name"] = candidate.get("display_name", "")
            feature["primary_external_id"] = candidate.get("primary_external_id", "")
            expanded = row.get("expanded_candidate") or {}
            seed_track = str(expanded.get("track") or "strict")
            seed_class = str(expanded.get("seed_class") or "strict_identity")
            weight_multiplier = maybe_number(expanded.get("weight_multiplier"))
            if weight_multiplier is None:
                weight_multiplier = 1.0
            weight_multiplier = max(0.0, float(weight_multiplier))
            raw_p_user = float(feature.get("p_user", 0.0) or 0.0)
            feature["p_user_raw"] = round(raw_p_user, 6)
            feature["p_user"] = round(raw_p_user * weight_multiplier, 6)
            feature["seed_track"] = seed_track
            feature["seed_class"] = seed_class
            feature["seed_weight_multiplier"] = round(weight_multiplier, 6)
            feature["seed_weight"] = round((1.0 + raw_p_user) * weight_multiplier, 6)
            feature["seed_weight_policy"] = "strict_identity=1.0; expanded classes use fractional multipliers split across candidates"
            if expanded:
                feature["expanded_candidate"] = dict(expanded)
            ratio_component = expanded.get("ratio_component") if expanded else {}
            if seed_class == "ratio_component" and ratio_component:
                signed_ratio_direction = str(ratio_component.get("signed_ratio_direction") or "unknown")
                feature["ratio_trait"] = dict(expanded.get("ratio_trait") or {})
                feature["ratio_component"] = dict(ratio_component)
                feature["ratio_effect_semantics"] = "relative_ratio_consistent_direction_not_measured_abundance"
                feature["ratio_original_direction"] = feature.get("direction", "unknown")
                feature["ratio_original_log2_fold_change"] = feature.get("log2_fold_change")
                if signed_ratio_direction == "opposite_to_ratio":
                    if feature.get("log2_fold_change") is not None:
                        feature["log2_fold_change"] = round(-float(feature["log2_fold_change"]), 6)
                    feature["direction"] = opposite_direction(feature.get("direction"))
                elif signed_ratio_direction == "same_as_ratio":
                    feature["direction"] = feature.get("direction", "unknown")
                else:
                    feature["direction"] = "unknown"
                feature.setdefault("score_components", {})
                feature["score_components"]["ratio_direction_basis"] = signed_ratio_direction
            class_seed = expanded.get("class_seed") if expanded else {}
            if seed_class == "class_or_pool" and class_seed:
                feature["class_seed"] = dict(class_seed)
                feature["class_effect_semantics"] = "class_or_pool_signal_not_unique_structural_isomer"
            identity_cluster = expanded.get("identity_cluster") if expanded else {}
            if seed_class == "soft_identity" and identity_cluster:
                feature["identity_cluster"] = dict(identity_cluster)
                feature["identity_effect_semantics"] = "identity_cluster_or_consensus_only"
            normalized_record = self.normalized_record_keys(row.get("record", {}))
            feature["input_name"] = str(
                first_present(normalized_record, ("name", "metabolite", "compound")) or feature["display_name"] or ""
            )
            feature["biological_entity_pools"] = self.biological_entity_pools_for_values(
                feature.get("input_name", ""),
                feature.get("display_name", ""),
                candidate.get("primary_external_id", ""),
            )
            feature.setdefault("score_components", {})
            feature["score_components"].update(
                {
                    "p_user_raw": round(raw_p_user, 6),
                    "p_user_weighted": feature["p_user"],
                    "seed_weight_multiplier": feature["seed_weight_multiplier"],
                    "seed_weight": feature["seed_weight"],
                }
            )
            features.append(feature)
            if metabolite_uid:
                by_uid[metabolite_uid].append(feature)
        return features, by_uid

    def feature_summary(self, features: list[dict[str, Any]]) -> dict[str, Any]:
        directions = {
            "up": sum(1 for feature in features if feature["direction"] == "up"),
            "down": sum(1 for feature in features if feature["direction"] == "down"),
            "unchanged": sum(1 for feature in features if feature["direction"] == "unchanged"),
            "unknown": sum(1 for feature in features if feature["direction"] == "unknown"),
        }
        seed_tracks = Counter(str(feature.get("seed_track") or "strict") for feature in features)
        seed_classes = Counter(str(feature.get("seed_class") or "strict_identity") for feature in features)
        weighted_seed_mass = sum(float(feature.get("seed_weight") or (1.0 + float(feature.get("p_user", 0.0) or 0.0))) for feature in features)
        return {
            "matched_records": len(features),
            "analysis_seed_records": len(features),
            "strict_seed_records": int(seed_tracks.get("strict", 0)),
            "expanded_seed_records": int(seed_tracks.get("expanded", 0)),
            "seed_tracks": dict(sorted(seed_tracks.items())),
            "seed_classes": dict(sorted(seed_classes.items())),
            "weighted_seed_mass": round(weighted_seed_mass, 6),
            "with_fold_change": sum(1 for feature in features if feature["log2_fold_change"] is not None),
            "with_p_value": sum(1 for feature in features if feature["p_value"] is not None or feature["adjusted_p_value"] is not None),
            "significant": sum(1 for feature in features if feature["significant"]),
            "directions": directions,
            "p_user": {
                "max": round(max((feature["p_user"] for feature in features), default=0.0), 6),
                "mean": round(sum(feature["p_user"] for feature in features) / max(1, len(features)), 6),
                "formula": "1-(1-p_value_support)(1-effect_support)",
            },
        }

    def seed_weights_from_features(self, features_by_uid: dict[str, list[dict[str, Any]]]) -> dict[str, float]:
        weights = {}
        for metabolite_uid in sorted(features_by_uid):
            weights[metabolite_uid] = round(
                sum(
                    float(feature.get("seed_weight"))
                    if feature.get("seed_weight") is not None
                    else 1.0 + float(feature.get("p_user", 0.0))
                    for feature in features_by_uid[metabolite_uid]
                ),
                6,
            )
        return weights

    def seed_presence_weights_from_features(self, features_by_uid: dict[str, list[dict[str, Any]]]) -> dict[str, float]:
        weights = {}
        for metabolite_uid, features in features_by_uid.items():
            value = sum(float(feature.get("seed_weight_multiplier") or 0.0) for feature in features)
            if value > 0.0:
                weights[metabolite_uid] = round(value, 6)
        return weights

    def directional_uids(self, features_by_uid: dict[str, list[dict[str, Any]]]) -> dict[str, list[str]]:
        grouped: dict[str, set[str]] = {"up": set(), "down": set(), "significant_up": set(), "significant_down": set()}
        for metabolite_uid, features in features_by_uid.items():
            for feature in features:
                direction = feature.get("direction")
                if direction in {"up", "down"}:
                    grouped[direction].add(metabolite_uid)
                    if feature.get("significant"):
                        grouped[f"significant_{direction}"].add(metabolite_uid)
        return {key: sorted(value) for key, value in grouped.items()}

    def attach_pathway_user_support(
        self,
        pathways: list[dict[str, Any]],
        features_by_uid: dict[str, list[dict[str, Any]]],
        input_records: list[Any] | None = None,
    ) -> list[dict[str, Any]]:
        theme_profile = self.input_theme_profile(features_by_uid, input_records=input_records)
        for row in pathways:
            counts = {"up": 0, "down": 0, "unchanged": 0, "unknown": 0, "significant": 0}
            p_user_sum = 0.0
            for metabolite_uid in row.get("matched_metabolite_uids", []):
                for feature in features_by_uid.get(metabolite_uid, []):
                    direction = feature.get("direction", "unknown")
                    counts[direction if direction in counts else "unknown"] += 1
                    counts["significant"] += 1 if feature.get("significant") else 0
                    p_user_sum += float(feature.get("p_user", 0.0))
            row["directional_support"] = {
                **counts,
                "p_user_sum": round(p_user_sum, 6),
            }
            row.setdefault("score_components", {})["user_support_sum"] = round(p_user_sum, 6)
            theme_boost, theme_hits = self.pathway_theme_boost(row, features_by_uid)
            row["directional_support"]["pathway_theme_hits"] = theme_hits
            row["score_components"]["pathway_theme_boost"] = theme_boost
            if theme_boost:
                row["score"] = round(float(row.get("score", 0.0) or 0.0) + theme_boost, 6)
            input_theme_boost, input_theme_hits = self.pathway_input_theme_boost(row, theme_profile)
            row["directional_support"]["input_theme_hits"] = input_theme_hits
            row["score_components"]["input_theme_boost"] = input_theme_boost
            if input_theme_boost:
                row["score"] = round(float(row.get("score", 0.0) or 0.0) + input_theme_boost, 6)
            lexical_boost, lexical_hits = self.pathway_lexical_input_boost(row, theme_profile)
            row["directional_support"]["input_lexical_hits"] = lexical_hits
            row["score_components"]["input_lexical_boost"] = lexical_boost
            if lexical_boost:
                row["score"] = round(float(row.get("score", 0.0) or 0.0) + lexical_boost, 6)
            pathway_penalty = self.pathway_context_penalty(row)
            row["score_components"]["pathway_context_penalty"] = pathway_penalty
            if pathway_penalty:
                row["score"] = round(max(0.0, float(row.get("score", 0.0) or 0.0) - pathway_penalty), 6)
            generic_penalty = self.generic_pathway_penalty(row)
            row["score_components"]["generic_pathway_penalty"] = generic_penalty
            if generic_penalty:
                row["score"] = round(max(0.0, float(row.get("score", 0.0) or 0.0) - generic_penalty), 6)
        return pathways

    def input_theme_profile(
        self,
        features_by_uid: dict[str, list[dict[str, Any]]],
        input_records: list[Any] | None = None,
    ) -> dict[str, dict[str, Any]]:
        profile: dict[str, dict[str, Any]] = {}
        lexical_terms: dict[str, dict[str, Any]] = {}

        def add_label(label: Any, support_uid: str, p_user: float = 0.0) -> None:
            normalized_label = normalize_lookup_key(label)
            if not normalized_label:
                return
            for theme_id, spec in BIOCHEMICAL_THEME_SPECS.items():
                matched_terms = [term for term in spec["input_terms"] if term in normalized_label]
                if not matched_terms:
                    continue
                theme_id = THEME_ID_ALIASES.get(theme_id, theme_id)
                row = profile.setdefault(
                    theme_id,
                    {"theme_id": theme_id, "metabolite_uids": set(), "input_terms": set(), "p_user_sum": 0.0},
                )
                if support_uid not in row["metabolite_uids"]:
                    row["metabolite_uids"].add(support_uid)
                    row["p_user_sum"] += p_user
                row["input_terms"].update(matched_terms)
            for term in self.biochemical_label_terms(normalized_label):
                row = lexical_terms.setdefault(term, {"term": term, "support_uids": set(), "p_user_sum": 0.0})
                if support_uid not in row["support_uids"]:
                    row["support_uids"].add(support_uid)
                    row["p_user_sum"] += p_user

        for metabolite_uid, features in features_by_uid.items():
            for feature in features:
                p_user = float(feature.get("p_user", 0.0) or 0.0)
                add_label(feature.get("input_name", ""), metabolite_uid, p_user)
                add_label(feature.get("display_name", ""), metabolite_uid, p_user)
        for index, record in enumerate(input_records or []):
            normalized_record = self.normalized_record_keys(record)
            label = first_present(normalized_record, ("name", "metabolite", "compound", "label"))
            add_label(label, f"input_row_{index}", 0.0)
        normalized = {}
        for theme_id, row in profile.items():
            support_count = len(row["metabolite_uids"])
            if support_count <= 0:
                continue
            normalized[theme_id] = {
                "theme_id": theme_id,
                "support_count": support_count,
                "input_terms": sorted(row["input_terms"]),
                "p_user_sum": round(float(row["p_user_sum"]), 6),
            }
        if lexical_terms:
            normalized["_lexical"] = {
                "theme_id": "_lexical",
                "terms": {
                    term: {
                        "support_count": len(row["support_uids"]),
                        "p_user_sum": round(float(row["p_user_sum"]), 6),
                    }
                    for term, row in sorted(lexical_terms.items())
                },
            }
        return normalized

    def biological_entity_pools_for_values(self, *values: Any) -> list[dict[str, Any]]:
        labels = [normalize_lookup_key(value) for value in values if value is not None]
        cache_key = tuple(label for label in labels if label)
        if cache_key in self._biological_entity_pool_cache:
            return [dict(row) for row in self._biological_entity_pool_cache[cache_key]]
        joined = " ".join(label for label in labels if label)
        if not joined:
            return []
        joined_spaced = joined.replace("-", " ")
        pools = []
        for pool_id, spec in sorted(BIOLOGICAL_ENTITY_POOLS.items()):
            matched_terms = sorted(
                {
                    term
                    for term in spec["terms"]
                    if normalized_text_has_term(joined, normalize_lookup_key(term))
                    or normalized_text_has_term(joined_spaced, normalize_lookup_key(term).replace("-", " "))
                }
            )
            if not matched_terms:
                continue
            pools.append(
                {
                    "pool_id": pool_id,
                    "display_name": spec["label"],
                    "matched_terms": matched_terms[:10],
                }
            )
        self._biological_entity_pool_cache[cache_key] = [dict(row) for row in pools]
        return pools

    def biological_interpretation_entities_for_resolution_row(self, row: dict[str, Any]) -> list[dict[str, Any]]:
        normalized_record = self.normalized_record_keys(row.get("record", {}))
        values = [
            row.get("query", ""),
            row.get("input_id", ""),
            first_present(normalized_record, ("name", "metabolite", "compound", "label")),
        ]
        for candidate in (row.get("resolution", {}).get("candidates") or [])[:5]:
            values.extend([candidate.get("display_name", ""), candidate.get("primary_external_id", "")])
        return self.biological_entity_pools_for_values(*values)

    def build_biological_entity_pool_summary(self, precheck: dict[str, Any]) -> list[dict[str, Any]]:
        summary: dict[str, dict[str, Any]] = {}
        for status in ("matched", "ambiguous", "unmatched", "invalid"):
            for row in precheck.get(status, []) or []:
                input_id = str(row.get("input_id") or row.get("query") or "")
                for pool in self.biological_interpretation_entities_for_resolution_row(row):
                    pool_id = pool["pool_id"]
                    entry = summary.setdefault(
                        pool_id,
                        {
                            "pool_id": pool_id,
                            "display_name": pool["display_name"],
                            "related_input_count": 0,
                            "matched_support_count": 0,
                            "ambiguous_support_count": 0,
                            "unmatched_support_count": 0,
                            "input_ids": set(),
                            "matched_terms": set(),
                        },
                    )
                    if input_id not in entry["input_ids"]:
                        entry["input_ids"].add(input_id)
                        entry["related_input_count"] += 1
                        if status == "matched":
                            entry["matched_support_count"] += 1
                        elif status == "ambiguous":
                            entry["ambiguous_support_count"] += 1
                        else:
                            entry["unmatched_support_count"] += 1
                    entry["matched_terms"].update(pool.get("matched_terms", []))
        rows = []
        for entry in summary.values():
            rows.append(
                {
                    "pool_id": entry["pool_id"],
                    "display_name": entry["display_name"],
                    "related_input_count": entry["related_input_count"],
                    "matched_support_count": entry["matched_support_count"],
                    "ambiguous_support_count": entry["ambiguous_support_count"],
                    "unmatched_support_count": entry["unmatched_support_count"],
                    "input_ids": sorted(entry["input_ids"])[:20],
                    "matched_terms": sorted(entry["matched_terms"])[:20],
                }
            )
        return sorted(rows, key=lambda row: (-row["related_input_count"], row["pool_id"]))

    def resolution_row_theme_hits(self, row: dict[str, Any]) -> dict[str, set[str]]:
        normalized_record = self.normalized_record_keys(row.get("record", {}))
        labels = [
            row.get("query", ""),
            first_present(normalized_record, ("name", "metabolite", "compound", "label")),
        ]
        for candidate in (row.get("resolution", {}).get("candidates") or [])[:5]:
            labels.append(candidate.get("display_name", ""))
        normalized_label = " ".join(normalize_lookup_key(label) for label in labels if label)
        normalized_label_spaced = normalized_label.replace("-", " ")
        pools = {pool["pool_id"] for pool in self.biological_entity_pools_for_values(*labels)}
        hits: dict[str, set[str]] = {}
        for theme_id, spec in BIOCHEMICAL_THEME_SPECS.items():
            matched_terms = {
                term
                for term in spec.get("input_terms", ())
                if normalized_text_has_term(normalized_label, normalize_lookup_key(term))
                or normalized_text_has_term(normalized_label_spaced, normalize_lookup_key(term).replace("-", " "))
            }
            matched_pools = pools & set(spec.get("pool_ids", ()))
            if not matched_terms and not matched_pools:
                continue
            canonical_theme_id = THEME_ID_ALIASES.get(theme_id, theme_id)
            theme_hits = hits.setdefault(canonical_theme_id, set())
            theme_hits.update(matched_terms)
            theme_hits.update(f"pool:{pool_id}" for pool_id in matched_pools)
        return hits

    def biological_pattern_consistency_for_theme(self, theme_id: str, term_directions: list[dict[str, Any]]) -> tuple[float | None, str]:
        observations: list[tuple[str, str]] = []
        for row in term_directions:
            direction = str(row.get("direction") or "unknown")
            if direction not in {"up", "down"}:
                continue
            for term in row.get("terms") or []:
                observations.append((normalize_lookup_key(term), direction))

        def has_any(terms: tuple[str, ...], direction: str) -> bool:
            return any(direction == observed_direction and any(term in observed_term for term in terms) for observed_term, observed_direction in observations)

        if theme_id == "glutathione_redox_stress":
            if has_any(("reduced glutathione", "gsh", "glutathione"), "down") and has_any(("oxidized glutathione", "gssg"), "up"):
                return 1.0, "redox_pair_consistent"
        if theme_id == "acylcarnitine_fatty_acid_oxidation_pressure":
            free_down = has_any(("carnitine", "acetylcarnitine"), "down")
            long_chain_up = has_any(("palmitoylcarnitine", "stearoylcarnitine", "long chain acylcarnitine", "long-chain acylcarnitine"), "up")
            if free_down and long_chain_up:
                return 1.0, "fatty_acid_oxidation_pressure_consistent"
        if theme_id == "arginine_no_metabolism":
            if has_any(("adma", "asymmetric dimethylarginine"), "up") and has_any(("arginine",), "down") and has_any(("nitrite",), "down"):
                return 1.0, "no_bioavailability_pressure_consistent"
        if theme_id == "idh_like_metabolic_pressure":
            if has_any(("2 hydroxyglutarate", "2-hydroxyglutarate", "hydroxyglutarate"), "up") and has_any(
                ("alpha ketoglutarate", "alpha-ketoglutarate", "2 oxoglutarate", "2-oxoglutarate", "akg"), "down"
            ):
                return 1.0, "idh_like_2hg_akg_pattern_consistent"
        return None, "not_pattern_evaluable"

    def build_theme_coverage(
        self,
        precheck: dict[str, Any],
        input_summary: dict[str, Any],
    ) -> dict[str, dict[str, Any]]:
        coverage: dict[str, dict[str, Any]] = {}
        for status in ("matched", "ambiguous", "unmatched", "invalid"):
            for row in precheck.get(status, []) or []:
                hits = self.resolution_row_theme_hits(row)
                if not hits:
                    continue
                feature = analysis_feature_from_record(row.get("record", {}), row.get("input_id", ""))
                direction = feature.get("direction", "unknown")
                input_id = str(row.get("input_id") or row.get("query") or "")
                for theme_id, matched_terms in hits.items():
                    entry = coverage.setdefault(
                        theme_id,
                        {
                            "theme_id": theme_id,
                            "related_input_count": 0,
                            "matched_support_count": 0,
                            "ambiguous_support_count": 0,
                            "unmatched_support_count": 0,
                            "significance_support_count": 0,
                            "input_terms": set(),
                            "input_ids": set(),
                            "directions": Counter(),
                            "term_directions": [],
                            "p_user_sum": 0.0,
                        },
                    )
                    if input_id in entry["input_ids"]:
                        entry["input_terms"].update(matched_terms)
                        continue
                    entry["input_ids"].add(input_id)
                    entry["related_input_count"] += 1
                    if status == "matched":
                        entry["matched_support_count"] += 1
                    elif status == "ambiguous":
                        entry["ambiguous_support_count"] += 1
                    else:
                        entry["unmatched_support_count"] += 1
                    if feature.get("significant"):
                        entry["significance_support_count"] += 1
                    entry["directions"][direction if direction in {"up", "down", "unchanged", "unknown"} else "unknown"] += 1
                    entry["term_directions"].append({"terms": sorted(matched_terms), "direction": direction})
                    entry["p_user_sum"] += float(feature.get("p_user", 0.0) or 0.0)
                    entry["input_terms"].update(matched_terms)

        result = {}
        confidence_weight = self.match_confidence_weight(input_summary)
        for theme_id, entry in coverage.items():
            direction_counts = dict(sorted(entry["directions"].items()))
            up = int(entry["directions"].get("up", 0))
            down = int(entry["directions"].get("down", 0))
            directed_total = up + down
            direction_consistency = round(max(up, down) / directed_total, 6) if directed_total else 0.0
            pattern_direction, biological_pattern = self.biological_pattern_consistency_for_theme(theme_id, entry.get("term_directions") or [])
            if pattern_direction is not None:
                direction_consistency = pattern_direction
            related = int(entry["related_input_count"])
            matched = int(entry["matched_support_count"])
            ambiguous = int(entry["ambiguous_support_count"])
            unmatched = int(entry["unmatched_support_count"])
            p_user_sum = float(entry["p_user_sum"] or 0.0)
            support_score = matched + 0.5 * ambiguous + 0.2 * unmatched
            score = min(0.9, 0.16 + 0.13 * support_score + 0.04 * min(5.0, p_user_sum))
            score *= max(0.55, confidence_weight)
            downgrade_reason = ""
            if matched == 0 and (ambiguous or unmatched):
                score = min(score, 0.35)
                downgrade_reason = "input_signal_without_stable_exact_match"
            elif ambiguous > matched:
                score = min(score, 0.44)
                downgrade_reason = "input_signal_mostly_ambiguous"
            elif unmatched and related >= 2:
                score = min(score, 0.62)
                downgrade_reason = "some_related_inputs_unmatched"
            elif direction_consistency and direction_consistency < 0.67:
                score = min(score, 0.62)
                downgrade_reason = "mixed_input_direction"
            if related >= 2 and score < 0.18:
                score = 0.18
            result[theme_id] = {
                "theme_id": theme_id,
                "coverage": {
                    "related_input_count": related,
                    "matched_support_count": matched,
                    "ambiguous_support_count": ambiguous,
                    "unmatched_support_count": unmatched,
                    "direction_consistency": direction_consistency,
                    "direction_summary": direction_counts,
                    "biological_pattern_consistency": biological_pattern,
                    "significance_support_count": int(entry["significance_support_count"]),
                    "downgrade_reason": downgrade_reason,
                },
                "input_terms": sorted(entry["input_terms"])[:20],
                "input_ids": sorted(entry["input_ids"])[:20],
                "p_user_sum": round(p_user_sum, 6),
                "calibrated_confidence": round(score, 6),
            }
        return result

    def merge_theme_coverages(self, coverage_by_theme: dict[str, dict[str, Any]], theme_ids: Iterable[str]) -> dict[str, Any]:
        merged = {
            "related_input_count": 0,
            "matched_support_count": 0,
            "ambiguous_support_count": 0,
            "unmatched_support_count": 0,
            "direction_consistency": 0.0,
            "direction_summary": Counter(),
            "significance_support_count": 0,
            "downgrade_reason": "",
        }
        input_ids: set[str] = set()
        reasons: set[str] = set()
        for theme_id in theme_ids:
            row = coverage_by_theme.get(theme_id, {})
            coverage = row.get("coverage") or {}
            input_ids.update(row.get("input_ids") or [])
            merged["matched_support_count"] += int(coverage.get("matched_support_count") or 0)
            merged["ambiguous_support_count"] += int(coverage.get("ambiguous_support_count") or 0)
            merged["unmatched_support_count"] += int(coverage.get("unmatched_support_count") or 0)
            merged["significance_support_count"] += int(coverage.get("significance_support_count") or 0)
            for direction, count in (coverage.get("direction_summary") or {}).items():
                merged["direction_summary"][direction] += int(count or 0)
            if coverage.get("downgrade_reason"):
                reasons.add(str(coverage["downgrade_reason"]))
        support_total = int(merged["matched_support_count"]) + int(merged["ambiguous_support_count"]) + int(merged["unmatched_support_count"])
        merged["related_input_count"] = max(len(input_ids), support_total)
        up = int(merged["direction_summary"].get("up", 0))
        down = int(merged["direction_summary"].get("down", 0))
        directed_total = up + down
        merged["direction_consistency"] = round(max(up, down) / directed_total, 6) if directed_total else 0.0
        merged["direction_summary"] = dict(sorted(merged["direction_summary"].items()))
        merged["downgrade_reason"] = ";".join(sorted(reasons))
        return merged

    def apply_theme_confidence_gates(self, score: float, coverage: dict[str, Any]) -> tuple[float, str, list[str]]:
        capped = float(score or 0.0)
        reasons: list[str] = []
        matched = int(coverage.get("matched_support_count") or 0)
        significant = int(coverage.get("significance_support_count") or 0)
        direction = float(coverage.get("direction_consistency") or 0.0)
        if matched <= 0:
            capped = min(capped, 0.44)
            reasons.append("no_matched_support")
        if significant <= 0:
            capped = min(capped, 0.74)
            reasons.append("no_significant_support")
        if direction == 0.0:
            capped = min(capped, 0.44)
            reasons.append("no_direction_consistency")
        elif direction < 0.5:
            capped = min(capped, 0.74)
            reasons.append("low_direction_consistency")
        return round(capped, 6), self.confidence_tier(capped), reasons

    def biochemical_label_terms(self, normalized_label: str) -> list[str]:
        tokens = [
            token
            for token in re.findall(r"[a-z][a-z0-9]+", normalized_label)
            if len(token) >= 4 and token not in BIOCHEMICAL_TOKEN_STOPWORDS
        ]
        terms = set(tokens)
        for index in range(len(tokens) - 1):
            phrase = f"{tokens[index]} {tokens[index + 1]}"
            if len(phrase) >= 9:
                terms.add(phrase)
        return sorted(terms)

    def pathway_input_theme_boost(
        self,
        pathway: dict[str, Any],
        theme_profile: dict[str, dict[str, Any]],
    ) -> tuple[float, list[dict[str, Any]]]:
        pathway_name = normalize_lookup_key(pathway.get("name") or pathway.get("display_name") or "")
        if not pathway_name or not theme_profile:
            return 0.0, []
        hits = []
        total_boost = 0.0
        for theme_id, theme in sorted(theme_profile.items()):
            support_count = int(theme.get("support_count") or 0)
            if support_count < 2:
                continue
            spec = BIOCHEMICAL_THEME_SPECS.get(theme_id, {})
            matched_terms = [term for term in spec.get("pathway_terms", ()) if normalized_text_has_term(pathway_name, term)]
            if not matched_terms:
                continue
            p_user_sum = float(theme.get("p_user_sum") or 0.0)
            boost = min(4.0, 0.8 + 0.55 * support_count + 0.2 * p_user_sum)
            total_boost = max(total_boost, boost)
            hits.append(
                {
                    "theme_id": theme_id,
                    "support_count": support_count,
                    "pathway_terms": matched_terms,
                    "input_terms": theme.get("input_terms", [])[:10],
                    "boost": round(boost, 6),
                }
            )
        return round(total_boost, 6), hits

    def pathway_lexical_input_boost(
        self,
        pathway: dict[str, Any],
        theme_profile: dict[str, dict[str, Any]],
    ) -> tuple[float, list[dict[str, Any]]]:
        lexical = theme_profile.get("_lexical", {})
        terms = lexical.get("terms", {}) if isinstance(lexical, dict) else {}
        pathway_name = normalize_lookup_key(pathway.get("name") or pathway.get("display_name") or "")
        if not terms or not pathway_name:
            return 0.0, []
        hits = []
        total = 0.0
        for term, meta in terms.items():
            if not normalized_text_has_term(pathway_name, term):
                continue
            support_count = int(meta.get("support_count") or 0)
            if support_count < 2 and " " not in term:
                continue
            p_user_sum = float(meta.get("p_user_sum") or 0.0)
            boost = min(1.25, 0.25 + 0.2 * support_count + 0.1 * p_user_sum)
            total += boost
            hits.append({"term": term, "support_count": support_count, "boost": round(boost, 6)})
        hits.sort(key=lambda row: (-float(row["boost"]), row["term"]))
        return round(min(2.5, total), 6), hits[:8]

    def generic_pathway_penalty(self, pathway: dict[str, Any]) -> float:
        score_components = pathway.get("score_components", {})
        if score_components.get("input_theme_boost") or score_components.get("pathway_theme_boost"):
            return 0.0
        overlap_count = int(pathway.get("overlap_count") or score_components.get("overlap_count") or 0)
        pathway_name = normalize_lookup_key(pathway.get("name") or pathway.get("display_name") or "")
        if overlap_count <= 1 and any(term in pathway_name for term in GENERIC_PATHWAY_TERMS):
            return 0.25
        return 0.0

    def pathway_context_penalty(self, pathway: dict[str, Any]) -> float:
        pathway_name = normalize_lookup_key(pathway.get("name") or pathway.get("display_name") or "")
        disease_model_terms = ("defective", "causes", "syndrome", "deficiency")
        return 0.75 if any(term in pathway_name for term in disease_model_terms) else 0.0

    def disease_or_model_specific_result(self, display: Any, result_type: str = "") -> dict[str, Any]:
        normalized = normalize_lookup_key(display)
        disease_terms = (
            "defective",
            "causes",
            "syndrome",
            "deficiency",
            "disease",
            "disorder",
            "congenital",
            "cancer cells",
            "carcinoma cells",
            "tumor cells",
            "triple negative breast cancer",
            "cell line",
        )
        if result_type == "disease" or any(term in normalized for term in disease_terms):
            return {
                "appendix": True,
                "appendix_reason": "disease_or_model_specific_result_not_primary_metabolic_theme",
                "downgrade_reason": "disease_or_model_specific_result_not_primary_metabolic_theme",
            }
        return {"appendix": False}

    def prediction_context_mismatch(
        self,
        display: Any,
        context: dict[str, Any] | None,
        result_type: str = "",
    ) -> dict[str, Any]:
        normalized_display = normalize_lookup_key(display)
        context_terms = set((context or {}).get("terms", []) or [])
        if not normalized_display:
            return {"context_mismatch": False}
        if not context_terms:
            mismatched = []
            matched_terms = []
            for group, terms in CONTEXT_MISMATCH_GROUPS.items():
                hits = [term for term in terms if term in normalized_display]
                if hits:
                    mismatched.append(group)
                    matched_terms.extend(hits)
            if mismatched:
                return {
                    "context_mismatch": True,
                    "context_mismatch_groups": sorted(set(mismatched)),
                    "context_mismatch_terms": sorted(set(matched_terms)),
                    "context_penalty_reason": "context_specific_result_without_matching_background",
                }
            return {"context_mismatch": False}
        if result_type == "disease":
            display_cancer_groups = {
                group
                for group, terms in CANCER_CONTEXT_GROUPS.items()
                if any(normalized_text_has_term(normalized_display, term) for term in terms)
            }
            context_cancer_groups = {
                group
                for group, terms in CANCER_CONTEXT_GROUPS.items()
                if any(
                    normalized_text_has_term(context_term, term) or normalized_text_has_term(term, context_term)
                    for context_term in context_terms
                    for term in terms
                )
            }
            if display_cancer_groups and context_cancer_groups and not display_cancer_groups.intersection(context_cancer_groups):
                return {
                    "context_mismatch": True,
                    "context_mismatch_groups": sorted(display_cancer_groups),
                    "context_mismatch_terms": [normalized_display],
                    "context_penalty_reason": "disease_context_mismatch_for_current_cancer_background",
                }
        epithelial_context = any(
            normalized_text_has_term(term, epithelial_term)
            for term in context_terms
            for epithelial_term in EPITHELIAL_CONTEXT_TERMS
        )
        brain_context = any(term in context_terms for term in ("brain", "glioma", "glioblastoma", "astrocyte", "astrocytic", "glial", "neuron", "neuronal"))
        immune_context = any(term in context_terms for term in ("immune", "lymphocyte", "t cell", "b cell", "macrophage", "monocyte", "neutrophil", "hematopoietic"))
        taste_context = any(term in context_terms for term in ("taste", "gustatory", "olfactory"))
        mismatched_groups = []
        matched_terms = []
        for group, terms in CONTEXT_MISMATCH_GROUPS.items():
            hits = [term for term in terms if term in normalized_display]
            if not hits:
                continue
            context_has_group = any(any(group_term in context_term for group_term in terms) for context_term in context_terms)
            group_allowed = context_has_group
            if group in {"astrocytic_glial", "neuronal"} and brain_context:
                group_allowed = True
            if group in {"immune_cell_specific", "hematopoietic"} and immune_context:
                group_allowed = True
            if group == "taste_perception" and taste_context:
                group_allowed = True
            if (epithelial_context or not group_allowed) and not group_allowed:
                mismatched_groups.append(group)
                matched_terms.extend(hits)
        if not mismatched_groups:
            return {"context_mismatch": False}
        return {
            "context_mismatch": True,
            "context_mismatch_groups": sorted(set(mismatched_groups)),
            "context_mismatch_terms": sorted(set(matched_terms)),
            "context_penalty_reason": "semantic_context_mismatch_for_current_background",
        }

    def classify_pathway_result_type(self, pathway_name: Any) -> str:
        normalized = normalize_lookup_key(pathway_name)
        if not normalized:
            return "unknown_mixed"
        for result_type, terms in PATHWAY_TYPE_KEYWORDS.items():
            if any(normalized_text_has_term(normalized, term) for term in terms):
                return result_type
        return "unknown_mixed"

    def prediction_task_for_result(self, id_key: str, result_type: str) -> str:
        if id_key == "target_uid":
            return "target_prediction"
        if id_key == "disease_uid":
            return "phenotype_context_prediction"
        if result_type == "transport":
            return "transport_process_prediction"
        if result_type in {"metabolic", "biochemical_process"}:
            return "pathway_prediction"
        return f"{result_type}_prediction"

    def graph_distance_for_ranking(
        self,
        id_key: str,
        source_row: dict[str, Any],
        path: dict[str, Any] | None,
    ) -> int:
        if id_key == "pathway_uid" and int(source_row.get("overlap_count") or 0) > 0:
            return 1
        if path and path.get("edges"):
            return len(path.get("edges") or [])
        if id_key == "pathway_uid" and source_row.get("matched_metabolite_uids"):
            return 1
        return 3 if id_key in {"target_uid", "disease_uid"} else 2

    def pathway_theme_support_count(self, row: dict[str, Any]) -> int:
        support_counts = []
        directional = row.get("directional_support") or {}
        for hit in directional.get("input_theme_hits", []) or []:
            if isinstance(hit, dict):
                support_counts.append(int(hit.get("support_count") or 0))
        for hit in directional.get("input_lexical_hits", []) or []:
            if isinstance(hit, dict):
                support_counts.append(int(hit.get("support_count") or 0))
        return max(support_counts, default=0)

    def path_seed_support_count(self, path: dict[str, Any] | None, features_by_uid: dict[str, list[dict[str, Any]]]) -> int:
        if not path:
            return 0
        seed_uids = set(features_by_uid)
        return len(seed_uids & set(path.get("node_uids") or []))

    def ranking_input_support_count(
        self,
        id_key: str,
        source_row: dict[str, Any],
        path: dict[str, Any] | None,
        features_by_uid: dict[str, list[dict[str, Any]]],
    ) -> int:
        direct_overlap = int(source_row.get("overlap_count") or 0)
        matched_uids = len(set(source_row.get("matched_metabolite_uids") or []))
        theme_support = self.pathway_theme_support_count(source_row) if id_key == "pathway_uid" else 0
        path_support = self.path_seed_support_count(path, features_by_uid)
        return max(direct_overlap, matched_uids, theme_support, path_support)

    def ranking_direction_consistency(self, row: dict[str, Any]) -> float:
        directional = row.get("directional_support") or {}
        up = int(directional.get("up") or 0)
        down = int(directional.get("down") or 0)
        known = up + down
        if known <= 0:
            return 0.0
        return round(max(up, down) / known, 6)

    def ranking_significant_support_count(
        self,
        source_row: dict[str, Any],
        features_by_uid: dict[str, list[dict[str, Any]]],
        path: dict[str, Any] | None = None,
    ) -> int:
        directional = source_row.get("directional_support") or {}
        if "significant" in directional:
            return int(directional.get("significant") or 0)
        support_uids = set(source_row.get("matched_metabolite_uids") or [])
        if path:
            support_uids.update(set(path.get("node_uids") or []) & set(features_by_uid))
        return sum(
            1
            for metabolite_uid in support_uids
            for feature in features_by_uid.get(metabolite_uid, [])
            if feature.get("significant")
        )

    def node_genericity_weight(self, node_uid: str) -> tuple[float, int]:
        node = self.get_node(node_uid)
        degree = 0
        if node and node.get("node_idx") is not None and (self.graph_dir / "sparse_edges.parquet").exists():
            degree = len(self.sparse_incident_edges(int(node["node_idx"])))
        elif node_uid:
            degree = len(self.incident_edges(node_uid))
        if degree <= 20:
            return 1.0, degree
        if degree <= 100:
            return 0.75, degree
        if degree <= 500:
            return 0.5, degree
        return 0.25, degree

    def distance_weight(self, graph_distance: int) -> float:
        return PREDICTION_DISTANCE_WEIGHTS.get(graph_distance, 0.25 if graph_distance >= 3 else 0.55)

    def match_confidence_weight(self, input_summary: dict[str, Any]) -> float:
        input_count = int(input_summary.get("input_count") or 0)
        matched = int(input_summary.get("matched_count") or 0)
        ambiguous = int(input_summary.get("ambiguous_count") or 0)
        unmatched = int(input_summary.get("unmatched_count") or 0) + int(input_summary.get("invalid_count") or 0)
        expanded_inputs = int(input_summary.get("expanded_candidate_input_count") or input_summary.get("expanded_input_count") or 0)
        analysis_seed_inputs = int(input_summary.get("analysis_seed_input_count") or matched + expanded_inputs)
        if input_count <= 0:
            return 0.0
        matched_ratio = matched / input_count
        ambiguous_ratio = ambiguous / input_count
        unmatched_ratio = unmatched / input_count
        seed_ratio = analysis_seed_inputs / input_count
        effective_ratio = (matched + 0.5 * expanded_inputs) / input_count
        if matched_ratio >= 0.8 and ambiguous_ratio <= 0.1 and unmatched_ratio <= 0.1:
            return 1.0
        if matched_ratio >= 0.6 and ambiguous_ratio <= 0.25:
            return 0.8
        if matched_ratio >= 0.4:
            return 0.55
        if seed_ratio >= 0.6 and effective_ratio >= 0.5 and unmatched_ratio <= 0.35:
            return 0.55
        if seed_ratio >= 0.45 and effective_ratio >= 0.35:
            return 0.4
        return 0.25

    def confidence_tier(self, score: float) -> str:
        for tier, threshold in CONFIDENCE_TIERS:
            if score >= threshold:
                return tier
        return "low"

    def bounded_confidence(self, score: Any, precision: int = 12) -> float:
        try:
            value = float(score or 0.0)
        except (TypeError, ValueError):
            value = 0.0
        return round(min(1.0, max(0.0, value)), precision)

    def calibration_status(
        self,
        confidence_tier: str,
        graph_distance: int,
        result_type: str,
        input_support_count: int,
        evidence_sources: list[str],
    ) -> str:
        if input_support_count <= 0:
            return "limited_input_support"
        if result_type in {"drug", "cell_context_overlay"} or "literature_overlay" in evidence_sources:
            return "overlay_informed"
        if graph_distance >= 3:
            return "graph_extrapolation"
        if confidence_tier in {"high", "medium"}:
            return "calibrated_in_scope"
        return "exploratory_in_scope"

    def prediction_boundary_note(
        self,
        confidence_tier: str,
        graph_distance: int,
        result_type: str,
        calibration_status: str,
    ) -> str:
        if confidence_tier == "high":
            return "High-confidence model prediction; still report as a prediction with evidence and context."
        if confidence_tier == "medium":
            return "Medium-confidence prediction; useful as a candidate mechanism and should be reviewed with context."
        if calibration_status in {"graph_extrapolation", "overlay_informed"} or graph_distance >= 2:
            return "Exploratory prediction from graph, literature, or overlay evidence; validate before use."
        return "Low-confidence prediction; treat as a weak lead or deprioritize."

    def build_prediction_text(self, display: str, prediction_task: str, confidence_tier: str) -> str:
        label = str(display or "result")
        if prediction_task == "target_prediction":
            return f"The model predicts {label} as a candidate explanatory target."
        if prediction_task == "phenotype_context_prediction":
            return f"The model predicts {label} as a possible phenotype or disease-context signal."
        return f"The model predicts {label} as a {confidence_tier}-confidence metabolic result."

    def calibrated_prediction_for_ranking_row(
        self,
        packed: dict[str, Any],
        source_row: dict[str, Any],
        id_key: str,
        path: dict[str, Any] | None,
        features_by_uid: dict[str, list[dict[str, Any]]],
        input_summary: dict[str, Any],
        context: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        node_uid = str(packed.get(id_key) or "")
        display = packed.get("display_name") or source_row.get("name") or source_row.get("display_name") or node_uid
        if id_key == "pathway_uid":
            result_type = self.classify_pathway_result_type(display)
            if result_type == "unknown_mixed" and any(
                source_row.get("score_components", {}).get(key)
                for key in ("input_theme_boost", "pathway_theme_boost", "input_lexical_boost")
            ):
                result_type = "biochemical_process"
        elif id_key == "target_uid":
            result_type = "target"
        elif id_key == "disease_uid":
            result_type = "disease"
        else:
            result_type = "unknown_mixed"
        prediction_task = self.prediction_task_for_result(id_key, result_type)
        graph_distance = self.graph_distance_for_ranking(id_key, source_row, path)
        input_support_count = self.ranking_input_support_count(id_key, source_row, path, features_by_uid)
        direct_overlap = int(source_row.get("overlap_count") or 0)
        is_directly_supported = id_key == "pathway_uid" and direct_overlap > 0 and graph_distance <= 1
        has_trace = bool(packed.get("claim_refs", {}).get("traceability_passed"))
        literature_count = int(packed.get("literature_support_count") or 0)
        support_classes = set(packed.get("support_classes") or [])
        score_components = packed.get("score_components") or {}
        propagation_score = float(score_components.get("propagation_score") or 0.0)
        evidence_sources = []
        if is_directly_supported:
            evidence_sources.extend(["input", "database_pathway"])
        elif input_support_count > 0 and id_key == "pathway_uid":
            evidence_sources.append("input_theme")
        if propagation_score > 0.0 or graph_distance >= 2:
            evidence_sources.append("graph_propagation")
        if literature_count:
            evidence_sources.append("literature_overlay" if "novel_candidate" in support_classes else "literature")
        if not evidence_sources and has_trace:
            evidence_sources.append("database")
        type_weight = PREDICTION_TYPE_WEIGHTS.get(result_type, 0.5)
        direction_consistency = self.ranking_direction_consistency(source_row)
        significant_support_count = self.ranking_significant_support_count(source_row, features_by_uid, path)
        genericity_weight, node_degree = self.node_genericity_weight(node_uid)
        input_support_weight = min(1.0, 0.2 + 0.2 * max(0, input_support_count))
        confidence_weight = self.match_confidence_weight(input_summary)
        direction_weight = 0.75 + 0.25 * direction_consistency
        distance_penalty = self.distance_weight(graph_distance)
        raw_score = float(packed.get("score") or 0.0)
        calibrated_confidence = self.bounded_confidence(
            raw_score
            * input_support_weight
            * confidence_weight
            * direction_weight
            * type_weight
            * distance_penalty
            * genericity_weight,
        )
        context_eval = self.prediction_context_mismatch(display, context, result_type)
        disease_eval = self.disease_or_model_specific_result(display, result_type)
        pre_context_confidence = calibrated_confidence
        gate_reasons: list[str] = []
        if input_support_count <= 0:
            calibrated_confidence = min(calibrated_confidence, 0.44)
            gate_reasons.append("no_matched_support")
        if significant_support_count <= 0:
            calibrated_confidence = min(calibrated_confidence, 0.74)
            gate_reasons.append("no_significant_support")
        if direction_consistency == 0.0:
            calibrated_confidence = min(calibrated_confidence, 0.44)
            gate_reasons.append("no_direction_consistency")
        elif direction_consistency < 0.5:
            calibrated_confidence = min(calibrated_confidence, 0.74)
            gate_reasons.append("low_direction_consistency")
        if (
            id_key == "pathway_uid"
            and input_summary.get("analysis_mode") == "two_group_trait_comparison_table"
            and (direct_overlap < 2 or significant_support_count < 2 or literature_count <= 0)
        ):
            calibrated_confidence = min(calibrated_confidence, 0.74)
            gate_reasons.append("trait_score_requires_multi_seed_or_literature_support")
        context_mismatch = bool(context_eval.get("context_mismatch"))
        if context_mismatch:
            calibrated_confidence = min(calibrated_confidence, 0.25)
        disease_model_specific = bool(disease_eval.get("appendix"))
        if disease_model_specific:
            calibrated_confidence = min(calibrated_confidence, 0.25)
        tier = self.confidence_tier(calibrated_confidence)
        if context_mismatch:
            status = "context_mismatch"
        elif disease_model_specific:
            status = "appendix_low"
        else:
            status = self.calibration_status(tier, graph_distance, result_type, input_support_count, sorted(set(evidence_sources)))
        downgrade_reasons = [
            *gate_reasons,
            context_eval.get("context_penalty_reason", "") if context_mismatch else "",
            disease_eval.get("downgrade_reason", "") if disease_model_specific else "",
        ]
        appendix = disease_model_specific or (result_type in {"target", "disease"} and (tier == "low" or confidence_weight < 0.6))
        appendix_reason = (
            disease_eval.get("appendix_reason", "")
            if disease_model_specific
            else ("low_or_input_limited_upstream_confidence" if appendix else "")
        )
        return {
            "prediction_id": f"{node_uid}:rank_{packed.get('rank', '')}",
            "prediction": self.build_prediction_text(str(display), prediction_task, tier),
            "prediction_task": prediction_task,
            "result_type": result_type,
            "confidence_tier": tier,
            "calibrated_confidence": calibrated_confidence,
            "calibration_status": status,
            "input_support_count": input_support_count,
            "matched_support_count": input_support_count,
            "significant_support_count": significant_support_count,
            "matched_input_count": int(input_summary.get("matched_count") or 0),
            "ambiguous_input_count": int(input_summary.get("ambiguous_count") or 0),
            "evidence_sources": sorted(set(evidence_sources)),
            "graph_distance": graph_distance,
            "direction_consistency": direction_consistency,
            **context_eval,
            "is_directly_supported": is_directly_supported,
            "is_extrapolated": graph_distance >= 2 or result_type in {"target", "disease", "drug", "cell_context_overlay"},
            "needs_validation": tier != "high" or graph_distance >= 2,
            "boundary": self.prediction_boundary_note(tier, graph_distance, result_type, status),
            "downgrade_reason": ";".join(reason for reason in downgrade_reasons if reason),
            "ranker_score": calibrated_confidence,
            "calibration_components": {
                "raw_score": raw_score,
                "pre_context_confidence": round(pre_context_confidence, 6),
                "input_support_weight": round(input_support_weight, 6),
                "match_confidence_weight": round(confidence_weight, 6),
                "direction_consistency_weight": round(direction_weight, 6),
                "significant_support_count": significant_support_count,
                "confidence_gate_reasons": gate_reasons,
                "result_type_weight": round(type_weight, 6),
                "distance_weight": round(distance_penalty, 6),
                "genericity_weight": round(genericity_weight, 6),
                "context_penalty_applied": context_mismatch,
                "node_degree": node_degree,
                "formula": "raw_score*input_support_weight*match_confidence_weight*direction_consistency_weight*result_type_weight*distance_weight*genericity_weight, bounded to [0,1], then context mismatch cap if applicable",
            },
            "appendix": appendix,
            "appendix_reason": appendix_reason,
        }

    def pathway_theme_boost(
        self,
        pathway: dict[str, Any],
        features_by_uid: dict[str, list[dict[str, Any]]],
    ) -> tuple[float, list[str]]:
        pathway_name = normalize_lookup_key(pathway.get("name") or pathway.get("display_name") or "")
        if not pathway_name:
            return 0.0, []
        hits: list[str] = []
        for metabolite_uid in pathway.get("matched_metabolite_uids", []):
            for feature in features_by_uid.get(metabolite_uid, []):
                for label in (
                    normalize_lookup_key(feature.get("input_name", "")),
                    normalize_lookup_key(feature.get("display_name", "")),
                ):
                    if not label:
                        continue
                    terms = PATHWAY_THEME_TERMS.get(label, [])
                    if terms and any(term in pathway_name for term in terms):
                        hits.append(label)
                        break
        hits = sorted(set(hits))
        return round(min(1.0, 0.85 * len(hits)), 6), hits

    def directional_pathway_enrichment(self, features_by_uid: dict[str, list[dict[str, Any]]]) -> dict[str, Any]:
        grouped = self.directional_uids(features_by_uid)
        seed_presence_weights = self.seed_presence_weights_from_features(features_by_uid)
        result = {
            "summary": {key: len(value) for key, value in grouped.items()},
            "up": self.attach_pathway_user_support(self.pathway_enrichment(grouped["up"], seed_presence_weights)[:50], features_by_uid),
            "down": self.attach_pathway_user_support(self.pathway_enrichment(grouped["down"], seed_presence_weights)[:50], features_by_uid),
            "significant_up": self.attach_pathway_user_support(
                self.pathway_enrichment(grouped["significant_up"], seed_presence_weights)[:50], features_by_uid
            ),
            "significant_down": self.attach_pathway_user_support(
                self.pathway_enrichment(grouped["significant_down"], seed_presence_weights)[:50], features_by_uid
            ),
        }
        return result

    def collect_propagation_graph(self, seed_weights: dict[str, float], max_hops: int | None = None) -> dict[str, Any]:
        max_hops = max_hops or self.config.max_hops
        sparse_path = self.graph_dir / "sparse_edges.parquet"
        if sparse_path.exists() and (self.graph_dir / "edge_type_index.parquet").exists():
            return self.collect_sparse_propagation_graph(seed_weights, max_hops)
        return self.collect_edge_propagation_graph(seed_weights, max_hops)

    def empty_compression_stats(self) -> dict[str, Any]:
        return {
            "enabled": True,
            "mode": "typed_beam_mass_retention",
            "retained_mass_threshold": self.config.propagation_retained_mass,
            "beam_per_type": self.config.propagation_beam_per_type,
            "max_edges_per_node": self.config.max_edges_per_node,
            "max_adaptive_edges_per_node": self.config.propagation_max_adaptive_edges_per_node,
            "degree_penalty": self.config.propagation_degree_penalty,
            "candidate_edge_count": 0,
            "retained_edge_count": 0,
            "omitted_edge_count": 0,
            "budget_limited_edge_count": 0,
            "candidate_node_count": 0,
            "retained_node_count": 0,
            "omitted_node_count": 0,
            "budget_limited_node_count": 0,
            "transition_retained_mass": 1.0,
            "frontier_retained_mass": 1.0,
            "retained_mass": 1.0,
            "compressed": False,
            "omitted_mass": 0.0,
            "hops": [],
        }

    def effective_propagation_weight(self, weight: Any, degree: int) -> float:
        base = clamp_probability(weight, 0.5)
        if not self.config.propagation_degree_penalty:
            return base
        return base / max(1.0, math.log2(2.0 + max(0, degree)))

    def propagation_edge_limit(self) -> int:
        return max(1, max(int(self.config.max_edges_per_node), int(self.config.propagation_max_adaptive_edges_per_node)))

    def finalize_compression_stats(self, stats: dict[str, Any]) -> dict[str, Any]:
        transition_total = float(stats.pop("_transition_total", 0.0))
        transition_retained = float(stats.pop("_transition_retained", 0.0))
        frontier_total = float(stats.pop("_frontier_total", 0.0))
        frontier_retained = float(stats.pop("_frontier_retained", 0.0))
        transition_ratio = 1.0 if transition_total <= 0 else transition_retained / transition_total
        frontier_ratio = 1.0 if frontier_total <= 0 else frontier_retained / frontier_total
        retained_mass = min(1.0, max(0.0, min(transition_ratio, frontier_ratio)))
        stats["transition_retained_mass"] = round(transition_ratio, 6)
        stats["frontier_retained_mass"] = round(frontier_ratio, 6)
        stats["retained_mass"] = round(retained_mass, 6)
        stats["omitted_mass"] = round(max(0.0, 1.0 - retained_mass), 6)
        stats["compressed"] = bool(stats["omitted_edge_count"] or stats["omitted_node_count"])
        stats["retained_mass_meets_threshold"] = retained_mass >= float(self.config.propagation_retained_mass)
        stats["quality_status"] = "passed" if stats["retained_mass_meets_threshold"] else "lossy_below_threshold"
        return stats

    def collect_sparse_propagation_graph(self, seed_weights: dict[str, float], max_hops: int) -> dict[str, Any]:
        seed_nodes = {uid: self.get_node(uid) for uid in sorted(seed_weights)}
        seed_idx_to_uid = {
            int(node["node_idx"]): uid
            for uid, node in seed_nodes.items()
            if node and node.get("node_idx") is not None
        }
        total_seed_weight = sum(float(weight) for weight in seed_weights.values()) or 1.0
        frontier_scores = {
            node_idx: float(seed_weights[uid]) / total_seed_weight
            for node_idx, uid in seed_idx_to_uid.items()
        }
        frontier = set(frontier_scores)
        seen_idxs = set(seed_idx_to_uid)
        seen_edges: dict[str, dict[str, Any]] = {}
        truncated = {"nodes": False, "edges": False}
        compression = self.empty_compression_stats()
        compression["_transition_total"] = 0.0
        compression["_transition_retained"] = 0.0
        compression["_frontier_total"] = 0.0
        compression["_frontier_retained"] = 0.0
        for hop in range(max_hops):
            if not frontier:
                break
            next_scores: dict[int, float] = defaultdict(float)
            incident_by_idx = self.sparse_incident_edges_many(list(frontier))
            hop_stats = {
                "hop": hop + 1,
                "frontier_count": len(frontier),
                "candidate_edge_count": 0,
                "retained_edge_count": 0,
                "omitted_edge_count": 0,
                "candidate_node_count": 0,
                "retained_node_count": 0,
                "omitted_node_count": 0,
                "transition_retained_mass": 1.0,
                "frontier_retained_mass": 1.0,
            }
            for node_idx in sorted(frontier):
                candidate_edges = []
                source_score = float(frontier_scores.get(node_idx, 0.0))
                if source_score <= 0.0:
                    continue
                incident = incident_by_idx.get(node_idx, [])
                degree = len(incident)
                for edge in incident:
                    current_type = edge["subject_type"] if edge["source_idx"] == node_idx else edge["object_type"]
                    neighbor_type = edge["object_type"] if edge["source_idx"] == node_idx else edge["subject_type"]
                    if not propagation_allows(current_type, edge["edge_type"], neighbor_type):
                        continue
                    effective_weight = self.effective_propagation_weight(edge["weight"], degree)
                    candidate_edges.append(
                        (
                            propagation_priority(current_type, edge["edge_type"], neighbor_type),
                            -effective_weight,
                            edge["edge_uid"],
                            neighbor_type,
                            {**edge, "effective_weight": effective_weight, "degree_penalty_applied": self.config.propagation_degree_penalty},
                        )
                    )
                if not candidate_edges:
                    continue
                candidate_edges = sorted(candidate_edges)
                total_weight = sum(-item[1] for item in candidate_edges) or 1.0
                retained = []
                retained_weight = 0.0
                adaptive_edge_limit = self.propagation_edge_limit()
                for item in candidate_edges:
                    if len(retained) >= adaptive_edge_limit:
                        break
                    retained.append(item)
                    retained_weight += -item[1]
                    if (
                        len(retained) >= self.config.max_edges_per_node
                        and retained_weight / total_weight >= self.config.propagation_retained_mass
                    ):
                        break
                omitted_count = max(0, len(candidate_edges) - len(retained))
                omitted_weight = max(0.0, total_weight - retained_weight)
                compression["candidate_edge_count"] += len(candidate_edges)
                compression["retained_edge_count"] += len(retained)
                compression["omitted_edge_count"] += omitted_count
                compression["_transition_total"] += source_score
                compression["_transition_retained"] += source_score * retained_weight / total_weight
                hop_stats["candidate_edge_count"] += len(candidate_edges)
                hop_stats["retained_edge_count"] += len(retained)
                hop_stats["omitted_edge_count"] += omitted_count
                terminal_node_limit = self.config.max_propagation_nodes + adaptive_edge_limit
                for _priority, neg_weight, _edge_uid, neighbor_type, edge in retained:
                    if edge["edge_uid"] not in seen_edges:
                        if len(seen_edges) >= self.config.max_propagation_edges:
                            compression["omitted_edge_count"] += 1
                            compression["budget_limited_edge_count"] += 1
                            hop_stats["omitted_edge_count"] += 1
                            continue
                        seen_edges[edge["edge_uid"]] = edge
                    neighbor_idx = edge["target_idx"] if edge["source_idx"] == node_idx else edge["source_idx"]
                    next_scores[neighbor_idx] += source_score * (-neg_weight) / total_weight
                if omitted_count or omitted_weight > 0:
                    hop_stats["transition_retained_mass"] = min(
                        hop_stats["transition_retained_mass"],
                        round(retained_weight / total_weight, 6),
                    )

            compression["candidate_node_count"] += len(next_scores)
            hop_stats["candidate_node_count"] = len(next_scores)
            by_type: dict[str, list[tuple[float, int]]] = defaultdict(list)
            node_rows_for_next = self.get_nodes_by_idx(list(next_scores))
            idx_to_type = {int(row["node_idx"]): row.get("node_type", "") for row in node_rows_for_next}
            for idx, score in next_scores.items():
                by_type[idx_to_type.get(idx, "")].append((float(score), int(idx)))

            next_frontier_scores: dict[int, float] = {}
            total_next_mass = sum(next_scores.values())
            retained_next_mass = 0.0
            for node_type in sorted(by_type):
                rows = sorted(by_type[node_type], key=lambda item: (-item[0], item[1]))
                retained_rows = []
                type_total = sum(score for score, _idx in rows) or 1.0
                type_retained = 0.0
                for score, idx in rows:
                    if (
                        len(retained_rows) >= max(1, self.config.propagation_beam_per_type)
                        and type_retained / type_total >= self.config.propagation_retained_mass
                    ):
                        break
                    retained_rows.append((score, idx))
                    type_retained += score
                for score, idx in retained_rows:
                    if idx not in seen_idxs:
                        if len(seen_idxs) >= self.config.max_propagation_nodes and node_type != "disease":
                            compression["omitted_node_count"] += 1
                            compression["budget_limited_node_count"] += 1
                            hop_stats["omitted_node_count"] += 1
                            continue
                        if len(seen_idxs) >= terminal_node_limit:
                            compression["omitted_node_count"] += 1
                            compression["budget_limited_node_count"] += 1
                            hop_stats["omitted_node_count"] += 1
                            continue
                        seen_idxs.add(idx)
                    next_frontier_scores[idx] = score
                    retained_next_mass += score
                omitted_rows = rows[len(retained_rows):]
                compression["omitted_node_count"] += len(omitted_rows)
                hop_stats["omitted_node_count"] += len(omitted_rows)

            retained_count = len(next_frontier_scores)
            compression["retained_node_count"] += retained_count
            compression["_frontier_total"] += total_next_mass
            compression["_frontier_retained"] += retained_next_mass
            hop_stats["retained_node_count"] = retained_count
            hop_stats["frontier_retained_mass"] = round(1.0 if total_next_mass <= 0 else retained_next_mass / total_next_mass, 6)
            compression["hops"].append(hop_stats)
            frontier_scores = next_frontier_scores
            frontier = set(frontier_scores)

        node_rows = self.get_nodes_by_idx(list(seen_idxs))
        idx_to_node = {int(row["node_idx"]): row for row in node_rows}
        adjacency: dict[str, list[tuple[str, float]]] = defaultdict(list)
        for edge in seen_edges.values():
            source = idx_to_node.get(edge["source_idx"])
            target = idx_to_node.get(edge["target_idx"])
            if not source or not target:
                continue
            source_uid = source["node_uid"]
            target_uid = target["node_uid"]
            weight = float(edge.get("effective_weight") or clamp_probability(edge["weight"], 0.5))
            adjacency[source_uid].append((target_uid, weight))
            adjacency[target_uid].append((source_uid, weight))
        compression = self.finalize_compression_stats(compression)
        return {
            "mode": "sparse_edges",
            "adjacency": adjacency,
            "node_uids": {row["node_uid"] for row in node_rows},
            "edge_count": len(seen_edges),
            "truncated": truncated,
            "compression": compression,
        }

    def collect_edge_propagation_graph(self, seed_weights: dict[str, float], max_hops: int) -> dict[str, Any]:
        total_seed_weight = sum(float(weight) for weight in seed_weights.values()) or 1.0
        frontier_scores = {uid: float(weight) / total_seed_weight for uid, weight in seed_weights.items()}
        frontier = set(frontier_scores)
        seen_uids = set(seed_weights)
        seen_edges: dict[str, dict[str, Any]] = {}
        adjacency: dict[str, list[tuple[str, float]]] = defaultdict(list)
        truncated = {"nodes": False, "edges": False}
        compression = self.empty_compression_stats()
        compression["_transition_total"] = 0.0
        compression["_transition_retained"] = 0.0
        compression["_frontier_total"] = 0.0
        compression["_frontier_retained"] = 0.0
        for hop in range(max_hops):
            if not frontier:
                break
            next_scores: dict[str, float] = defaultdict(float)
            next_types: dict[str, str] = {}
            hop_stats = {
                "hop": hop + 1,
                "frontier_count": len(frontier),
                "candidate_edge_count": 0,
                "retained_edge_count": 0,
                "omitted_edge_count": 0,
                "candidate_node_count": 0,
                "retained_node_count": 0,
                "omitted_node_count": 0,
                "transition_retained_mass": 1.0,
                "frontier_retained_mass": 1.0,
            }
            for node_uid in sorted(frontier):
                current_node = self.get_node(node_uid) or {}
                current_type = str(current_node.get("node_type") or "")
                candidate_edges = []
                source_score = float(frontier_scores.get(node_uid, 0.0))
                if source_score <= 0.0:
                    continue
                incident = self.incident_edges(node_uid)
                degree = len(incident)
                for edge in incident:
                    if edge.get("edge_type") not in ALLOWED_EXPLANATION_EDGE_TYPES:
                        continue
                    neighbor = edge["object_uid"] if edge["subject_uid"] == node_uid else edge["subject_uid"]
                    neighbor_type = edge["object_type"] if edge["subject_uid"] == node_uid else edge["subject_type"]
                    if not propagation_allows(current_type, edge["edge_type"], neighbor_type):
                        continue
                    effective_weight = self.effective_propagation_weight(edge.get("p_final"), degree)
                    candidate_edges.append(
                        (
                            propagation_priority(current_type, edge["edge_type"], neighbor_type),
                            -effective_weight,
                            edge["edge_uid"],
                            neighbor,
                            neighbor_type,
                            {**edge, "effective_weight": effective_weight, "degree_penalty_applied": self.config.propagation_degree_penalty},
                        )
                    )
                if not candidate_edges:
                    continue
                candidate_edges = sorted(candidate_edges)
                total_weight = sum(-item[1] for item in candidate_edges) or 1.0
                retained = []
                retained_weight = 0.0
                adaptive_edge_limit = self.propagation_edge_limit()
                for item in candidate_edges:
                    if len(retained) >= adaptive_edge_limit:
                        break
                    retained.append(item)
                    retained_weight += -item[1]
                    if (
                        len(retained) >= self.config.max_edges_per_node
                        and retained_weight / total_weight >= self.config.propagation_retained_mass
                    ):
                        break
                omitted_count = max(0, len(candidate_edges) - len(retained))
                compression["candidate_edge_count"] += len(candidate_edges)
                compression["retained_edge_count"] += len(retained)
                compression["omitted_edge_count"] += omitted_count
                compression["_transition_total"] += source_score
                compression["_transition_retained"] += source_score * retained_weight / total_weight
                hop_stats["candidate_edge_count"] += len(candidate_edges)
                hop_stats["retained_edge_count"] += len(retained)
                hop_stats["omitted_edge_count"] += omitted_count
                terminal_node_limit = self.config.max_propagation_nodes + adaptive_edge_limit
                for _priority, neg_weight, _edge_uid, neighbor, neighbor_type, edge in retained:
                    is_new_edge = edge["edge_uid"] not in seen_edges
                    if is_new_edge:
                        if len(seen_edges) >= self.config.max_propagation_edges:
                            compression["omitted_edge_count"] += 1
                            compression["budget_limited_edge_count"] += 1
                            hop_stats["omitted_edge_count"] += 1
                            continue
                        seen_edges[edge["edge_uid"]] = edge
                    if is_new_edge:
                        weight = float(edge.get("effective_weight") or clamp_probability(edge.get("p_final"), 0.5))
                        adjacency[node_uid].append((neighbor, weight))
                        adjacency[neighbor].append((node_uid, weight))
                    next_scores[neighbor] += source_score * (-neg_weight) / total_weight
                    next_types[neighbor] = neighbor_type
                if omitted_count:
                    hop_stats["transition_retained_mass"] = min(
                        hop_stats["transition_retained_mass"],
                        round(retained_weight / total_weight, 6),
                    )

            compression["candidate_node_count"] += len(next_scores)
            hop_stats["candidate_node_count"] = len(next_scores)
            by_type: dict[str, list[tuple[float, str]]] = defaultdict(list)
            for uid, score in next_scores.items():
                by_type[next_types.get(uid, "")].append((float(score), uid))

            next_frontier_scores: dict[str, float] = {}
            total_next_mass = sum(next_scores.values())
            retained_next_mass = 0.0
            for node_type in sorted(by_type):
                rows = sorted(by_type[node_type], key=lambda item: (-item[0], item[1]))
                retained_rows = []
                type_total = sum(score for score, _uid in rows) or 1.0
                type_retained = 0.0
                for score, uid in rows:
                    if (
                        len(retained_rows) >= max(1, self.config.propagation_beam_per_type)
                        and type_retained / type_total >= self.config.propagation_retained_mass
                    ):
                        break
                    retained_rows.append((score, uid))
                    type_retained += score
                for score, uid in retained_rows:
                    if uid not in seen_uids:
                        if len(seen_uids) >= self.config.max_propagation_nodes and node_type != "disease":
                            compression["omitted_node_count"] += 1
                            compression["budget_limited_node_count"] += 1
                            hop_stats["omitted_node_count"] += 1
                            continue
                        if len(seen_uids) >= terminal_node_limit:
                            compression["omitted_node_count"] += 1
                            compression["budget_limited_node_count"] += 1
                            hop_stats["omitted_node_count"] += 1
                            continue
                        seen_uids.add(uid)
                    next_frontier_scores[uid] = score
                    retained_next_mass += score
                omitted_rows = rows[len(retained_rows):]
                compression["omitted_node_count"] += len(omitted_rows)
                hop_stats["omitted_node_count"] += len(omitted_rows)

            retained_count = len(next_frontier_scores)
            compression["retained_node_count"] += retained_count
            compression["_frontier_total"] += total_next_mass
            compression["_frontier_retained"] += retained_next_mass
            hop_stats["retained_node_count"] = retained_count
            hop_stats["frontier_retained_mass"] = round(1.0 if total_next_mass <= 0 else retained_next_mass / total_next_mass, 6)
            compression["hops"].append(hop_stats)
            frontier_scores = next_frontier_scores
            frontier = set(frontier_scores)
        compression = self.finalize_compression_stats(compression)
        return {
            "mode": "edges",
            "adjacency": adjacency,
            "node_uids": seen_uids,
            "edge_count": len(seen_edges),
            "truncated": truncated,
            "compression": compression,
        }

    def propagate_scores(self, seed_weights: dict[str, float], max_hops: int | None = None) -> dict[str, Any]:
        seed_weights = {uid: float(weight) for uid, weight in seed_weights.items() if weight and weight > 0}
        requested_hops = max_hops or self.config.max_hops
        if not seed_weights:
            return {
                "formula": "pi=(1-alpha)y+alpha*W*pi",
                "alpha": self.config.propagation_alpha,
                "iterations": 0,
                "converged": True,
                "graph": {
                    "mode": "empty",
                    "node_count": 0,
                    "edge_count": 0,
                    "truncated": {"nodes": False, "edges": False},
                    "compression": self.empty_compression_stats(),
                    "compressed_not_truncated": False,
                },
                "nodes": [],
            }
        graph = self.collect_propagation_graph(seed_weights, max_hops=requested_hops)
        node_uids = set(graph["node_uids"]) | set(seed_weights)
        total_seed_weight = sum(seed_weights.values()) or 1.0
        y = {uid: seed_weights.get(uid, 0.0) / total_seed_weight for uid in node_uids}
        scores = dict(y)
        alpha = clamp_unit(self.config.propagation_alpha, DEFAULT_PROPAGATION_ALPHA)
        adjacency = graph["adjacency"]
        converged = False
        iterations = 0
        for iteration in range(1, self.config.propagation_iterations + 1):
            next_scores = {uid: (1.0 - alpha) * y.get(uid, 0.0) for uid in node_uids}
            for source in sorted(node_uids):
                neighbors = adjacency.get(source, [])
                if not neighbors:
                    continue
                source_score = scores.get(source, 0.0)
                if source_score == 0.0:
                    continue
                denominator = sum(weight for _neighbor, weight in neighbors) or 1.0
                for neighbor, weight in neighbors:
                    next_scores[neighbor] = next_scores.get(neighbor, 0.0) + alpha * source_score * weight / denominator
            delta = max((abs(next_scores.get(uid, 0.0) - scores.get(uid, 0.0)) for uid in node_uids), default=0.0)
            scores = next_scores
            iterations = iteration
            if delta <= self.config.propagation_tolerance:
                converged = True
                break
        if not converged and iterations >= self.config.propagation_iterations:
            converged = False

        terminal_disease_scores: dict[str, float] = {}
        if requested_hops >= 4:
            terminal_disease_scores = self.terminal_disease_expansion(scores, node_uids, alpha)
            for disease_uid, score in terminal_disease_scores.items():
                scores[disease_uid] = scores.get(disease_uid, 0.0) + score
                y.setdefault(disease_uid, 0.0)
                node_uids.add(disease_uid)

        nodes = []
        for node in self.get_nodes(sorted(node_uids)):
            node_uid = node["node_uid"]
            score = scores.get(node_uid, 0.0)
            if score <= 0.0:
                continue
            nodes.append(
                {
                    "node_uid": node_uid,
                    "node_type": node.get("node_type", ""),
                    "display_name": node.get("display_name", ""),
                    "primary_external_id": node.get("primary_external_id", ""),
                    "propagation_score": round(score, 12),
                    "seed_score": round(y.get(node_uid, 0.0), 12),
                    "seed_weight": round(seed_weights.get(node_uid, 0.0), 6),
                }
            )
        nodes.sort(key=lambda row: (-row["propagation_score"], row["node_type"], row["node_uid"]))
        return {
            "formula": "pi=(1-alpha)y+alpha*W*pi",
            "alpha": alpha,
            "iterations": iterations,
            "converged": converged,
            "graph": {
                "mode": f"compressed_{graph['mode']}" if graph.get("compression", {}).get("compressed") else graph["mode"],
                "node_count": len(node_uids),
                "edge_count": graph["edge_count"],
                "truncated": graph["truncated"],
                "compression": graph.get("compression", {}),
                "compressed_not_truncated": bool(
                    graph.get("compression", {}).get("compressed")
                    and not any(graph.get("truncated", {}).values())
                ),
                "compression_quality_warning": bool(
                    graph.get("compression", {}).get("compressed")
                    and not graph.get("compression", {}).get("retained_mass_meets_threshold", True)
                ),
                "terminal_disease_expansion_count": len(terminal_disease_scores),
            },
            "nodes": nodes,
        }

    def terminal_disease_expansion(self, scores: dict[str, float], node_uids: set[str], alpha: float) -> dict[str, float]:
        target_nodes = [
            node
            for node in self.get_nodes(sorted(node_uids))
            if node.get("node_type") == "target" and scores.get(node["node_uid"], 0.0) > 0.0
        ]
        target_nodes.sort(key=lambda row: (-scores.get(row["node_uid"], 0.0), row["node_uid"]))
        target_nodes = target_nodes[: self.config.max_candidates * 10]
        disease_scores: dict[str, float] = defaultdict(float)
        sparse_path = self.graph_dir / "sparse_edges.parquet"
        if sparse_path.exists() and (self.graph_dir / "edge_type_index.parquet").exists():
            disease_idx_scores: dict[int, float] = defaultdict(float)
            incident_by_idx = self.sparse_incident_edges_many([int(target["node_idx"]) for target in target_nodes])
            for target in target_nodes:
                target_idx = int(target["node_idx"])
                edges = []
                for edge in incident_by_idx.get(target_idx, []):
                    current_type = edge["subject_type"] if edge["source_idx"] == target_idx else edge["object_type"]
                    neighbor_type = edge["object_type"] if edge["source_idx"] == target_idx else edge["subject_type"]
                    if (current_type, edge["edge_type"], neighbor_type) != ("target", "target_associated_with_disease", "disease"):
                        continue
                    edges.append((-float(edge["weight"]), edge["edge_uid"], edge))
                edges = sorted(edges)[: self.propagation_edge_limit()]
                denominator = sum(float(edge["weight"]) for _weight, _uid, edge in edges) or 1.0
                target_score = scores.get(target["node_uid"], 0.0)
                for _weight, _uid, edge in edges:
                    disease_idx = edge["target_idx"] if edge["source_idx"] == target_idx else edge["source_idx"]
                    disease_idx_scores[disease_idx] += alpha * target_score * float(edge["weight"]) / denominator
            disease_nodes = self.get_nodes_by_idx(list(disease_idx_scores))
            idx_to_uid = {int(node["node_idx"]): node["node_uid"] for node in disease_nodes}
            for disease_idx, score in disease_idx_scores.items():
                disease_uid = idx_to_uid.get(disease_idx)
                if disease_uid and score > 0.0:
                    disease_scores[disease_uid] += score
            return dict(disease_scores)

        for target in target_nodes:
            edges = []
            for edge in self.incident_edges(target["node_uid"]):
                neighbor = edge["object_uid"] if edge["subject_uid"] == target["node_uid"] else edge["subject_uid"]
                neighbor_type = edge["object_type"] if edge["subject_uid"] == target["node_uid"] else edge["subject_type"]
                if (target.get("node_type"), edge["edge_type"], neighbor_type) != ("target", "target_associated_with_disease", "disease"):
                    continue
                edges.append((-float(edge.get("p_final") or 0.0), edge["edge_uid"], neighbor, edge))
            edges = sorted(edges)[: self.propagation_edge_limit()]
            denominator = sum(float(edge.get("p_final") or 0.0) for _weight, _uid, _neighbor, edge in edges) or 1.0
            target_score = scores.get(target["node_uid"], 0.0)
            for _weight, _uid, neighbor, edge in edges:
                disease_scores[neighbor] += alpha * target_score * float(edge.get("p_final") or 0.0) / denominator
        return dict(disease_scores)

    def best_paths_by_terminal(self, paths: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
        best: dict[str, dict[str, Any]] = {}
        for path in paths:
            node_uids = path.get("node_uids") or []
            if not node_uids:
                continue
            terminal_uid = node_uids[-1]
            current = best.get(terminal_uid)
            if current is None or (path.get("cost", math.inf), path.get("path_id", "")) < (current.get("cost", math.inf), current.get("path_id", "")):
                best[terminal_uid] = path
        return best

    def literature_support_components_from_edges(self, edges: list[dict[str, Any]]) -> dict[str, Any]:
        support_by_uid: dict[str, dict[str, Any]] = {}
        pmids: set[str] = set()
        support_classes: set[str] = set()
        refs = []
        for edge in edges:
            edge_uid = edge.get("edge_uid", "")
            for support in edge.get("literature_support") or []:
                support_class = support.get("support_class", "")
                if support_class not in {"confirm", "support_direction", "novel_candidate"}:
                    continue
                support_uid = support.get("support_uid", "")
                if support_uid in support_by_uid:
                    continue
                support_by_uid[support_uid] = support
                p_literature = round(float(support.get("p_literature") or 0.0), 6)
                row_pmids = parse_list(support.get("pmids"))
                row_pmcids = parse_list(support.get("pmcids"))
                pmids.update(row_pmids)
                support_classes.add(support_class)
                refs.append(
                    {
                        "edge_uid": edge_uid,
                        "support_uid": support_uid,
                        "support_class": support_class,
                        "p_literature": p_literature,
                        "pmids": row_pmids[:5],
                        "pmcids": row_pmcids[:5],
                        "relation_uids": parse_list(support.get("evidence_relation_uids"))[:5],
                        "sentence_uids": parse_list(support.get("sentence_uids"))[:5],
                        "license_id": support.get("license_id", ""),
                        "source_release": support.get("source_release", ""),
                        "parser_hash": support.get("parser_hash", ""),
                        "config_hash": support.get("config_hash", ""),
                    }
                )
        refs.sort(key=lambda row: (-float(row.get("p_literature") or 0.0), row.get("edge_uid", ""), row.get("support_uid", "")))
        p_values = [float(row.get("p_literature") or 0.0) for row in support_by_uid.values()]
        return {
            "literature_support_count": len(support_by_uid),
            "max_p_literature": round(max(p_values), 6) if p_values else 0.0,
            "supported_pmids": sorted(pmids)[:20],
            "support_classes": sorted(support_classes),
            "evidence_refs": refs[:10],
        }

    def literature_support_components_from_path(self, path: dict[str, Any] | None, terminal_node_uid: str = "") -> dict[str, Any]:
        edge_rows = []
        seen_edge_uids: set[str] = set()
        for edge in (path or {}).get("edges", []):
            edge_rows.append(edge)
            if edge.get("edge_uid"):
                seen_edge_uids.add(edge["edge_uid"])
        if terminal_node_uid:
            self.load_literature_support()
            by_entity = self._literature_support_by_entity or {}
            novel_fallbacks = []
            for support in by_entity.get(terminal_node_uid, []):
                edge_uids = parse_list(support.get("supported_existing_edge_uids"))
                if edge_uids:
                    edge_uid = edge_uids[0]
                    if edge_uid in seen_edge_uids:
                        continue
                    edge_rows.append({"edge_uid": edge_uid, "literature_support": [support]})
                    seen_edge_uids.add(edge_uid)
                    continue
                if (
                    support.get("support_class") == "novel_candidate"
                    and float(support.get("p_literature") or 0.0) >= self.config.min_literature_overlay_prob
                ):
                    novel_fallbacks.append(support)
            novel_fallbacks.sort(key=lambda row: (-float(row.get("p_literature") or 0.0), row.get("support_uid", "")))
            for support in novel_fallbacks[:5]:
                edge = self.literature_overlay_edge_from_support(support)
                if edge and edge.get("edge_uid") not in seen_edge_uids:
                    edge_rows.append(edge)
                    seen_edge_uids.add(edge.get("edge_uid", ""))
        return self.literature_support_components_from_edges(edge_rows)

    def literature_ranking_boost(self, literature_components: dict[str, Any]) -> float:
        support_count = int(literature_components.get("literature_support_count") or 0)
        max_p_literature = float(literature_components.get("max_p_literature") or 0.0)
        if support_count <= 0 or max_p_literature <= 0.0:
            return 0.0
        return round(0.01 * max_p_literature * math.log2(support_count + 1.0), 12)

    def propagation_ranking(
        self,
        propagation: dict[str, Any],
        node_type: str,
        id_key: str,
        best_paths: dict[str, dict[str, Any]],
    ) -> list[dict[str, Any]]:
        rows = []
        for node in propagation.get("nodes", []):
            if node.get("node_type") != node_type:
                continue
            terminal_path = best_paths.get(node["node_uid"])
            literature_components = self.literature_support_components_from_path(terminal_path, node["node_uid"])
            propagation_score = float(node["propagation_score"])
            raw_literature_boost = self.literature_ranking_boost(literature_components)
            specificity_multiplier = self.ranking_specificity_multiplier(
                node_type,
                propagation_score,
                terminal_path,
                literature_components,
            )
            literature_boost = round(raw_literature_boost * specificity_multiplier, 12)
            row = {
                id_key: node["node_uid"],
                "display_name": node.get("display_name", ""),
                "primary_external_id": node.get("primary_external_id", ""),
                "score": round(propagation_score + literature_boost, 12),
                "score_components": {
                    "propagation_score": node["propagation_score"],
                    "literature_boost": literature_boost,
                    "raw_literature_boost": raw_literature_boost,
                    "specificity_multiplier": specificity_multiplier,
                    "seed_score": node["seed_score"],
                    "seed_weight": node["seed_weight"],
                    "best_path_confidence": terminal_path.get("path_confidence", 0.0) if terminal_path else 0.0,
                    "formula": "propagation_score + 0.01*max_p_literature*log2(literature_support_count+1)",
                    "literature_support_count": literature_components["literature_support_count"],
                    "max_p_literature": literature_components["max_p_literature"],
                },
                "literature_support_count": literature_components["literature_support_count"],
                "max_p_literature": literature_components["max_p_literature"],
                "supported_pmids": literature_components["supported_pmids"],
                "support_classes": literature_components["support_classes"],
                "evidence_refs": literature_components["evidence_refs"],
            }
            if terminal_path:
                row.update(
                    {
                        "best_path_id": terminal_path["path_id"],
                        "best_path_cost": terminal_path["cost"],
                        "path_confidence": terminal_path["path_confidence"],
                    }
                )
            rows.append(row)
        return sorted(rows, key=lambda row: (-row["score"], row.get("best_path_cost", math.inf), row[id_key]))

    def ranking_specificity_multiplier(
        self,
        node_type: str,
        propagation_score: float,
        terminal_path: dict[str, Any] | None,
        literature_components: dict[str, Any],
    ) -> float:
        support_count = int(literature_components.get("literature_support_count") or 0)
        if support_count <= 0:
            return 1.0
        if terminal_path and float(terminal_path.get("path_confidence") or 0.0) > 0.0:
            return 1.0
        if node_type == "disease":
            return round(max(0.12, min(0.45, 0.12 + math.sqrt(max(0.0, propagation_score)) * 8.0)), 6)
        if node_type == "target" and propagation_score < 0.00005:
            return 0.7
        return 1.0

    def merge_pathway_rankings(
        self,
        enrichment_rows: list[dict[str, Any]],
        propagation: dict[str, Any],
        best_paths: dict[str, dict[str, Any]],
        features_by_uid: dict[str, list[dict[str, Any]]],
        input_records: list[Any] | None = None,
    ) -> list[dict[str, Any]]:
        rows_by_uid = {row["pathway_uid"]: dict(row) for row in enrichment_rows}
        propagation_pathways = {node["node_uid"]: node for node in propagation.get("nodes", []) if node.get("node_type") == "pathway"}
        for pathway_uid, node in propagation_pathways.items():
            rows_by_uid.setdefault(
                pathway_uid,
                {
                    "pathway_uid": pathway_uid,
                    "name": node.get("display_name", ""),
                    "primary_external_id": node.get("primary_external_id", ""),
                    "source_name": "",
                    "species": "",
                    "overlap_count": 0,
                    "input_count": 0,
                    "pathway_size": 0,
                    "universe_size": 0,
                    "coverage": 0.0,
                    "p_value": 1.0,
                    "fdr": 1.0,
                    "matched_metabolite_uids": [],
                    "evidence_refs": [],
                    "score": 0.0,
                    "score_components": {"overlap_count": 0, "coverage": 0.0, "p_value": 1.0, "fdr": 1.0},
                },
            )
        theme_profile = self.input_theme_profile(features_by_uid, input_records=input_records)
        large_seed_set = len(features_by_uid) > 30
        if theme_profile and not large_seed_set:
            lexical_tokens = (theme_profile.get("_lexical") or {}).get("tokens") or []
            pathway_names = self.build_pathway_index().get("pathway_names", {})
            for pathway_uid, pathway in pathway_names.items():
                pathway_name = normalize_lookup_key(pathway.get("name") or "")
                if not pathway_name:
                    continue
                theme_match = False
                for theme_id, theme in theme_profile.items():
                    if theme_id == "_lexical":
                        if not lexical_tokens:
                            continue
                        if self.pathway_lexical_input_boost(pathway, theme_profile)[0] > 0.0:
                            theme_match = True
                            break
                        continue
                    if int(theme.get("support_count") or 0) < 2:
                        continue
                    spec = BIOCHEMICAL_THEME_SPECS.get(theme_id, {})
                    if any(normalized_text_has_term(pathway_name, term) for term in spec.get("pathway_terms", ())):
                        theme_match = True
                        break
                if not theme_match:
                    continue
                rows_by_uid.setdefault(
                    pathway_uid,
                    {
                        "pathway_uid": pathway_uid,
                        "name": pathway.get("name", ""),
                        "primary_external_id": pathway.get("primary_external_id", ""),
                        "source_name": pathway.get("source_name", ""),
                        "species": pathway.get("species", ""),
                        "overlap_count": 0,
                        "input_count": 0,
                        "pathway_size": 0,
                        "universe_size": 0,
                        "coverage": 0.0,
                        "p_value": 1.0,
                        "fdr": 1.0,
                        "matched_metabolite_uids": [],
                        "evidence_refs": [],
                        "score": 0.0,
                        "score_components": {
                            "overlap_count": 0,
                            "coverage": 0.0,
                            "p_value": 1.0,
                            "fdr": 1.0,
                            "theme_candidate": True,
                        },
                    },
                )
        rows = list(rows_by_uid.values())
        self.attach_pathway_user_support(rows, features_by_uid, input_records=input_records)
        for row in rows:
            row["enrichment_evidence_refs"] = list(row.get("evidence_refs", []))
            node = propagation_pathways.get(row["pathway_uid"], {})
            terminal_path = best_paths.get(row["pathway_uid"])
            literature_components = self.literature_support_components_from_path(terminal_path, row["pathway_uid"])
            enrichment_score = float(row.get("score", 0.0) or 0.0)
            propagation_score = float(node.get("propagation_score", 0.0) or 0.0)
            literature_boost = self.literature_ranking_boost(literature_components)
            row["score"] = round(enrichment_score + propagation_score + literature_boost, 12)
            row.setdefault("score_components", {})
            row["score_components"].update(
                {
                    "enrichment_score": round(enrichment_score, 12),
                    "propagation_score": round(propagation_score, 12),
                    "literature_boost": literature_boost,
                    "seed_score": node.get("seed_score", 0.0),
                    "best_path_confidence": terminal_path.get("path_confidence", 0.0) if terminal_path else 0.0,
                    "combined_score_formula": "enrichment_score + propagation_score + literature_boost",
                    "literature_support_count": literature_components["literature_support_count"],
                    "max_p_literature": literature_components["max_p_literature"],
                }
            )
            row["literature_support_count"] = literature_components["literature_support_count"]
            row["max_p_literature"] = literature_components["max_p_literature"]
            row["supported_pmids"] = literature_components["supported_pmids"]
            row["support_classes"] = literature_components["support_classes"]
            row["evidence_refs"] = literature_components["evidence_refs"]
            if terminal_path:
                row["best_path_id"] = terminal_path["path_id"]
                row["best_path_cost"] = terminal_path["cost"]
                row["path_confidence"] = terminal_path["path_confidence"]
        return sorted(rows, key=lambda row: (-row["score"], row.get("fdr", 1.0), row["pathway_uid"]))

    def path_search_incident_edges(self, node_uid: str, limit: int | None = None) -> list[dict[str, Any]]:
        node = self.get_node(node_uid)
        sparse_path = self.graph_dir / "sparse_edges.parquet"
        if node and node.get("node_idx") is not None and sparse_path.exists() and (self.graph_dir / "edge_type_index.parquet").exists():
            raw_edges = self.sparse_incident_edges(int(node["node_idx"]))
            if limit is not None:
                raw_edges = raw_edges[:limit]
            endpoint_idxs = sorted(
                {
                    int(idx)
                    for edge in raw_edges
                    for idx in (edge.get("source_idx"), edge.get("target_idx"))
                    if idx is not None
                }
            )
            idx_to_node = {int(row["node_idx"]): row for row in self.get_nodes_by_idx(endpoint_idxs)}
            edges: list[dict[str, Any]] = []
            for edge in raw_edges:
                source = idx_to_node.get(int(edge.get("source_idx")))
                target = idx_to_node.get(int(edge.get("target_idx")))
                if not source or not target:
                    continue
                p_final = clamp_probability(edge.get("p_final"), edge.get("weight", 0.5))
                edges.append(
                    {
                        **edge,
                        "subject_uid": source.get("node_uid", ""),
                        "object_uid": target.get("node_uid", ""),
                        "subject_type": edge.get("subject_type", source.get("node_type", "")),
                        "object_type": edge.get("object_type", target.get("node_type", "")),
                        "p_final": p_final,
                        "weight": p_final,
                        "cost": round(-math.log(p_final), 6),
                        "source_table": edge.get("source_table", "sparse_edges"),
                    }
                )
            return edges
        return self.incident_edges(node_uid, limit=limit)

    def enumerate_paths(
        self,
        seed_ids: list[str],
        max_hops: int | None = None,
        max_paths: int | None = None,
        target_types: set[str] | None = None,
    ) -> list[dict[str, Any]]:
        max_hops = max_hops or self.config.max_hops
        max_paths = max_paths or self.config.max_paths
        target_types = target_types or TARGET_ANALYSIS_TYPES
        heap: list[tuple[float, int, str, list[str], list[dict[str, Any]]]] = []
        counter = 0
        seed_set = set(seed_ids)
        for seed_id in sorted(seed_set):
            heapq.heappush(heap, (0.0, counter, seed_id, [seed_id], []))
            counter += 1
        paths = []
        expansions = 0
        while heap and len(paths) < max_paths and expansions < 2000:
            cost, _order, current, node_path, edge_path = heapq.heappop(heap)
            expansions += 1
            current_node = self.get_node(current)
            if edge_path and current_node and current_node.get("node_type") in target_types:
                paths.append(self.path_contract(node_path, edge_path, cost))
                if len(paths) >= max_paths:
                    break
            if len(edge_path) >= max_hops:
                continue
            for edge in self.path_search_incident_edges(current, limit=self.config.max_edges_per_node):
                if edge.get("edge_type") not in ALLOWED_EXPLANATION_EDGE_TYPES:
                    continue
                neighbor = edge["object_uid"] if edge["subject_uid"] == current else edge["subject_uid"]
                if neighbor in node_path:
                    continue
                next_node = self.get_node(neighbor)
                if not next_node:
                    continue
                next_cost = cost + float(edge["cost"])
                heapq.heappush(heap, (next_cost, counter, neighbor, node_path + [neighbor], edge_path + [edge]))
                counter += 1
        return paths

    def path_contract(self, node_path: list[str], edge_path: list[dict[str, Any]], cost: float) -> dict[str, Any]:
        nodes = self.get_nodes(node_path)
        edge_contracts = [self.edge_contract(edge) for edge in edge_path]
        return {
            "path_id": content_hash({"nodes": node_path, "edges": [edge["edge_uid"] for edge in edge_path]})[:20],
            "cost": round(cost, 6),
            "path_confidence": round(math.exp(-cost), 12),
            "node_uids": node_path,
            "nodes": nodes,
            "edges": edge_contracts,
            "score_components": {
                "cost_formula": "sum(-log(p_final(edge)))",
                "edge_p_final": [edge["p_final"] for edge in edge_contracts],
            },
        }

    def compress_explanation_paths(self, paths: list[dict[str, Any]]) -> dict[str, Any]:
        groups: dict[str, dict[str, Any]] = {}
        for path in paths:
            node_uids = path.get("node_uids") or []
            ordered_nodes = [self.get_node(uid) for uid in node_uids]
            node_types = [str(node.get("node_type") or "unknown") for node in ordered_nodes if node]
            edge_types = [edge.get("edge_type", "") for edge in path.get("edges", [])]
            signature = " -> ".join(node_types)
            edge_signature = " -> ".join(edge_types)
            key = content_hash({"node_types": node_types, "edge_types": edge_types})[:20]
            group = groups.setdefault(
                key,
                {
                    "signature_id": key,
                    "node_type_signature": signature,
                    "edge_type_signature": edge_signature,
                    "path_count": 0,
                    "aggregated_confidence": 0.0,
                    "representative_path": None,
                    "_confidence_remaining": 1.0,
                    "_intermediate_nodes": Counter(),
                },
            )
            group["path_count"] += 1
            confidence = clamp_unit(path.get("path_confidence", 0.0))
            group["_confidence_remaining"] *= 1.0 - confidence
            current_representative = group.get("representative_path")
            if current_representative is None or (
                path.get("cost", math.inf),
                path.get("path_id", ""),
            ) < (
                current_representative.get("cost", math.inf),
                current_representative.get("path_id", ""),
            ):
                group["representative_path"] = {
                    "path_id": path.get("path_id", ""),
                    "cost": path.get("cost", 0.0),
                    "path_confidence": path.get("path_confidence", 0.0),
                    "node_uids": node_uids,
                }
            for node in ordered_nodes[1:-1]:
                if not node:
                    continue
                label = node.get("display_name") or node.get("node_uid", "")
                group["_intermediate_nodes"][f"{node.get('node_type', '')}:{label}"] += 1

        compressed = []
        for group in groups.values():
            confidence = 1.0 - float(group.pop("_confidence_remaining"))
            intermediates = group.pop("_intermediate_nodes")
            group["aggregated_confidence"] = round(confidence, 12)
            group["top_intermediate_nodes"] = [
                {"node": label, "count": count}
                for label, count in sorted(intermediates.items(), key=lambda item: (-item[1], item[0]))[:10]
            ]
            group["omitted_count"] = max(0, group["path_count"] - 1)
            compressed.append(group)
        compressed.sort(
            key=lambda row: (
                -row["aggregated_confidence"],
                row["representative_path"]["cost"] if row.get("representative_path") else math.inf,
                row["signature_id"],
            )
        )
        return {
            "mode": "path_signature_compression",
            "input_path_count": len(paths),
            "signature_count": len(compressed),
            "signatures": compressed,
        }

    def analysis_pack_candidate(self, candidate: dict[str, Any]) -> dict[str, Any]:
        return {
            "entity_uid": candidate.get("entity_uid", ""),
            "entity_type": candidate.get("entity_type", ""),
            "display_name": candidate.get("display_name", ""),
            "primary_external_id": candidate.get("primary_external_id", ""),
            "score": candidate.get("score", 0.0),
            "score_components": candidate.get("score_components", {}),
        }

    def analysis_pack_expanded_seed_index(self, precheck: dict[str, Any]) -> dict[str, dict[str, Any]]:
        index: dict[str, dict[str, Any]] = {}
        for row in precheck.get("expanded_candidates", []) or []:
            expanded = row.get("expanded_candidate") or {}
            original_input_id = str(expanded.get("original_input_id") or row.get("input_id") or "")
            if not original_input_id:
                continue
            entry = index.setdefault(
                original_input_id,
                {
                    "candidate_count": 0,
                    "seed_classes": Counter(),
                    "seed_weight_sum": 0.0,
                    "reasons": Counter(),
                },
            )
            entry["candidate_count"] += 1
            seed_class = str(expanded.get("seed_class") or "expanded")
            reason = str(expanded.get("reason") or "")
            entry["seed_classes"][seed_class] += 1
            if reason:
                entry["reasons"][reason] += 1
            entry["seed_weight_sum"] += float(expanded.get("weight_multiplier") or 0.0)
        return {
            input_id: {
                "candidate_count": entry["candidate_count"],
                "seed_classes": dict(sorted(entry["seed_classes"].items())),
                "seed_weight_sum": round(entry["seed_weight_sum"], 6),
                "reasons": dict(sorted(entry["reasons"].items())),
            }
            for input_id, entry in index.items()
        }

    def analysis_pack_genetic_exposure_row(
        self,
        row: dict[str, Any],
        resolution_status: str,
        expanded_seed: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        evidence = self.european_trait_source_evidence_for_record(row.get("record", {}))
        if not evidence.get("accession_id") and not evidence.get("summary_statistics_url"):
            return {}
        if resolution_status == "matched":
            chemical_identity_track = "strict_identity"
            seed_classes = {"strict_identity": 1}
            seed_weight_sum = 1.0
        elif expanded_seed:
            chemical_identity_track = "expanded_mechanism_seed"
            seed_classes = expanded_seed.get("seed_classes", {})
            seed_weight_sum = expanded_seed.get("seed_weight_sum", 0.0)
        else:
            chemical_identity_track = "review_only"
            seed_classes = {}
            seed_weight_sum = 0.0
        return {
            "input_id": row.get("input_id", ""),
            "exposure_id": evidence.get("accession_id", ""),
            "trait_description": evidence.get("reported_trait", ""),
            "evidence_type": "gwas_summary_statistics",
            "genetic_exposure_status": "available" if evidence.get("summary_statistics_url") else "metadata_only",
            "genetic_confidence_scope": "high_for_gcst_exposure_not_chemical_identity",
            "allowed_uses": [
                "genetic_overlap",
                "lead_locus_mapping",
                "colocalization",
                "mendelian_randomization_exposure",
            ],
            "identity_boundary": (
                "GCST summary statistics define a GWAS phenotype/exposure; they do not by themselves "
                "assign a unique HMDB/PubChem/ChEBI chemical identity."
            ),
            "chemical_identity_track": chemical_identity_track,
            "chemical_resolution_status": resolution_status,
            "chemical_seed_classes": seed_classes,
            "chemical_seed_weight_sum": seed_weight_sum,
            "expanded_candidate_count": int((expanded_seed or {}).get("candidate_count") or 0),
            "source_trait_evidence": evidence,
        }

    def analysis_pack_genetic_exposures(self, precheck: dict[str, Any]) -> list[dict[str, Any]]:
        expanded_index = self.analysis_pack_expanded_seed_index(precheck)
        exposures: list[dict[str, Any]] = []
        for bucket, status in (
            ("matched", "matched"),
            ("ambiguous", "ambiguous"),
            ("unmatched", "unmatched"),
            ("invalid", "invalid"),
        ):
            for row in precheck.get(bucket, []) or []:
                input_id = str(row.get("input_id") or "")
                exposure = self.analysis_pack_genetic_exposure_row(row, status, expanded_index.get(input_id))
                if exposure:
                    exposures.append(exposure)
        exposures.sort(key=lambda row: (str(row.get("exposure_id") or ""), str(row.get("input_id") or "")))
        return exposures

    def analysis_pack_four_track_summary(
        self,
        input_summary: dict[str, Any],
        genetic_exposures: list[dict[str, Any]],
        expanded_summary: dict[str, Any],
    ) -> dict[str, Any]:
        by_chemical_track = Counter(row.get("chemical_identity_track", "unknown") for row in genetic_exposures)
        return {
            "genetic_exposure_count": len(genetic_exposures),
            "genetic_exposure_available_count": sum(
                1 for row in genetic_exposures if row.get("genetic_exposure_status") == "available"
            ),
            "strict_identity_count": int(input_summary.get("strict_matched_count") or input_summary.get("matched_count") or 0),
            "expanded_mechanism_input_count": int(
                input_summary.get("expanded_candidate_input_count") or input_summary.get("expanded_input_count") or 0
            ),
            "expanded_mechanism_candidate_count": int(input_summary.get("expanded_candidate_count") or 0),
            "unresolved_chemical_identity_count": int(input_summary.get("unresolved_input_count") or 0),
            "analysis_seed_input_count": int(input_summary.get("analysis_seed_input_count") or 0),
            "weighted_seed_mass": input_summary.get("weighted_seed_mass", 0.0),
            "genetic_exposure_by_chemical_track": dict(sorted(by_chemical_track.items())),
            "expanded_by_class": expanded_summary.get("expanded_by_class", {}),
            "interpretation": {
                "genetic_exposure": "Use GCST-level summary statistics for genetic overlap, locus mapping, colocalization, and MR-style exposure analyses.",
                "mechanism_seed": "Use reportedTrait-derived strict/soft/ratio/class seeds for pathway and mechanism hypotheses with class-specific weights.",
                "strict_identity": "Use only strict identities for exact compound details, exact literature claims, and stronger pathway/target wording.",
                "review_only": "Keep unresolved, platform-code, isomeric, or pool-like traits out of exact entity scoring until curated identifiers are supplied.",
            },
        }

    def analysis_pack_resolution_row(self, row: dict[str, Any], status: str) -> dict[str, Any]:
        resolution = row.get("resolution", {})
        candidates = resolution.get("candidates") or []
        expanded = row.get("expanded_candidate") or {}
        match_status = "matched" if status == "matched" else ("ambiguous" if status == "ambiguous" else "unmatched")
        if status == "expanded_candidate" or expanded:
            match_status = "expanded"
        contract = {
            "input_id": row.get("input_id", ""),
            "query": row.get("query", ""),
            "raw_input": row.get("query", ""),
            "query_field": row.get("query_field", ""),
            "record": row.get("record", {}),
            "resolution_status": status,
            "match_status": match_status,
            "top_score": resolution.get("top_score", 0.0),
            "match_score": resolution.get("top_score", 0.0),
            "top_margin": resolution.get("top_margin", 0.0),
            "margin": resolution.get("top_margin", 0.0),
            "identity_review_reasons": resolution.get("identity_review_reasons", []),
            "query_notes": normalize_query_notes(resolution.get("query_notes", [])),
            "candidates": [self.analysis_pack_candidate(candidate) for candidate in candidates[:5]],
            "candidate_entities": [self.analysis_pack_candidate(candidate) for candidate in candidates[:5]],
            "biological_interpretation_entities": self.biological_interpretation_entities_for_resolution_row(row),
            "identity_resolution_v2": row.get("identity_resolution_v2") or self.identity_resolution_v2_for_row(row),
        }
        if expanded:
            contract["expanded_candidate"] = dict(expanded)
            contract["seed_track"] = expanded.get("track", "expanded")
            contract["seed_class"] = expanded.get("seed_class", "")
            contract["seed_weight_multiplier"] = expanded.get("weight_multiplier", 0.0)
        else:
            contract["seed_track"] = "strict" if status == "matched" else "review"
            contract["seed_class"] = "strict_identity" if status == "matched" else ""
            contract["seed_weight_multiplier"] = 1.0 if status == "matched" else 0.0
        if row.get("metabolite_uid"):
            contract["metabolite_uid"] = row.get("metabolite_uid", "")
            contract["normalized_entity_id"] = row.get("metabolite_uid", "")
        else:
            contract["normalized_entity_id"] = ""
        contract["biological_entity_pool"] = contract["biological_interpretation_entities"]
        source_trait_evidence = self.european_trait_source_evidence_for_record(row.get("record", {}))
        if source_trait_evidence:
            contract["source_trait_evidence"] = source_trait_evidence
        if row.get("invalid_reason"):
            contract["invalid_reason"] = row.get("invalid_reason", "")
        return contract

    def analysis_pack_quality_warnings(
        self,
        precheck: dict[str, Any],
        analysis_features: list[dict[str, Any]],
        features_by_uid: dict[str, list[dict[str, Any]]],
        propagation: dict[str, Any],
        paths: list[dict[str, Any]],
    ) -> list[dict[str, Any]]:
        warnings = []

        def add(code: str, severity: str, message: str, details: dict[str, Any]) -> None:
            warnings.append({"code": code, "severity": severity, "message": message, "details": details})

        if not self.database_accuracy_store_available():
            add(
                "database_accuracy_store_missing",
                "info",
                "Database accuracy v2 fact store is not enabled for this release; legacy graph rankings remain available.",
                {"release_id": self.release_id, "expected_manifest": str(self.database_accuracy_manifest_path)},
            )

        summary = precheck.get("summary", {})
        for warning in (precheck.get("input_normalization") or {}).get("warnings", []) or []:
            if isinstance(warning, dict):
                add(
                    str(warning.get("code") or "input_normalization_warning"),
                    str(warning.get("severity") or "warning"),
                    str(warning.get("message") or "Input normalization warning."),
                    {key: value for key, value in warning.items() if key not in {"code", "severity", "message"}},
                )
        if int(summary.get("ambiguous", 0) or 0):
            add(
                "ambiguous_inputs_present",
                "warning",
                "Some input records did not pass strict identity resolution; eligible candidates are used only as low-weight expanded seeds.",
                {
                    "count": int(summary.get("ambiguous", 0) or 0),
                    "expanded_summary": precheck.get("expanded_summary", {}),
                },
            )
        if int(summary.get("unmatched", 0) or 0):
            add(
                "unmatched_inputs_present",
                "warning",
                "Some input records could not be matched to the frozen release.",
                {"count": int(summary.get("unmatched", 0) or 0)},
            )
        if int(summary.get("invalid", 0) or 0):
            add(
                "invalid_inputs_present",
                "warning",
                "Some input records did not contain a supported identifier or feature set.",
                {"count": int(summary.get("invalid", 0) or 0)},
            )

        duplicates = [
            {"metabolite_uid": metabolite_uid, "input_ids": sorted(feature.get("input_id", "") for feature in features)}
            for metabolite_uid, features in sorted(features_by_uid.items())
            if len(features) > 1
        ]
        if duplicates:
            add(
                "duplicate_matched_metabolites",
                "warning",
                "Multiple input records resolved to the same metabolite and were retained as separate seed support.",
                {"duplicates": duplicates[:20], "duplicate_metabolite_count": len(duplicates)},
            )

        feature_summary = self.feature_summary(analysis_features)
        if analysis_features and not feature_summary.get("with_fold_change") and not feature_summary.get("with_p_value"):
            add(
                "metabolite_list_mode",
                "info",
                "No differential statistics were detected; rankings use matched metabolite presence as seed support.",
                {"matched_records": len(analysis_features)},
            )

        graph = propagation.get("graph", {})
        truncated = graph.get("truncated", {})
        if truncated.get("nodes") or truncated.get("edges"):
            add(
                "propagation_graph_truncated",
                "warning",
                "Propagation hit graph size limits.",
                {"truncated": truncated},
            )
        compression = graph.get("compression", {})
        retained_mass = maybe_number(compression.get("retained_mass"))
        if compression.get("quality_status") == "lossy_below_threshold" or (
            retained_mass is not None and retained_mass < float(self.config.propagation_retained_mass)
        ):
            add(
                "propagation_compression_lossy",
                "warning",
                "Typed beam compression retained less graph mass than the configured threshold.",
                {
                    "retained_mass": compression.get("retained_mass"),
                    "retained_mass_threshold": self.config.propagation_retained_mass,
                },
            )
        if analysis_features and not paths:
            add(
                "no_explanation_paths",
                "warning",
                "Matched seeds produced no stable stepwise graph paths under the current hop and confidence limits; ranking evidence fallback chains may still be shown as low-confidence audit support.",
                {
                    "max_hops": self.config.max_hops,
                    "max_paths": self.config.max_paths,
                    "fallback_chain_policy": "ranking_evidence_fallback_only_when_top_explanation_paths_are_empty",
                },
            )
        return sorted(warnings, key=lambda row: (row["severity"], row["code"]))

    def analysis_pack_blocked_reasons(
        self,
        precheck: dict[str, Any],
        pathways: list[dict[str, Any]],
        targets: list[dict[str, Any]],
        diseases: list[dict[str, Any]],
    ) -> list[dict[str, Any]]:
        reasons = []
        matched_count = int(precheck.get("summary", {}).get("matched", 0) or 0)
        if matched_count == 0:
            reasons.append(
                {
                    "code": "no_matched_metabolites",
                    "message": "Analysis cannot score rankings until at least one input resolves to a metabolite.",
                }
            )
        elif not pathways and not targets and not diseases:
            reasons.append(
                {
                    "code": "empty_rankings",
                    "message": "Matched inputs did not connect to pathway, target, or disease ranking claims in this release.",
                }
            )
        return reasons

    def database_accuracy_identity_rows(self, precheck: dict[str, Any]) -> list[dict[str, Any]]:
        rows = []
        for bucket in ("matched", "ambiguous", "unmatched", "invalid"):
            for row in precheck.get(bucket, []) or []:
                identity = row.get("identity_resolution_v2") or {}
                if identity:
                    rows.append({"bucket": bucket, **identity})
        return rows

    def database_accuracy_mechanism_facts(
        self,
        accepted_exact_uids: set[str],
        limit: int = 50,
        context: dict[str, Any] | None = None,
    ) -> list[dict[str, Any]]:
        if not self.database_accuracy_store_available() or not accepted_exact_uids:
            return []
        facts = []
        for row in self.database_accuracy_rows("analysis_view", "mechanism_ready_facts"):
            if str(row.get("subject_uid") or "") not in accepted_exact_uids:
                continue
            metadata = safe_metadata(row.get("metadata_json"))
            source_record_uids = parse_list(row.get("source_record_uids"))
            assertion_uids = parse_list(row.get("supporting_assertion_uids"))
            fact_payload = {
                    "fact_uid": row.get("fact_uid", ""),
                    "fact_type": row.get("fact_type", ""),
                    "subject_uid": row.get("subject_uid", ""),
                    "subject_type": row.get("subject_type", ""),
                    "subject_name": row.get("subject_name", "") or metadata.get("physical_entity_name", ""),
                    "predicate": row.get("predicate", ""),
                    "object_uid": row.get("object_uid", ""),
                    "object_type": row.get("object_type", ""),
                    "object_name": row.get("object_name", "") or metadata.get("reaction_name", ""),
                    "direction": row.get("direction", ""),
                    "role": row.get("role", ""),
                    "readiness_tier": row.get("readiness_tier", ""),
                    "allowed_claim_scope": row.get("allowed_claim_scope", ""),
                    "mechanism_scope": row.get("allowed_claim_scope", ""),
                    "support_status": "support",
                    "evidence_assertion_uid": assertion_uids[0] if assertion_uids else "",
                    "supporting_assertion_uids": assertion_uids,
                    "source_record_uid": source_record_uids[0] if source_record_uids else "",
                    "source_record_uids": source_record_uids,
                    "identity_decision_uids": parse_list(row.get("identity_decision_uids")),
                    "boundary_text": row.get("boundary_text", ""),
                    "metadata": metadata,
                }
            display_text = " ".join(
                str(fact_payload.get(key) or "")
                for key in ("subject_name", "predicate", "object_name", "role", "boundary_text")
            )
            context_eval = self.prediction_context_mismatch(display_text, context, "mechanism_fact")
            if context_eval.get("context_mismatch"):
                fact_payload.update(
                    {
                        "appendix": True,
                        "appendix_reason": context_eval.get("context_penalty_reason", "mechanism_fact_context_mismatch"),
                        **context_eval,
                    }
                )
            else:
                fact_payload.setdefault("appendix", False)
                fact_payload.setdefault("context_mismatch", False)
            facts.append(fact_payload)
        scope_rank = {
            "directional_reaction_fact": 0,
            "enzyme_reaction_fact": 1,
            "transporter_reaction_fact": 1,
            "bidirectional_reaction_fact": 2,
            "exact_reaction_fact": 3,
            "role_unknown_reaction_fact": 4,
        }
        relation_source_rank = {
            "rhea_rdf_explicit_side": 0,
            "rhea_biopax_explicit_side": 0,
            "reactome_sbml_species_reference": 1,
            "reactome_biopax_conversion": 1,
            "rhea_smiles_exact_component": 2,
            "reactome_pe_participation": 4,
        }
        return sorted(
            facts,
            key=lambda row: (
                scope_rank.get(str(row.get("allowed_claim_scope") or ""), 9),
                relation_source_rank.get(str((row.get("metadata") or {}).get("relation_source") or ""), 9),
                str(row.get("subject_uid") or ""),
                str(row.get("object_uid") or ""),
                str(row.get("fact_uid") or ""),
            ),
        )[:limit]

    def database_accuracy_pack(
        self,
        precheck: dict[str, Any],
        features_by_uid: dict[str, list[dict[str, Any]]],
        context: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        identity_rows = self.database_accuracy_identity_rows(precheck)
        decisions = [row.get("decision") or {} for row in identity_rows]
        decision_summary = Counter(str(decision.get("decision_status") or "unknown") for decision in decisions)
        feature_summary = Counter(
            str((row.get("input_feature") or {}).get("feature_type") or "unknown")
            for row in identity_rows
        )
        accepted_exact_uids = {
            str(decision.get("accepted_entity_uid") or "")
            for decision in decisions
            if decision.get("decision_status") == "accepted_exact" and decision.get("accepted_entity_uid")
        }
        accepted_exact_uids.update(features_by_uid)
        trait_observations = []
        class_observations = []
        exact_input_observations = []
        blocked_inputs = []
        for row in identity_rows:
            feature = row.get("input_feature") or {}
            decision = row.get("decision") or {}
            accepted_uid = str(decision.get("accepted_entity_uid") or "")
            accepted_candidate = next(
                (
                    candidate
                    for candidate in row.get("candidates", []) or []
                    if str(candidate.get("candidate_entity_uid") or "") == accepted_uid
                ),
                {},
            )
            accepted_name = str(accepted_candidate.get("candidate_name") or "").strip()
            observation = {
                "input_id": feature.get("input_row_id", ""),
                "feature_label": feature.get("feature_label", ""),
                "raw_feature_label": feature.get("feature_label", ""),
                "feature_type": feature.get("feature_type", ""),
                "decision_status": decision.get("decision_status", ""),
                "accepted_entity_uid": accepted_uid,
                "accepted_entity_name": accepted_name,
                "decision_rule": decision.get("decision_rule", ""),
                "direction": feature.get("direction", "unknown"),
                "effect_value": feature.get("effect_value"),
                "effect_label": feature.get("effect_label", ""),
                "allowed_claim_scope": "input_observation_only",
                "appendix": True,
            }
            if decision.get("decision_status") == "accepted_trait":
                trait_observations.append(observation)
            elif decision.get("decision_status") == "accepted_class":
                class_observations.append(
                    {
                        **observation,
                        "allowed_claim_scope": "input_class_observation",
                        "appendix": False,
                        "input_trace": True,
                    }
                )
            elif decision.get("decision_status") == "accepted_exact":
                exact_input_observations.append(
                    {
                        **observation,
                        "feature_label": accepted_name or observation["feature_label"],
                        "allowed_claim_scope": "input_exact_observation",
                        "appendix": False,
                        "input_trace": True,
                    }
                )
            elif decision.get("decision_status") in {"ambiguous", "unmatched", "rejected"}:
                blocked_inputs.append({**observation, "rejected_reasons": decision.get("rejected_reasons", [])})
        evidence_summary: dict[str, Any] = {}
        if self.database_accuracy_store_available():
            evidence_rows = self.database_accuracy_rows("evidence_store", "evidence_assertions")
            evidence_summary = {
                "assertion_count": len(evidence_rows),
                "by_support_status": dict(sorted(Counter(str(row.get("support_status") or "unknown") for row in evidence_rows).items())),
                "by_evidence_class": dict(sorted(Counter(str(row.get("evidence_class") or "unknown") for row in evidence_rows).items())),
            }
        manifest = self.database_accuracy_manifest()
        return {
            "contract_version": DATABASE_ACCURACY_CONTRACT_VERSION,
            "store_available": self.database_accuracy_store_available(),
            "manifest": {
                "release_id": manifest.get("release_id", self.release_id),
                "table_count": manifest.get("table_count", 0),
                "total_rows": manifest.get("total_rows", 0),
                "parser_hash": manifest.get("parser_hash", ""),
                "manifest_hash": manifest.get("manifest_hash", ""),
            },
            "identity_decision_summary": dict(sorted(decision_summary.items())),
            "input_feature_type_summary": dict(sorted(feature_summary.items())),
            "mechanism_ready_facts": self.database_accuracy_mechanism_facts(accepted_exact_uids, context=context),
            "exact_input_observations": exact_input_observations[:200],
            "trait_observations": trait_observations[:100],
            "class_observations": class_observations[:100],
            "blocked_inputs": blocked_inputs[:200],
            "evidence_assertion_summary": evidence_summary,
        }

    def gold_standard_sources(self) -> list[Path]:
        return [self.workspace / relative_path for relative_path in DEFAULT_GOLD_STANDARD_PATHS]

    def load_gold_standard_conclusions(self) -> dict[str, Any]:
        if self._gold_standard_cache is not None:
            return self._gold_standard_cache
        entries: list[dict[str, Any]] = []
        files: list[str] = []
        for path in self.gold_standard_sources():
            if not path.exists():
                continue
            files.append(str(path))
            if path.suffix.casefold() == ".csv":
                raw_rows = load_records_file(path)
            else:
                parsed = read_json_file(path)
                if isinstance(parsed, dict):
                    raw_rows = parsed.get("conclusions", parsed.get("entries", []))
                else:
                    raw_rows = parsed
            if not isinstance(raw_rows, list):
                continue
            for index, row in enumerate(raw_rows, start=1):
                if not isinstance(row, dict):
                    continue
                entry = self.normalize_gold_standard_entry(row, path, index)
                if entry:
                    entries.append(entry)
        entries.sort(key=lambda row: (str(row.get("gold_id") or ""), str(row.get("_source_file") or ""), int(row.get("_source_row") or 0)))
        self._gold_standard_cache = {
            "contract_version": CONCLUSION_EVALUATION_CONTRACT_VERSION,
            "files": sorted(set(files)),
            "entries": entries,
            "entry_count": len(entries),
        }
        return self._gold_standard_cache

    def normalize_gold_standard_entry(self, row: dict[str, Any], path: Path, row_number: int) -> dict[str, Any]:
        gold_id = str(row.get("gold_id") or row.get("id") or "").strip()
        if not gold_id:
            gold_id = f"gold_{content_hash({'path': str(path), 'row': row_number, 'payload': row})[:16]}"
        polarity = normalize_lookup_key(row.get("polarity") or row.get("label") or "positive")
        is_negative = bool(row.get("is_negative_control") is True or str(row.get("is_negative_control") or "").casefold() in {"1", "true", "yes"})
        is_background = polarity in {"background", "context", "review"}
        context_terms = []
        for key in ("cancer_type", "disease", "tissue", "cell_type", "cell_state", "comparison", "species"):
            context_terms.extend(prediction_terms_from_value(row.get(key)))
        context_terms.extend(prediction_terms_from_value(row.get("context_terms")))
        expected_entities = [
            term
            for term in prediction_terms_from_value(row.get("expected_entities") or row.get("entities"))
            if term not in GOLD_GENERIC_TERMS
        ]
        mechanism_terms = [
            term
            for term in prediction_terms_from_value(row.get("mechanism_axis") or row.get("mechanism") or row.get("pathway") or row.get("theme"))
            if term not in GOLD_GENERIC_TERMS
        ]
        trap_terms = exact_prediction_terms_from_value(row.get("negative_trap_terms") or row.get("trap_terms"))
        if is_negative and not trap_terms:
            trap_terms = sorted(
                term
                for term in set(
                    exact_prediction_terms_from_value(row.get("mechanism_axis") or row.get("mechanism") or row.get("pathway") or row.get("theme"))
                    + exact_prediction_terms_from_value(row.get("expected_entities") or row.get("entities"))
                )
                if term not in GOLD_GENERIC_TERMS
            )
        if not mechanism_terms and not expected_entities and not trap_terms:
            return {}
        return {
            "gold_id": gold_id,
            "claim_type": str(row.get("claim_type") or "mechanism_axis"),
            "mechanism_axis": str(row.get("mechanism_axis") or row.get("mechanism") or row.get("pathway") or row.get("theme") or ""),
            "mechanism_terms": mechanism_terms,
            "expected_entities": expected_entities,
            "expected_direction": normalize_lookup_key(row.get("expected_direction") or row.get("direction") or ""),
            "required_evidence_type": str(row.get("required_evidence_type") or row.get("evidence_type") or ""),
            "polarity": "negative" if is_negative else ("background" if is_background else "positive"),
            "is_negative_control": is_negative,
            "context_terms": sorted(set(context_terms)),
            "negative_trap_terms": trap_terms,
            "source_refs": parse_list(row.get("source_refs") or row.get("source_pmid_or_doi") or row.get("pmid") or row.get("doi")),
            "notes": str(row.get("notes") or ""),
            "_source_file": str(path),
            "_source_row": row_number,
        }

    def conclusion_candidate_rows(self, analysis_pack: dict[str, Any]) -> list[dict[str, Any]]:
        rows: list[dict[str, Any]] = []

        def add_candidate(source: str, label: str, payload: dict[str, Any], rank: int = 0) -> None:
            text_parts = [label]
            for key in ("result_type", "boundary", "predicate", "subject_name", "object_name", "direction", "role"):
                if payload.get(key) not in {None, ""}:
                    text_parts.append(str(payload.get(key)))
            allowed_claim_scope = str(payload.get("allowed_claim_scope") or "")
            observation_only = source in {"trait_observation"} or allowed_claim_scope == "input_observation_only"
            rows.append(
                {
                    "candidate_id": f"{source}:{payload.get('prediction_id') or payload.get('fact_uid') or payload.get('pathway_uid') or payload.get('target_uid') or payload.get('disease_uid') or rank}",
                    "source": source,
                    "rank": rank,
                    "label": label,
                    "normalized_text": normalize_lookup_key(" ".join(text_parts)),
                    "has_evidence_refs": bool(payload.get("evidence_refs") or payload.get("source_record_uids") or payload.get("supporting_assertion_uids")),
                    "has_mechanism_trace": bool(payload.get("fact_uid") or payload.get("claim_refs") or payload.get("source_record_uids")),
                    "has_input_trace": bool(payload.get("input_trace")),
                    "appendix": bool(payload.get("appendix") or payload.get("context_mismatch") or observation_only),
                    "allowed_claim_scope": allowed_claim_scope or ("input_observation_only" if observation_only else ""),
                    "context_mismatch": bool(payload.get("context_mismatch")),
                    "payload": payload,
                }
            )

        for source, key in (("pathway_ranking", "pathway_rankings"), ("target_ranking", "target_rankings"), ("disease_ranking", "disease_rankings")):
            for rank, row in enumerate(analysis_pack.get(key, []) or [], start=1):
                label = str(row.get("label") or row.get("display_name") or row.get("name") or "")
                add_candidate(source, label, row, rank)
        for rank, fact in enumerate(((analysis_pack.get("database_accuracy") or {}).get("mechanism_ready_facts") or []), start=1):
            label = " ".join(str(fact.get(key) or "") for key in ("subject_name", "predicate", "object_name"))
            add_candidate("mechanism_ready_fact", label, fact, rank)
        for rank, row in enumerate(
            ((analysis_pack.get("database_accuracy") or {}).get("exact_input_observations") or [])[:CONCLUSION_OBSERVATION_CANDIDATE_LIMIT],
            start=1,
        ):
            add_candidate("input_exact_observation", str(row.get("feature_label") or ""), row, rank)
        for source, key in (("trait_observation", "trait_observations"), ("class_observation", "class_observations")):
            for rank, row in enumerate(
                ((analysis_pack.get("database_accuracy") or {}).get(key) or [])[:CONCLUSION_OBSERVATION_CANDIDATE_LIMIT],
                start=1,
            ):
                add_candidate(source, str(row.get("feature_label") or ""), row, rank)
        return rows

    def gold_context_status(self, gold: dict[str, Any], context: dict[str, Any] | None) -> str:
        required = set(gold.get("context_terms") or [])
        if not required:
            return "no_context_requirement"
        requested = set((context or {}).get("terms", []) or [])
        if not requested:
            return "context_not_provided"
        matched = {term for term in required if any(normalized_text_has_term(ctx, term) for ctx in requested)}
        if matched and matched <= BROAD_CONTEXT_TERMS and any(term not in BROAD_CONTEXT_TERMS for term in required):
            return "context_mismatch"
        if matched == required:
            return "exact_context"
        if matched:
            return "compatible_context"
        return "context_mismatch"

    def gold_context_match_terms(self, gold: dict[str, Any], context: dict[str, Any] | None) -> list[str]:
        required = sorted(set(gold.get("context_terms") or []))
        requested = set((context or {}).get("terms", []) or [])
        if not required or not requested:
            return []
        return [
            term
            for term in required
            if any(normalized_text_has_term(ctx, term) for ctx in requested)
        ]

    def gold_candidate_match(self, gold: dict[str, Any], candidate: dict[str, Any]) -> dict[str, Any]:
        text = str(candidate.get("normalized_text") or "")
        axis_hits = [term for term in gold.get("mechanism_terms", []) or [] if normalized_text_has_term(text, term)]
        entity_hits = [term for term in gold.get("expected_entities", []) or [] if normalized_text_has_term(text, term)]
        trap_hits = [term for term in gold.get("negative_trap_terms", []) or [] if normalized_text_has_term(text, term)]
        score = 0
        if bool(gold.get("is_negative_control") or gold.get("polarity") == "negative"):
            if trap_hits:
                score += 3
        else:
            if axis_hits:
                score += 2
            if entity_hits:
                score += 1
        return {
            "candidate_id": candidate.get("candidate_id", ""),
            "source": candidate.get("source", ""),
            "label": candidate.get("label", ""),
            "axis_hits": axis_hits,
            "entity_hits": entity_hits,
            "trap_hits": trap_hits,
            "match_score": score,
            "has_evidence_refs": bool(candidate.get("has_evidence_refs")),
            "has_mechanism_trace": bool(candidate.get("has_mechanism_trace")),
            "has_input_trace": bool(candidate.get("has_input_trace")),
            "appendix": bool(candidate.get("appendix")),
            "allowed_claim_scope": str(candidate.get("allowed_claim_scope") or ""),
            "context_mismatch": bool(candidate.get("context_mismatch")),
        }

    def build_conclusion_evaluation(self, analysis_pack: dict[str, Any], context: dict[str, Any] | None = None) -> dict[str, Any]:
        gold_payload = self.load_gold_standard_conclusions()
        gold_entries = gold_payload.get("entries", []) or []
        candidates = self.conclusion_candidate_rows(analysis_pack)
        input_count = int((analysis_pack.get("input_summary") or {}).get("input_count") or 0)
        source_priority = {
            "mechanism_ready_fact": 0,
            "pathway_ranking": 1,
            "target_ranking": 2,
            "disease_ranking": 3,
            "input_exact_observation": 5,
            "class_observation": 6,
            "trait_observation": 8,
        }

        def candidate_rank_key(row: dict[str, Any]) -> tuple[int, int, str]:
            return (
                int(source_priority.get(str(row.get("source") or ""), 9)),
                int(row.get("rank") or 0),
                str(row.get("label") or ""),
            )

        if input_count > 50 and len(candidates) > 120:
            candidates = sorted(candidates, key=candidate_rank_key)[:120]

        def match_rank_key(row: dict[str, Any]) -> tuple[int, int, int, str, str]:
            return (
                -int(row.get("match_score") or 0),
                int(row.get("appendix") or 0),
                int(source_priority.get(str(row.get("source") or ""), 9)),
                str(row.get("source") or ""),
                str(row.get("label") or ""),
            )

        def class_observation_allowed_for_gold(gold: dict[str, Any]) -> bool:
            evidence_type = normalize_prediction_text(gold.get("required_evidence_type") or "")
            mechanism_axis = normalize_prediction_text(gold.get("mechanism_axis") or "")
            return any(
                token in evidence_type or token in mechanism_axis
                for token in ("class", "lipid", "fatty acid")
            )

        per_gold: list[dict[str, Any]] = []
        bug_queue: list[dict[str, Any]] = []
        positive_total = 0
        positive_matched = 0
        negative_total = 0
        negative_hits = 0
        matched_with_evidence = 0
        matched_total = 0
        context_compatible_matches = 0
        unsupported_top_claims = 0
        positive_gold_for_ranking: list[dict[str, Any]] = []
        for gold in gold_entries:
            context_status = self.gold_context_status(gold, context)
            matched_context_terms = self.gold_context_match_terms(gold, context)
            matches = [self.gold_candidate_match(gold, candidate) for candidate in candidates]
            matches = [match for match in matches if match["match_score"] > 0]
            matches.sort(key=match_rank_key)
            is_negative = bool(gold.get("is_negative_control") or gold.get("polarity") == "negative")
            positive_context_applicable = context_status != "context_mismatch"
            if is_negative:
                negative_total += 1
            elif gold.get("polarity") == "positive" and positive_context_applicable:
                positive_total += 1
                positive_gold_for_ranking.append(gold)
            status = "unsupported"
            if not is_negative and gold.get("polarity") == "positive" and not positive_context_applicable:
                status = "not_applicable_context"
            elif matches:
                top = matches[0]
                if is_negative and not top["appendix"]:
                    status = "false_positive_trap"
                    negative_hits += 1
                    bug_queue.append(
                        {
                            "code": "negative_trap_hit",
                            "severity": "error",
                            "gold_id": gold.get("gold_id", ""),
                            "matched_candidate": top,
                            "recommended_action": "Inspect schema, resolver, relation parser, literature parser, and ranking weights before release.",
                        }
                    )
                elif is_negative:
                    status = "negative_control_appendix_only"
                elif top.get("source") == "class_observation" and not class_observation_allowed_for_gold(gold):
                    status = "appendix_only"
                elif top["appendix"]:
                    status = "appendix_only"
                elif context_status == "context_mismatch" or top["context_mismatch"]:
                    status = "context_mismatch"
                elif top["axis_hits"] and top["entity_hits"] and top.get("source") != "input_exact_observation":
                    status = "exact_match"
                    positive_matched += 1
                else:
                    status = "partial_match"
                    positive_matched += 1
                if not is_negative and status in {"exact_match", "partial_match", "context_mismatch"}:
                    matched_total += 1
                    if top["has_evidence_refs"] or top["has_mechanism_trace"] or top.get("has_input_trace"):
                        matched_with_evidence += 1
                    if context_status in {"exact_context", "compatible_context", "no_context_requirement"} and not top["context_mismatch"]:
                        context_compatible_matches += 1
            elif is_negative:
                status = "negative_control_not_triggered"
            per_gold.append(
                {
                    "gold_id": gold.get("gold_id", ""),
                    "polarity": gold.get("polarity", "positive"),
                    "mechanism_axis": gold.get("mechanism_axis", ""),
                    "expected_entities": gold.get("expected_entities", []),
                    "context_status": context_status,
                    "matched_context_terms": matched_context_terms,
                    "gold_match_status": status,
                    "top_matches": matches[:20],
                    "source_refs": gold.get("source_refs", []),
                }
            )

        rankable_gold_statuses = {"exact_match", "partial_match", "context_mismatch"}
        context_priority = {
            "exact_context": 0,
            "compatible_context": 1,
            "no_context_requirement": 2,
            "context_not_provided": 3,
            "context_mismatch": 4,
        }
        status_priority = {
            "exact_match": 0,
            "partial_match": 1,
            "context_mismatch": 2,
        }

        def gold_conclusion_rank_key(row: dict[str, Any]) -> tuple[int, int, int, int, int, str]:
            top_match = (row.get("top_matches") or [{}])[0] if isinstance(row.get("top_matches"), list) else {}
            matched_context_terms = [term for term in row.get("matched_context_terms", []) or [] if term not in BROAD_CONTEXT_TERMS]
            evidence_trace = bool(
                top_match.get("has_evidence_refs")
                or top_match.get("has_mechanism_trace")
                or top_match.get("has_input_trace")
            )
            return (
                int(context_priority.get(str(row.get("context_status") or ""), 9)),
                -len(matched_context_terms),
                int(status_priority.get(str(row.get("gold_match_status") or ""), 9)),
                -int(top_match.get("match_score") or 0),
                0 if evidence_trace else 1,
                str(row.get("gold_id") or ""),
            )

        rankable_gold_conclusions = sorted(
            [
                row
                for row in per_gold
                if row.get("polarity") == "positive"
                and row.get("gold_match_status") in rankable_gold_statuses
            ],
            key=gold_conclusion_rank_key,
        )
        rankable_candidates = sorted([candidate for candidate in candidates if not candidate.get("appendix")], key=candidate_rank_key)
        top_candidates = rankable_candidates[:10]
        for candidate in top_candidates:
            if not candidate.get("has_evidence_refs") and not candidate.get("has_mechanism_trace") and not candidate.get("has_input_trace"):
                unsupported_top_claims += 1
        topk_metrics: dict[str, float | int | None] = {"ranked_candidate_count": len(rankable_candidates)}
        for k in (1, 3, 5):
            window = rankable_candidates[:k]
            candidate_hits = 0
            matched_gold_ids: set[str] = set()
            for candidate in window:
                candidate_matched = False
                for gold in positive_gold_for_ranking:
                    if candidate.get("source") == "class_observation" and not class_observation_allowed_for_gold(gold):
                        continue
                    match = self.gold_candidate_match(gold, candidate)
                    if int(match.get("match_score") or 0) <= 0:
                        continue
                    candidate_matched = True
                    matched_gold_ids.add(str(gold.get("gold_id") or ""))
                if candidate_matched:
                    candidate_hits += 1
            denominator = min(k, len(window))
            topk_metrics[f"precision_at_{k}"] = round(candidate_hits / denominator, 6) if denominator else None
            topk_metrics[f"recall_at_{k}"] = round(len(matched_gold_ids) / positive_total, 6) if positive_total else None
        metrics = {
            "gold_positive_count": positive_total,
            "gold_positive_matched": positive_matched,
            "gold_recall": round(positive_matched / positive_total, 6) if positive_total else None,
            "negative_trap_count": negative_total,
            "negative_trap_hits": negative_hits,
            "negative_trap_rate": round(negative_hits / negative_total, 6) if negative_total else None,
            "context_accuracy": round(context_compatible_matches / matched_total, 6) if matched_total else None,
            "evidence_usefulness": round(matched_with_evidence / matched_total, 6) if matched_total else None,
            "unsupported_top_claim_rate": round(unsupported_top_claims / len(top_candidates), 6) if top_candidates else None,
        }
        metrics.update(topk_metrics)
        return {
            "contract_version": CONCLUSION_EVALUATION_CONTRACT_VERSION,
            "gold_standard": {
                "entry_count": gold_payload.get("entry_count", 0),
                "files": gold_payload.get("files", []),
            },
            "candidate_count": len(candidates),
            "top_ranked_candidates": [
                {
                    "candidate_id": candidate.get("candidate_id", ""),
                    "source": candidate.get("source", ""),
                    "rank": candidate.get("rank", 0),
                    "label": candidate.get("label", ""),
                    "has_evidence_refs": bool(candidate.get("has_evidence_refs")),
                    "has_mechanism_trace": bool(candidate.get("has_mechanism_trace")),
                    "has_input_trace": bool(candidate.get("has_input_trace")),
                }
                for candidate in rankable_candidates[:10]
            ],
            "top_ranked_gold_conclusions": [
                {
                    "gold_id": row.get("gold_id", ""),
                    "status": row.get("gold_match_status", ""),
                    "mechanism_axis": row.get("mechanism_axis", ""),
                    "context_status": row.get("context_status", ""),
                    "matched_context_terms": row.get("matched_context_terms", []),
                    "top_match": (row.get("top_matches") or [{}])[0],
                }
                for row in rankable_gold_conclusions[:10]
            ],
            "per_gold": per_gold,
            "metrics": metrics,
            "bug_queue": bug_queue,
            "release_gate": {
                "passed": not bug_queue,
                "blocking_issue_count": len(bug_queue),
            },
        }

    def analysis_pack_edge_ref(self, edge: dict[str, Any]) -> dict[str, Any]:
        return {
            "ref_type": "edge",
            "edge_uid": edge.get("edge_uid", ""),
            "edge_type": edge.get("edge_type", ""),
            "subject_uid": edge.get("subject_uid", ""),
            "predicate": edge.get("predicate", ""),
            "object_uid": edge.get("object_uid", ""),
            "source_name": edge.get("source_name", ""),
            "source_record_id": edge.get("source_record_id", ""),
            "source_release": edge.get("source_release", ""),
            "license_id": edge.get("license_id", ""),
            "p_final": edge.get("p_final", 0.0),
        }

    def analysis_pack_sparse_edge_refs_for_node(self, node_uid: str, limit: int = 5) -> list[dict[str, Any]]:
        node = self.get_node(node_uid)
        if not node or node.get("node_idx") is None or not (self.graph_dir / "sparse_edges.parquet").exists():
            return []
        rows = self.sparse_incident_edges(int(node["node_idx"]))[:limit]
        node_idxs = sorted({idx for row in rows for idx in (row.get("source_idx"), row.get("target_idx")) if idx is not None})
        nodes_by_idx = {int(row["node_idx"]): row for row in self.get_nodes_by_idx(node_idxs)}
        refs = []
        for row in rows:
            source = nodes_by_idx.get(int(row.get("source_idx", -1)), {})
            target = nodes_by_idx.get(int(row.get("target_idx", -1)), {})
            refs.append(
                {
                    "ref_type": "edge",
                    "edge_uid": row.get("edge_uid", ""),
                    "edge_type": row.get("edge_type", ""),
                    "subject_uid": source.get("node_uid", ""),
                    "predicate": row.get("edge_type", ""),
                    "object_uid": target.get("node_uid", ""),
                    "source_name": "",
                    "source_record_id": "",
                    "source_release": self.release_id,
                    "license_id": "",
                    "p_final": row.get("weight", 0.0),
                }
            )
        return refs

    def analysis_pack_literature_ref(self, ref: dict[str, Any], edge_uid: str = "") -> dict[str, Any]:
        support_uid = ref.get("support_uid", "")
        supported_edges = parse_list(ref.get("supported_existing_edge_uids"))
        return {
            "ref_type": "literature_support",
            "edge_uid": edge_uid or ref.get("edge_uid", "") or (supported_edges[0] if supported_edges else ""),
            "support_uid": support_uid,
            "support_class": ref.get("support_class", ""),
            "p_literature": round(float(ref.get("p_literature") or 0.0), 6),
            "pmids": parse_list(ref.get("pmids"))[:10],
            "pmcids": parse_list(ref.get("pmcids"))[:10],
            "relation_uids": parse_list(ref.get("relation_uids") or ref.get("evidence_relation_uids"))[:10],
            "sentence_uids": parse_list(ref.get("sentence_uids"))[:10],
            "license_id": ref.get("license_id", ""),
            "source_release": ref.get("source_release", ""),
            "parser_hash": ref.get("parser_hash", ""),
            "config_hash": ref.get("config_hash", ""),
        }

    def analysis_pack_dedupe_refs(self, refs: list[dict[str, Any]], limit: int = 100) -> list[dict[str, Any]]:
        unique: dict[tuple[Any, ...], dict[str, Any]] = {}
        for ref in refs:
            if not ref:
                continue
            key = (
                ref.get("ref_type", ""),
                ref.get("edge_uid", ""),
                ref.get("support_uid", ""),
                ref.get("source_record_id", ""),
                tuple(ref.get("sentence_uids", []) or []),
            )
            unique.setdefault(key, ref)
        return sorted(
            unique.values(),
            key=lambda row: (
                0 if row.get("ref_type") == "literature_support" else 1,
                row.get("ref_type", ""),
                row.get("edge_uid", ""),
                row.get("support_uid", ""),
                row.get("source_record_id", ""),
            ),
        )[:limit]

    def analysis_pack_path_refs(self, path: dict[str, Any] | None) -> list[dict[str, Any]]:
        refs = []
        for edge in (path or {}).get("edges", []):
            refs.append(self.analysis_pack_edge_ref(edge))
            refs.extend(
                self.analysis_pack_literature_ref(support, edge_uid=edge.get("edge_uid", ""))
                for support in edge.get("literature_support", [])
            )
        return self.analysis_pack_dedupe_refs(refs)

    def analysis_pack_claim_refs(self, refs: list[dict[str, Any]], path_id: str = "") -> dict[str, Any]:
        edge_uids = sorted({ref.get("edge_uid", "") for ref in refs if ref.get("edge_uid")})
        evidence_ref_uids = sorted(
            {ref.get("support_uid", "") for ref in refs if ref.get("ref_type") == "literature_support" and ref.get("support_uid")}
        )
        source_records = sorted(
            {
                f"{ref.get('source_name', '')}:{ref.get('source_record_id', '')}"
                for ref in refs
                if ref.get("source_name") or ref.get("source_record_id")
            }
        )
        return {
            "path_id": path_id,
            "edge_uids": edge_uids,
            "evidence_ref_uids": evidence_ref_uids,
            "source_records": source_records[:20],
            "traceability_passed": bool(edge_uids or evidence_ref_uids or source_records),
        }

    def analysis_pack_ranking_rows(
        self,
        rows: list[dict[str, Any]],
        id_key: str,
        paths_by_id: dict[str, dict[str, Any]],
        features_by_uid: dict[str, list[dict[str, Any]]] | None = None,
        input_summary: dict[str, Any] | None = None,
        context: dict[str, Any] | None = None,
        limit: int = ANALYSIS_PACK_RANKING_LIMIT,
    ) -> list[dict[str, Any]]:
        ranking_rows = []
        features_by_uid = features_by_uid or {}
        input_summary = input_summary or {}
        for rank, row in enumerate(rows[:limit], start=1):
            path_id = row.get("best_path_id", "")
            path = paths_by_id.get(path_id)
            refs = []
            refs.extend(self.analysis_pack_path_refs(path))
            refs.extend(self.analysis_pack_edge_ref(ref) for ref in row.get("enrichment_evidence_refs", []) if ref.get("edge_uid"))
            refs.extend(self.analysis_pack_literature_ref(ref) for ref in row.get("evidence_refs", []) if ref.get("support_uid"))
            if not refs and row.get(id_key):
                refs.extend(self.analysis_pack_sparse_edge_refs_for_node(row.get(id_key, ""), limit=5))
            if not refs and row.get(id_key):
                refs.extend(self.analysis_pack_edge_ref(edge) for edge in self.incident_edges(row.get(id_key, ""), limit=5))
            refs = self.analysis_pack_dedupe_refs(refs)
            packed = {
                "rank": rank,
                id_key: row.get(id_key, ""),
                "display_name": row.get("name") or row.get("display_name", ""),
                "primary_external_id": row.get("primary_external_id", ""),
                "score": row.get("score", 0.0),
                "score_components": row.get("score_components", {}),
                "best_path_id": path_id,
                "best_path_cost": row.get("best_path_cost"),
                "path_confidence": row.get("path_confidence"),
                "literature_support_count": row.get("literature_support_count", 0),
                "max_p_literature": row.get("max_p_literature", 0.0),
                "supported_pmids": row.get("supported_pmids", []),
                "support_classes": row.get("support_classes", []),
                "claim_refs": self.analysis_pack_claim_refs(refs, path_id=path_id),
                "evidence_refs": refs,
            }
            if id_key == "pathway_uid":
                packed["overlap_count"] = row.get("overlap_count", 0)
                packed["matched_metabolite_uids"] = row.get("matched_metabolite_uids", [])
            packed.update(self.calibrated_prediction_for_ranking_row(packed, row, id_key, path, features_by_uid, input_summary, context))
            if packed.get("context_mismatch") and not packed.get("appendix"):
                packed["appendix"] = True
                packed["appendix_reason"] = packed.get("context_penalty_reason") or "context_mismatch_appendix_only"
            ranking_rows.append(packed)
        return ranking_rows

    def analysis_pack_path_row(self, path: dict[str, Any]) -> dict[str, Any]:
        nodes = path.get("nodes", [])
        terminal = nodes[-1] if nodes else {}
        refs = self.analysis_pack_path_refs(path)
        return {
            "path_id": path.get("path_id", ""),
            "cost": path.get("cost", 0.0),
            "path_confidence": path.get("path_confidence", 0.0),
            "terminal_node_uid": terminal.get("node_uid", ""),
            "terminal_node_type": terminal.get("node_type", ""),
            "terminal_display_name": terminal.get("display_name", ""),
            "node_uids": path.get("node_uids", []),
            "edge_uids": [edge.get("edge_uid", "") for edge in path.get("edges", [])],
            "claim_refs": self.analysis_pack_claim_refs(refs, path_id=path.get("path_id", "")),
            "evidence_refs": refs,
        }

    def select_analysis_pack_top_paths(
        self,
        paths: list[dict[str, Any]],
        ranking_sections: list[list[dict[str, Any]]],
        limit: int = ANALYSIS_PACK_PATH_LIMIT,
    ) -> list[dict[str, Any]]:
        paths_by_id = {path.get("path_id", ""): path for path in paths if path.get("path_id")}
        selected: list[dict[str, Any]] = []
        seen: set[str] = set()

        def add_path(path: dict[str, Any] | None) -> None:
            if not path or len(selected) >= limit:
                return
            path_id = str(path.get("path_id") or "")
            if not path_id or path_id in seen:
                return
            selected.append(path)
            seen.add(path_id)

        for rows in ranking_sections:
            for row in rows[:limit]:
                add_path(paths_by_id.get(str(row.get("best_path_id") or "")))
                if len(selected) >= limit:
                    return selected
        for path in paths:
            add_path(path)
            if len(selected) >= limit:
                break
        return selected

    def analysis_pack_fallback_explanation_chains(
        self,
        pathway_rankings: list[dict[str, Any]],
        target_rankings: list[dict[str, Any]],
        disease_rankings: list[dict[str, Any]],
        input_summary: dict[str, Any],
        limit: int = ANALYSIS_PACK_PATH_LIMIT,
    ) -> list[dict[str, Any]]:
        chains = []
        ranking_groups = (
            ("pathway_rankings", "pathway_uid", "pathway", pathway_rankings),
            ("target_rankings", "target_uid", "target", target_rankings),
            ("disease_rankings", "disease_uid", "disease", disease_rankings),
        )
        downgrade_reasons = [
            "no_direct_graph_path_under_current_hop_and_confidence_limits",
            "ranking_row_uses_propagation_or_literature_support_not_stepwise_causal_path",
        ]
        if int(input_summary.get("ambiguous_count", 0) or 0) or int(input_summary.get("expanded_candidate_count", 0) or 0):
            downgrade_reasons.append("ambiguous_or_expanded_seed_support_present")
        if str(input_summary.get("analysis_mode") or "") == "differential_table":
            downgrade_reasons.append("differential_table_input_is_precomputed_score_not_raw_abundance")

        for section, uid_key, terminal_type, rows in ranking_groups:
            for row in rows[:5]:
                if len(chains) >= limit:
                    return chains
                terminal_uid = str(row.get(uid_key) or "")
                if not terminal_uid:
                    continue
                refs = self.analysis_pack_dedupe_refs(row.get("evidence_refs", []) or [])
                claim_refs = row.get("claim_refs") or self.analysis_pack_claim_refs(refs)
                if not claim_refs.get("traceability_passed"):
                    continue
                seed_uids = [
                    str(uid)
                    for uid in row.get("matched_metabolite_uids", []) or []
                    if uid
                ][:10]
                if not seed_uids:
                    seed_uids = sorted(
                        {
                            str(ref.get("subject_uid") or "")
                            for ref in refs
                            if str(ref.get("subject_uid") or "").startswith("met_")
                        }
                        | {
                            str(ref.get("object_uid") or "")
                            for ref in refs
                            if str(ref.get("object_uid") or "").startswith("met_")
                        }
                    )[:10]
                confidence_source = (
                    row.get("calibrated_confidence")
                    if row.get("calibrated_confidence") is not None
                    else row.get("path_confidence")
                    if row.get("path_confidence") is not None
                    else row.get("max_p_literature")
                )
                confidence = min(0.25, clamp_unit(confidence_source if confidence_source is not None else 0.05))
                chain_id = content_hash(
                    {
                        "chain_type": "ranking_evidence_fallback",
                        "section": section,
                        "terminal_uid": terminal_uid,
                        "rank": row.get("rank", 0),
                        "evidence_ref_uids": claim_refs.get("evidence_ref_uids", []),
                        "edge_uids": claim_refs.get("edge_uids", []),
                    }
                )[:20]
                chains.append(
                    {
                        "chain_id": chain_id,
                        "path_id": chain_id,
                        "chain_type": "ranking_evidence_fallback",
                        "is_fallback_chain": True,
                        "source_ranking": section,
                        "source_rank": row.get("rank", 0),
                        "terminal_node_uid": terminal_uid,
                        "terminal_node_type": terminal_type,
                        "terminal_display_name": row.get("display_name", ""),
                        "node_uids": [*seed_uids[:5], terminal_uid] if seed_uids else [terminal_uid],
                        "edge_uids": claim_refs.get("edge_uids", []),
                        "path_confidence": round(confidence, 6),
                        "score": row.get("score", 0.0),
                        "score_components": row.get("score_components", {}),
                        "claim_refs": claim_refs,
                        "evidence_refs": refs,
                        "downgrade_reasons": downgrade_reasons,
                        "boundary": (
                            "Fallback explanation chain from ranking evidence. It is auditable support for a "
                            "research-priority row, not a stable stepwise graph path or causal conclusion."
                        ),
                    }
                )
        return chains

    def analysis_pack_literature_summary(self, refs: list[dict[str, Any]]) -> dict[str, Any]:
        literature_refs = [ref for ref in refs if ref.get("ref_type") == "literature_support"]
        p_values = [float(ref.get("p_literature") or 0.0) for ref in literature_refs]
        pmids = sorted({pmid for ref in literature_refs for pmid in ref.get("pmids", [])})
        support_classes = sorted({ref.get("support_class", "") for ref in literature_refs if ref.get("support_class")})
        return {
            "support_count": len({ref.get("support_uid", "") for ref in literature_refs if ref.get("support_uid")}),
            "max_p_literature": round(max(p_values), 6) if p_values else 0.0,
            "supported_pmids": pmids[:50],
            "support_classes": support_classes,
            "evidence_refs": literature_refs[:50],
        }

    def prediction_payload_from_ranking_row(self, row: dict[str, Any], uid_key: str) -> dict[str, Any]:
        return {
            "prediction_id": row.get("prediction_id", ""),
            "prediction": row.get("prediction", ""),
            "prediction_task": row.get("prediction_task", ""),
            "result_type": row.get("result_type", ""),
            "confidence_tier": row.get("confidence_tier", ""),
            "calibrated_confidence": row.get("calibrated_confidence", 0.0),
            "calibration_status": row.get("calibration_status", ""),
            "score": row.get("score", 0.0),
            "ranker_score": row.get("ranker_score", row.get("calibrated_confidence", 0.0)),
            "calibration_components": row.get("calibration_components", {}),
            "input_support_count": row.get("input_support_count", 0),
            "direction_consistency": row.get("direction_consistency", 0.0),
            "graph_distance": row.get("graph_distance"),
            "evidence_sources": row.get("evidence_sources", []),
            "is_extrapolated": bool(row.get("is_extrapolated")),
            "needs_validation": bool(row.get("needs_validation")),
            "boundary": row.get("boundary", ""),
            "display_name": row.get("display_name", ""),
            "result_uid": row.get(uid_key, ""),
            "primary_external_id": row.get("primary_external_id", ""),
            "claim_refs": row.get("claim_refs", {}),
            "coverage": row.get("coverage", {}),
            "theme_id": row.get("theme_id", ""),
            "context_mismatch": bool(row.get("context_mismatch")),
            "context_mismatch_groups": row.get("context_mismatch_groups", []),
            "downgrade_reason": row.get("downgrade_reason", ""),
            "appendix": bool(row.get("appendix")),
            "appendix_reason": row.get("appendix_reason", ""),
        }

    def prediction_model_assessment(
        self,
        precheck: dict[str, Any],
        predictions: list[dict[str, Any]],
        quality_warnings: list[dict[str, Any]] | None = None,
    ) -> dict[str, Any]:
        input_count = int(precheck.get("input_count") or 0)
        summary = precheck.get("summary") or {}
        matched = int(summary.get("matched") or 0)
        ambiguous = int(summary.get("ambiguous") or 0)
        unmatched = int(summary.get("unmatched") or 0) + int(summary.get("invalid") or 0)
        expanded_summary = precheck.get("expanded_summary") or {}
        expanded_inputs = int(expanded_summary.get("expanded_candidate_input_count") or expanded_summary.get("expanded_input_count") or 0)
        analysis_seed_inputs = int(expanded_summary.get("analysis_seed_input_count") or matched + expanded_inputs)
        denominator = max(1, input_count)
        matched_ratio = matched / denominator
        ambiguous_ratio = ambiguous / denominator
        unmatched_ratio = unmatched / denominator
        expanded_ratio = expanded_inputs / denominator
        analysis_seed_ratio = analysis_seed_inputs / denominator
        effective_input_ratio = (matched + 0.5 * expanded_inputs) / denominator
        tier_counts = Counter(row.get("confidence_tier", "low") for row in predictions)
        high_medium = tier_counts.get("high", 0) + tier_counts.get("medium", 0)
        warning_codes = {row.get("code", "") for row in quality_warnings or []}
        graph_lossy = "propagation_compression_lossy" in warning_codes or "propagation_graph_truncated" in warning_codes

        def graded_assessment(grade: str, labels: dict[str, tuple[str, str]]) -> dict[str, Any]:
            level, message = labels[grade]
            return {"grade": grade, "level": level, "message": message}

        theme_predictions = [row for row in predictions if row.get("result_type") in {"metabolic_theme", "metabolic_state"}]
        theme_tier_counts = Counter(row.get("confidence_tier", "low") for row in theme_predictions)
        has_medium_theme = theme_tier_counts.get("high", 0) + theme_tier_counts.get("medium", 0) > 0
        if not theme_predictions or (matched == 0 and expanded_inputs == 0):
            theme_grade = "D"
        elif matched_ratio >= 0.8 and ambiguous_ratio <= 0.1 and has_medium_theme:
            theme_grade = "A"
        elif matched_ratio >= 0.6 and has_medium_theme:
            theme_grade = "B"
        elif analysis_seed_ratio >= 0.6 and effective_input_ratio >= 0.45 and has_medium_theme:
            theme_grade = "C"
        else:
            theme_grade = "C"

        pathway_graph_grade = "B"
        if matched == 0 and expanded_inputs == 0:
            pathway_graph_grade = "D"
        elif matched_ratio < 0.4 and analysis_seed_ratio < 0.6:
            pathway_graph_grade = "D"
        elif graph_lossy:
            pathway_graph_grade = "C"
        elif matched_ratio >= 0.8 and ambiguous_ratio <= 0.1:
            pathway_graph_grade = "A"
        if matched == 0 and expanded_inputs == 0:
            grade = "D"
        elif matched_ratio < 0.4 and analysis_seed_ratio < 0.6:
            grade = "D"
        elif graph_lossy and matched_ratio >= 0.8 and ambiguous_ratio <= 0.1 and high_medium > 0:
            grade = "B"
        elif matched_ratio >= 0.8 and ambiguous_ratio <= 0.1 and unmatched_ratio <= 0.1 and high_medium > 0:
            grade = "A"
        elif matched_ratio >= 0.6 and ambiguous_ratio <= 0.25 and high_medium > 0:
            grade = "B"
        elif analysis_seed_ratio >= 0.6 and effective_input_ratio >= 0.45 and high_medium > 0:
            grade = "C"
        else:
            grade = "C"
        labels = {
            "A": ("well_calibrated", "Prediction quality is strong enough to prioritize high-confidence mechanisms."),
            "B": ("usable_with_caution", "Theme predictions are usable, but graph, context, or appendix outputs require caution."),
            "C": ("exploratory_usable", "Strict identities are limited, but expanded seeds provide usable research-level theme signals."),
            "D": ("input_limited", "Input support is too weak for reliable prediction; improve identifiers or keep results as a review queue."),
        }
        label, message = labels[grade]
        return {
            "grade": grade,
            "level": label,
            "message": message,
            "input_count": input_count,
            "matched_count": matched,
            "ambiguous_count": ambiguous,
            "expanded_candidate_input_count": expanded_inputs,
            "analysis_seed_input_count": analysis_seed_inputs,
            "unmatched_or_invalid_count": unmatched,
            "matched_ratio": round(matched_ratio, 6),
            "ambiguous_ratio": round(ambiguous_ratio, 6),
            "expanded_candidate_input_ratio": round(expanded_ratio, 6),
            "analysis_seed_input_ratio": round(analysis_seed_ratio, 6),
            "effective_input_support_ratio": round(effective_input_ratio, 6),
            "unmatched_or_invalid_ratio": round(unmatched_ratio, 6),
            "confidence_tier_counts": dict(sorted(tier_counts.items())),
            "high_or_medium_prediction_count": high_medium,
            "theme_calibration": graded_assessment(theme_grade, labels),
            "pathway_graph_calibration": graded_assessment(pathway_graph_grade, labels),
            "overall_calibration": graded_assessment(grade, labels),
            "graph_lossy": graph_lossy,
            "rules": [
                "high confidence requires enough matched input support, short graph distance, specific result type, and low genericity",
                "medium confidence can support candidate mechanisms but should be reviewed with biological context",
                "exploratory and low-confidence predictions remain useful leads but require validation",
            ],
        }

    def build_metabolic_theme_predictions(
        self,
        features_by_uid: dict[str, list[dict[str, Any]]],
        input_records: list[Any],
        input_summary: dict[str, Any],
        precheck: dict[str, Any] | None = None,
        context: dict[str, Any] | None = None,
    ) -> list[dict[str, Any]]:
        profile = self.input_theme_profile(features_by_uid, input_records=input_records)
        coverage_by_theme = self.build_theme_coverage(precheck or {}, input_summary)
        predictions = []
        theme_ids = sorted({theme_id for theme_id in profile if theme_id != "_lexical"} | set(coverage_by_theme))
        context_vector = (context or {}).get("context_vector") or {}
        has_context = bool((context or {}).get("has_context"))
        for theme_id in theme_ids:
            theme = profile.get(theme_id, {})
            coverage_row = coverage_by_theme.get(theme_id, {})
            if theme_id == "_lexical":
                continue
            coverage = coverage_row.get("coverage") or {}
            support_count = int(max(int(theme.get("support_count") or 0), int(coverage.get("related_input_count") or 0)))
            p_user_sum = float(max(float(theme.get("p_user_sum") or 0.0), float(coverage_row.get("p_user_sum") or 0.0)))
            if coverage_row:
                score = float(coverage_row.get("calibrated_confidence") or 0.0)
            else:
                raw_score = min(1.0, 0.25 + 0.18 * support_count + 0.08 * p_user_sum)
                score = round(raw_score * self.match_confidence_weight(input_summary), 6)
            context_boost = float(context_vector.get(theme_id) or 0.0)
            if context_boost and support_count > 0:
                score = round(min(0.95, score + context_boost), 6)
            coverage = coverage or {
                "related_input_count": support_count,
                "matched_support_count": support_count,
                "ambiguous_support_count": 0,
                "unmatched_support_count": 0,
                "direction_consistency": 0.0,
                "direction_summary": {},
                "biological_pattern_consistency": "not_pattern_evaluable",
                "significance_support_count": 0,
                "downgrade_reason": "",
            }
            score, tier, gate_reasons = self.apply_theme_confidence_gates(score, coverage)
            downgrade_reason = ";".join(
                reason for reason in [str(coverage.get("downgrade_reason") or ""), *gate_reasons] if reason
            )
            calibration_status = "calibrated_in_scope" if tier in {"high", "medium"} else "exploratory_in_scope"
            if downgrade_reason:
                calibration_status = "coverage_downgraded"
            input_terms = sorted(set(theme.get("input_terms", []) or []) | set(coverage_row.get("input_terms", []) or []))
            evidence_sources = ["input_theme", "biological_entity_pool"]
            if context_boost:
                evidence_sources.append("context_encoder")
            context_fit = "context_boosted" if context_boost else ("context_neutral" if has_context else "generalized")
            prediction = {
                "prediction_id": f"theme:{theme_id}",
                "prediction": f"The model predicts the {theme_id} metabolic theme from the input metabolite set.",
                "prediction_task": "metabolic_theme_prediction",
                "result_type": "metabolic_theme",
                "theme_id": theme_id,
                "display_name": theme_id.replace("_", " "),
                "confidence_tier": tier,
                "calibrated_confidence": score,
                "calibration_status": calibration_status,
                "input_support_count": support_count,
                "input_terms": input_terms[:20],
                "biological_entity_pool_priority": True,
                "coverage": {**coverage, "downgrade_reason": downgrade_reason},
                "context_fit": context_fit,
                "context_fit_score": round(context_boost, 6),
                "evidence_sources": evidence_sources,
                "graph_distance": 0,
                "is_extrapolated": False,
                "needs_validation": tier != "high",
                "downgrade_reason": downgrade_reason,
                "boundary": self.prediction_boundary_note(tier, 0, "metabolic", "calibrated_in_scope"),
            }
            predictions.append(prediction)
        for state_id, spec in METABOLIC_STATE_SIGNATURES.items():
            support = 0
            matched_terms = set()
            for features in features_by_uid.values():
                for feature in features:
                    labels = " ".join(
                        normalize_lookup_key(value)
                        for value in (feature.get("input_name", ""), feature.get("display_name", ""))
                    )
                    for term in spec["terms"]:
                        if normalized_text_has_term(labels, term):
                            support += 1
                            matched_terms.add(term)
            if support <= 0:
                continue
            score = round(min(0.85, 0.2 + 0.12 * support) * self.match_confidence_weight(input_summary), 6)
            context_boost = float((context_vector.get(state_id) or 0.0))
            if context_boost and support > 0:
                score = round(min(0.95, score + context_boost), 6)
            coverage = self.merge_theme_coverages(coverage_by_theme, spec.get("theme_ids", ()))
            if not coverage.get("related_input_count"):
                coverage = {
                    "related_input_count": support,
                    "matched_support_count": support,
                    "ambiguous_support_count": 0,
                    "unmatched_support_count": 0,
                    "direction_consistency": 0.0,
                    "direction_summary": {},
                    "biological_pattern_consistency": "not_pattern_evaluable",
                    "significance_support_count": 0,
                    "downgrade_reason": "state_signature_without_pool_coverage",
                }
            score, tier, gate_reasons = self.apply_theme_confidence_gates(score, coverage)
            downgrade_reason = ";".join(
                reason for reason in [str(coverage.get("downgrade_reason") or ""), *gate_reasons] if reason
            )
            prediction = {
                "prediction_id": f"state:{state_id}",
                "prediction": f"The model predicts {spec['label']} as a possible metabolic state.",
                "prediction_task": "metabolic_state_prediction",
                "result_type": "metabolic_state",
                "state_id": state_id,
                "display_name": spec["label"],
                "confidence_tier": tier,
                "calibrated_confidence": score,
                "calibration_status": "exploratory_in_scope" if tier == "exploratory" else "calibrated_in_scope",
                "input_support_count": support,
                "input_terms": sorted(matched_terms),
                "coverage": {**coverage, "downgrade_reason": downgrade_reason},
                "context_fit": "context_boosted" if context_boost else ("context_neutral" if has_context else "generalized"),
                "context_fit_score": round(context_boost, 6),
                "evidence_sources": ["input_state_signature"] + (["context_encoder"] if context_boost else []),
                "graph_distance": 0,
                "is_extrapolated": False,
                "needs_validation": tier != "high",
                "downgrade_reason": downgrade_reason,
                "boundary": self.prediction_boundary_note(tier, 0, "metabolic", "calibrated_in_scope"),
            }
            if downgrade_reason:
                prediction["calibration_status"] = "coverage_downgraded"
            predictions.append(prediction)
        predictions.sort(key=lambda row: (-float(row.get("calibrated_confidence") or 0.0), row.get("prediction_id", "")))
        return predictions[:20]

    def build_prediction_model_pack(
        self,
        precheck: dict[str, Any],
        features_by_uid: dict[str, list[dict[str, Any]]],
        input_records: list[Any],
        input_summary: dict[str, Any],
        pathway_rankings: list[dict[str, Any]],
        target_rankings: list[dict[str, Any]],
        disease_rankings: list[dict[str, Any]],
        quality_warnings: list[dict[str, Any]] | None = None,
        context: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        theme_predictions = self.build_metabolic_theme_predictions(features_by_uid, input_records, input_summary, precheck, context)
        theme_tiers = [row.get("confidence_tier", "low") for row in theme_predictions]
        tier_rank = {"low": 0, "exploratory": 1, "medium": 2, "high": 3}
        best_theme_tier = max(theme_tiers, key=lambda tier: tier_rank.get(tier, 0), default="low")
        generalized_mode = not bool((context or {}).get("has_context"))
        large_input = int(input_summary.get("input_count") or 0) > 50
        if large_input:
            pathway_rankings = pathway_rankings[:10]
            target_rankings = target_rankings[:10]
            disease_rankings = disease_rankings[:10]
        adjusted_disease_rankings = []
        for row in disease_rankings:
            adjusted = dict(row)
            if tier_rank.get(best_theme_tier, 0) <= tier_rank["exploratory"] and tier_rank.get(adjusted.get("confidence_tier", "low"), 0) > tier_rank["exploratory"]:
                adjusted["upstream_confidence_tier"] = best_theme_tier
                adjusted["calibrated_confidence"] = min(float(adjusted.get("calibrated_confidence") or 0.0), 0.25)
                adjusted["confidence_tier"] = "exploratory"
                adjusted["calibration_status"] = "upstream_confidence_limited"
                adjusted["downgrade_reason"] = "disease_association_inherits_exploratory_metabolic_theme"
                adjusted["needs_validation"] = True
                adjusted["appendix"] = True
                adjusted["appendix_reason"] = "disease_association_inherits_exploratory_metabolic_theme"
            if generalized_mode:
                adjusted["calibrated_confidence"] = min(float(adjusted.get("calibrated_confidence") or 0.0), 0.25)
                adjusted["confidence_tier"] = self.confidence_tier(float(adjusted["calibrated_confidence"]))
                adjusted["calibration_status"] = "appendix_low"
                adjusted["downgrade_reason"] = ";".join(
                    reason
                    for reason in [str(adjusted.get("downgrade_reason") or ""), "generalized_mode_disease_prediction_appendix"]
                    if reason
                )
                adjusted["needs_validation"] = True
                adjusted["appendix"] = True
                adjusted["appendix_reason"] = "generalized_mode_disease_prediction_appendix"
            adjusted_disease_rankings.append(adjusted)
        predictions = [
            *theme_predictions,
            *(self.prediction_payload_from_ranking_row(row, "pathway_uid") for row in pathway_rankings),
            *(self.prediction_payload_from_ranking_row(row, "target_uid") for row in target_rankings),
            *(self.prediction_payload_from_ranking_row(row, "disease_uid") for row in adjusted_disease_rankings),
        ]
        predictions.sort(
            key=lambda row: (
                {"high": 0, "medium": 1, "exploratory": 2, "low": 3}.get(row.get("confidence_tier", ""), 4),
                -float(row.get("calibrated_confidence") or 0.0),
                row.get("prediction_id", ""),
            )
        )
        context_mismatch = [row for row in predictions if row.get("context_mismatch")]
        appendix = [
            row
            for row in predictions
            if row.get("appendix") or row.get("confidence_tier") == "low"
            if row.get("confidence_tier") != "high"
        ]
        degraded_with_coverage = [
            row
            for row in theme_predictions
            if row.get("coverage", {}).get("related_input_count", 0)
            and row.get("confidence_tier") in {"exploratory", "low"}
        ]
        return {
            "contract_version": "metabolic_prediction_model.v1",
            "model_name": "precision_metabolic_predictor_rule_mvp",
            "assessment": self.prediction_model_assessment(precheck, predictions, quality_warnings),
            "predictions": predictions[:100],
            "metabolic_themes": theme_predictions[:30],
            "high_confidence_metabolic_themes": [row for row in theme_predictions if row.get("confidence_tier") == "high"][:10],
            "medium_confidence_metabolic_themes": [row for row in theme_predictions if row.get("confidence_tier") == "medium"][:15],
            "degraded_covered_metabolic_themes": degraded_with_coverage[:20],
            "context_mismatch": context_mismatch[:20],
            "low_confidence_appendix": appendix[:30],
            "high_confidence": [row for row in predictions if row.get("confidence_tier") == "high"][:10],
            "medium_confidence": [row for row in predictions if row.get("confidence_tier") == "medium"][:15],
            "exploratory": [row for row in predictions if row.get("confidence_tier") == "exploratory"][:20],
            "low_confidence": [row for row in predictions if row.get("confidence_tier") == "low"][:20],
            "tasks": {
                "entity_precision": "handled by resolver/precheck rows with candidate scores and ambiguity detection",
                "metabolic_theme_prediction": "rule MVP from metabolite set themes and state signatures",
                "pathway_prediction": "calibrated reranking over enrichment, graph, direction, evidence, type, distance, and genericity",
                "target_prediction": "calibrated graph/literature candidate target ranking",
                "phenotype_context_prediction": "calibrated disease or phenotype-context association ranking",
            },
        }

    def structured_theme_row(self, row: dict[str, Any]) -> dict[str, Any]:
        coverage = row.get("coverage") or {}
        matched_support_count = int(coverage.get("matched_support_count") or 0)
        related_input_count = max(int(coverage.get("related_input_count") or row.get("input_support_count") or 0), matched_support_count)
        return {
            "theme_id": row.get("theme_id") or row.get("state_id") or "",
            "display_name": row.get("display_name", ""),
            "confidence": row.get("confidence_tier", "low"),
            "confidence_tier": row.get("confidence_tier", "low"),
            "related_input_count": related_input_count,
            "matched_support_count": matched_support_count,
            "ambiguous_support_count": int(coverage.get("ambiguous_support_count") or 0),
            "unmatched_support_count": int(coverage.get("unmatched_support_count") or 0),
            "significant_support_count": int(coverage.get("significance_support_count") or row.get("significant_support_count") or 0),
            "direction_consistency": coverage.get("direction_consistency", row.get("direction_consistency", 0.0)),
            "biological_pattern_consistency": coverage.get("biological_pattern_consistency", "not_pattern_evaluable"),
            "context_fit": row.get("context_fit", "generalized"),
            "context_fit_score": row.get("context_fit_score", 0.0),
            "evidence_sources": row.get("evidence_sources", []),
            "evidence_refs": row.get("evidence_refs", []),
            "claim_refs": row.get("claim_refs", {}),
            "downgrade_reason": coverage.get("downgrade_reason") or row.get("downgrade_reason", ""),
            "appendix": bool(row.get("appendix")),
            "appendix_reason": row.get("appendix_reason", ""),
            "research_only": bool(row.get("research_only")),
            "review_required": bool(row.get("needs_validation") or row.get("confidence_tier") in {"exploratory", "low"}),
        }

    def pathway_linked_theme(self, row: dict[str, Any]) -> str:
        directional = row.get("directional_support") or {}
        for field in ("input_theme_hits", "pathway_theme_hits"):
            for hit in directional.get(field, []) or []:
                if isinstance(hit, dict) and hit.get("theme_id"):
                    return str(hit["theme_id"])
        return ""

    def structured_prediction_json(
        self,
        analysis_pack: dict[str, Any],
        prediction_pack: dict[str, Any],
        context: dict[str, Any],
    ) -> dict[str, Any]:
        input_summary = analysis_pack.get("input_summary") or {}
        model = analysis_pack.get("prediction_model") or {}
        assessment = model.get("assessment") or {}

        def calibration_label(key: str) -> str:
            row = assessment.get(key) or {}
            if not row:
                return ""
            return " ".join(part for part in (row.get("grade", ""), row.get("level", "")) if part)

        themes = [self.structured_theme_row(row) for row in model.get("metabolic_themes", []) or []]
        high = [row for row in themes if row["confidence"] == "high"]
        medium = [row for row in themes if row["confidence"] == "medium"]
        downgraded_ids = {
            row.get("theme_id") or row.get("state_id") or ""
            for row in model.get("degraded_covered_metabolic_themes", []) or []
        }
        exploratory = [row for row in themes if row["confidence"] in {"exploratory", "low"} and row["theme_id"] not in downgraded_ids]
        downgraded = [row for row in themes if row["theme_id"] in downgraded_ids]

        pathway_evidence = []
        for row in analysis_pack.get("pathway_rankings", []) or []:
            pathway_evidence.append(
                {
                    "pathway_id": row.get("pathway_uid", ""),
                    "pathway_name": row.get("display_name") or row.get("name", ""),
                    "source_database": row.get("source_name", ""),
                    "linked_theme": self.pathway_linked_theme(row),
                    "input_support": row.get("input_support_count", 0),
                    "graph_distance": row.get("graph_distance", 0),
                    "raw_graph_score": (row.get("calibration_components") or {}).get("raw_score", row.get("score", 0.0)),
                    "adjusted_score": row.get("calibrated_confidence", row.get("score", 0.0)),
                    "confidence": row.get("confidence_tier", "low"),
                    "confidence_tier": row.get("confidence_tier", "low"),
                    "calibration": row.get("calibration_status", ""),
                    "context_mismatch": bool(row.get("context_mismatch")),
                    "appendix": bool(row.get("appendix")),
                    "appendix_reason": row.get("appendix_reason", ""),
                    "research_only": bool(row.get("research_only") or row.get("appendix")),
                    "evidence_refs": row.get("evidence_refs", []),
                    "claim_refs": row.get("claim_refs", {}),
                }
            )

        def excluded_from_primary_reason(row: dict[str, Any], result_type: str) -> str:
            tier = str(row.get("confidence_tier") or row.get("display_confidence") or "low")
            if row.get("context_mismatch"):
                return "context_mismatch"
            if row.get("research_only"):
                return "research_only"
            if row.get("appendix"):
                return str(row.get("appendix_reason") or "appendix")
            if result_type in {"drug", "disease"} and (context.get("has_context") is False or not context.get("has_context")):
                return "generalized_mode_research_only"
            if tier in {"exploratory", "low"}:
                return "low_or_exploratory_confidence"
            return ""

        def prediction_brief(row: dict[str, Any], uid_key: str, result_type: str) -> dict[str, Any]:
            reason = excluded_from_primary_reason(row, result_type)
            return {
                "prediction_id": row.get("prediction_id") or row.get(uid_key, "") or row.get("result_uid", ""),
                "result_id": row.get(uid_key, "") or row.get("result_uid", ""),
                "result_type": result_type,
                "display_name": row.get("display_name", ""),
                "raw_overlay_score": (row.get("calibration_components") or {}).get("raw_score", row.get("score", 0.0)),
                "display_confidence": row.get("confidence_tier", "low"),
                "confidence": row.get("confidence_tier", "low"),
                "confidence_tier": row.get("confidence_tier", "low"),
                "calibration": row.get("calibration_status", ""),
                "evidence_sources": row.get("evidence_sources", []),
                "evidence_refs": row.get("evidence_refs", []),
                "claim_refs": row.get("claim_refs", {}),
                "downgrade_reason": row.get("downgrade_reason", ""),
                "appendix": bool(row.get("appendix")),
                "appendix_reason": row.get("appendix_reason", ""),
                "research_only": bool(row.get("research_only")),
                "excluded_from_primary_reason": reason,
                "review_flag": bool(row.get("needs_validation") or row.get("appendix")),
            }

        prediction_rows = model.get("predictions", []) or []
        target_predictions = [prediction_brief(row, "target_uid", "target") for row in prediction_rows if row.get("result_type") == "target"]
        disease_predictions = [prediction_brief(row, "disease_uid", "disease") for row in prediction_rows if row.get("result_type") == "disease"]
        drug_hypotheses = []
        for row in prediction_pack.get("drug_rankings", []) or []:
            reason = excluded_from_primary_reason({**row, "research_only": True, "appendix": True}, "drug")
            drug_hypotheses.append(
                {
                    "prediction_id": row.get("prediction_id") or row.get("drug_id") or row.get("drug_name", ""),
                    "result_type": "drug",
                    "drug_name": row.get("drug_name", ""),
                    "raw_overlay_score": row.get("raw_overlay_score", 0.0),
                    "display_confidence": row.get("display_confidence", row.get("confidence_tier", "low")),
                    "confidence_tier": row.get("confidence_tier", "low"),
                    "upstream_targets": row.get("matched_targets", []),
                    "upstream_target_confidence": row.get("upstream_confidence_tier", "low"),
                    "downgrade_reason": row.get("downgrade_reason", ""),
                    "clinical_warning": "Research-only hypothesis; not a treatment recommendation.",
                    "research_only": True,
                    "appendix": True,
                    "excluded_from_primary_reason": reason or "clinical_warning_or_research_only",
                    "evidence_refs": row.get("evidence_refs", []),
                    "claim_refs": row.get("claim_refs", {}),
                }
            )

        review_queue = []
        for bucket in ("ambiguous", "unmatched"):
            for row in analysis_pack.get(bucket, []) or []:
                review_queue.append({"type": f"entity_{bucket}", "raw_input": row.get("raw_input") or row.get("query", ""), "reason": row.get("resolution_status", bucket)})
        for row in [*target_predictions, *disease_predictions, *drug_hypotheses]:
            if row.get("review_flag") or row.get("appendix") or row.get("research_only"):
                review_queue.append({"type": "prediction_review", "display_name": row.get("display_name") or row.get("drug_name", ""), "reason": row.get("downgrade_reason", "") or "appendix_or_research_only"})

        def layer_item(row: dict[str, Any], result_type: str) -> dict[str, Any]:
            return {
                "prediction_id": row.get("prediction_id") or row.get("result_id") or row.get("pathway_id") or row.get("theme_id") or row.get("drug_name", ""),
                "result_type": result_type,
                "display_name": row.get("display_name") or row.get("drug_name") or row.get("pathway_name", ""),
                "confidence_tier": row.get("confidence_tier") or row.get("confidence") or row.get("display_confidence", "low"),
                "appendix": bool(row.get("appendix")),
                "research_only": bool(row.get("research_only")),
                "excluded_from_primary_reason": row.get("excluded_from_primary_reason", ""),
                "evidence_refs": row.get("evidence_refs", []),
                "claim_refs": row.get("claim_refs", {}),
            }

        theme_items = [layer_item(row, "metabolic_theme") for row in [*high, *medium, *exploratory, *downgraded]]
        pathway_items = [
            layer_item({**row, "prediction_id": row.get("pathway_id"), "display_name": row.get("pathway_name")}, "pathway")
            for row in pathway_evidence
        ]
        target_items = [layer_item(row, "target") for row in target_predictions]
        disease_items = [layer_item(row, "disease") for row in disease_predictions]
        drug_items = [layer_item(row, "drug") for row in drug_hypotheses]
        all_prediction_items = [*theme_items, *pathway_items, *target_items, *disease_items, *drug_items]
        primary_research_candidates = [
            row
            for row in all_prediction_items
            if row["result_type"] in {"metabolic_theme", "pathway", "target"}
            and row["confidence_tier"] in {"high", "medium"}
            and not row["appendix"]
            and not row["research_only"]
            and not row["excluded_from_primary_reason"]
        ]
        mechanistic_support = [
            row
            for row in [*theme_items, *pathway_items, *target_items]
            if row["evidence_refs"] or row["claim_refs"]
        ]
        appendix_overlay_only = [
            row
            for row in all_prediction_items
            if row["appendix"] and row["result_type"] in {"target", "disease", "drug", "pathway"}
        ]
        clinical_warning_or_research_only = [
            row for row in all_prediction_items if row["research_only"] or row["result_type"] == "drug"
        ]
        excluded_from_primary = [
            row for row in all_prediction_items if row["excluded_from_primary_reason"] or row["appendix"] or row["research_only"]
        ]
        product_layers = {
            "primary_research_candidates": primary_research_candidates[:30],
            "mechanistic_support": mechanistic_support[:50],
            "appendix_overlay_only": appendix_overlay_only[:50],
            "context_mismatch": [layer_item(row, str(row.get("result_type") or "prediction")) for row in model.get("context_mismatch", []) or []][:30],
            "clinical_warning_or_research_only": clinical_warning_or_research_only[:50],
            "excluded_from_primary_reason": excluded_from_primary[:100],
        }

        return {
            "mode": "context_aware" if context.get("has_context") else "generalized",
            "context": context,
            "input_quality": {
                "total_inputs": int(input_summary.get("input_count") or 0),
                "matched": int(input_summary.get("matched_count") or 0),
                "ambiguous": int(input_summary.get("ambiguous_count") or 0),
                "unmatched": int(input_summary.get("unmatched_count") or 0) + int(input_summary.get("invalid_count") or 0),
                "genetic_exposures": int(input_summary.get("genetic_exposure_count") or 0),
                "biological_entity_pools": len(analysis_pack.get("biological_entity_pools") or []),
            },
            "calibration": {
                "theme": calibration_label("theme_calibration"),
                "pathway_graph": calibration_label("pathway_graph_calibration"),
                "target_drug_disease": calibration_label("overall_calibration"),
                "overall": calibration_label("overall_calibration"),
            },
            "high_confidence_themes": high,
            "medium_confidence_themes": medium,
            "exploratory_themes": exploratory,
            "downgraded_but_supported_themes": downgraded,
            "pathway_evidence": pathway_evidence,
            "target_predictions": target_predictions,
            "disease_predictions": disease_predictions,
            "drug_hypotheses": drug_hypotheses,
            "context_mismatch_results": model.get("context_mismatch", []) or [],
            "product_layers": product_layers,
            "review_queue": review_queue[:100],
            "warnings": (analysis_pack.get("quality_warnings") or []) + (prediction_pack.get("blocked_reasons") or []),
            "evidence_refs": analysis_pack.get("evidence_refs", []) or [],
        }

    def interpretation_input_label(self, row: dict[str, Any]) -> str:
        record = row.get("record") or {}
        for key in ("name", "metabolite", "compound", "HMDB", "ChEBI", "PubChem CID", "InChIKey", "query"):
            value = record.get(key) if isinstance(record, dict) else None
            if value:
                return str(value)
        return str(row.get("query") or row.get("raw_input") or row.get("input_id") or "input")

    def interpretation_candidate_label(self, row: dict[str, Any]) -> str:
        candidates = row.get("candidates") or row.get("candidate_entities") or []
        if candidates:
            candidate = candidates[0] or {}
            return str(candidate.get("display_name") or candidate.get("canonical_name") or candidate.get("entity_uid") or "")
        return str(row.get("display_name") or row.get("metabolite_uid") or "")

    def interpretation_input_support_index(self, analysis_pack: dict[str, Any]) -> dict[str, list[dict[str, Any]]]:
        by_uid: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for row in analysis_pack.get("matched", []) or []:
            uid = str(row.get("metabolite_uid") or row.get("normalized_entity_id") or "")
            if uid:
                by_uid[uid].append(row)
        return by_uid

    def interpretation_supporting_inputs(
        self,
        metabolite_uids: list[str],
        support_index: dict[str, list[dict[str, Any]]],
        limit: int = 8,
    ) -> list[dict[str, Any]]:
        inputs = []
        for uid in metabolite_uids:
            for row in support_index.get(uid, []):
                record = row.get("record") or {}
                effect_source = str(record.get("comparison_effect_source") or "")
                effect_label = str(record.get("effect_label") or (comparison_effect_display_label(effect_source) if effect_source else "log2FC"))
                effect_value = record.get("effect_value", record.get("log2FC", record.get("log2fc", "")))
                inputs.append(
                    {
                        "input_id": row.get("input_id", ""),
                        "input_name": self.interpretation_input_label(row),
                        "matched_name": self.interpretation_candidate_label(row),
                        "metabolite_uid": uid,
                        "direction": record.get("direction", ""),
                        "log2FC": record.get("log2FC", record.get("log2fc", "")),
                        "effect_label": effect_label,
                        "effect_value": effect_value,
                        "effect_is_surrogate": bool(record.get("log2FC_semantics") or record.get("log2fc_semantics")),
                        "pvalue": record.get("pvalue", record.get("p", "")),
                        "padj": record.get("padj", record.get("qvalue", "")),
                        "match_score": row.get("match_score", row.get("top_score", 0.0)),
                        "match_status": row.get("match_status", row.get("resolution_status", "")),
                    }
                )
        return inputs[:limit]

    def interpretation_confidence_reasons(
        self,
        row: dict[str, Any],
        support_inputs: list[dict[str, Any]],
        input_summary: dict[str, Any],
        source_kind: str,
    ) -> dict[str, list[str]]:
        positive: list[str] = []
        downgrade: list[str] = []
        evidence_refs = row.get("evidence_refs") or []
        claim_refs = row.get("claim_refs") or {}
        support_count = int(row.get("supporting_input_count") or row.get("input_support_count") or row.get("overlap_count") or len(support_inputs) or 0)
        if support_count > 1:
            positive.append(f"{support_count} 个已匹配输入共同支持该条结论")
        elif support_count == 1:
            positive.append("1 个已匹配输入支持该条结论")
        if claim_refs.get("traceability_passed"):
            positive.append("结论可回链到图谱边、文献支持或来源记录")
        if evidence_refs:
            positive.append(f"包含 {len(evidence_refs)} 个可追溯证据引用")
        if row.get("max_p_literature"):
            positive.append("存在句级文献证据支持")
        if row.get("confidence_tier") in {"high", "medium"}:
            positive.append(f"校准层级为 {row.get('confidence_tier')}")
        support_classes = {str(ref.get("support_class") or "") for ref in evidence_refs}
        if int(input_summary.get("ambiguous_count") or 0):
            downgrade.append("部分输入存在歧义，不能支撑精确化学身份的强结论")
        if int(input_summary.get("unmatched_count") or 0) or int(input_summary.get("invalid_count") or 0):
            downgrade.append("部分输入未能匹配到当前 release")
        if row.get("context_mismatch"):
            downgrade.append("与用户给定背景存在不匹配，已降级")
        if row.get("appendix"):
            downgrade.append(str(row.get("appendix_reason") or "仅作为附录结果"))
        if row.get("downgrade_reason"):
            downgrade.append(str(row.get("downgrade_reason")))
        if support_classes & {"novel_candidate", "conflict_candidate"}:
            downgrade.append("部分文献证据属于 novel/conflict overlay，只能作为探索性线索")
        if source_kind in {"target", "disease", "drug"}:
            downgrade.append("靶点、疾病和药物条目属于研究假设，不是直接因果结论")
        if not evidence_refs and not claim_refs.get("traceability_passed"):
            downgrade.append("可追溯证据较弱或缺失")
        return {
            "positive_factors": positive[:6],
            "downgrade_factors": list(dict.fromkeys(downgrade))[:6],
        }

    def interpretation_boundary(self, row: dict[str, Any], claim_type: str) -> str:
        if claim_type in {"metabolic_theme", "pathway"}:
            return "支持代谢主题或通路层面的研究信号；不能单独证明样本中通路活性增强、减弱或存在因果关系。"
        if claim_type == "target":
            return "可作为由图谱或文献关联得到的机制候选；不等同于证明该靶点被输入代谢物直接调控。"
        if claim_type == "disease":
            return "可作为疾病或表型背景关联阅读；不是诊断、预测或临床判断。"
        if claim_type == "drug":
            return "仅作为研究用药物或扰动假设；不是治疗建议。"
        if row.get("boundary"):
            return str(row["boundary"])
        return "仅用于研究解释；在作为生物学事实使用前需要独立验证。"

    def interpretation_next_validation(self, claim_type: str) -> list[str]:
        if claim_type in {"metabolic_theme", "pathway"}:
            return [
                "检查原始数据中同一通路的更多代谢物是否呈一致变化。",
                "结合关键酶的转录、蛋白或活性数据进行验证。",
                "优先用靶向 LC-MS/MS 复核代表性代谢物。",
            ]
        if claim_type == "target":
            return [
                "检查同一生物背景下该靶点的基因或蛋白表达。",
                "查找扰动、knockdown、抑制剂或依赖性数据。",
                "确认该靶点关系是否为短路径且有证据支持。",
            ]
        return [
            "先作为文献检索或候选排序线索使用。",
            "在独立队列、实验或背景匹配的公开数据中验证。",
        ]

    def interpretation_theme_claims(
        self,
        analysis_pack: dict[str, Any],
        support_index: dict[str, list[dict[str, Any]]],
        input_summary: dict[str, Any],
    ) -> list[dict[str, Any]]:
        model = analysis_pack.get("prediction_model") or {}
        rows = [
            *(model.get("high_confidence_metabolic_themes") or []),
            *(model.get("medium_confidence_metabolic_themes") or []),
            *(model.get("degraded_covered_metabolic_themes") or []),
        ]
        claims = []
        for row in rows:
            coverage = row.get("coverage") or {}
            support_uids = list(coverage.get("matched_metabolite_uids") or row.get("matched_metabolite_uids") or [])
            support_inputs = self.interpretation_supporting_inputs(support_uids, support_index)
            label = str(row.get("display_name") or row.get("theme_id") or "metabolic theme")
            claim_type = "metabolic_theme"
            confidence_tier = str(row.get("confidence_tier") or row.get("confidence") or "exploratory")
            claim_id = f"claim_theme_{len(claims) + 1:03d}"
            support_count = int(
                coverage.get("related_input_count")
                or row.get("input_support_count")
                or row.get("overlap_count")
                or len(support_inputs)
                or 0
            )
            row_for_confidence = {**row, "supporting_input_count": support_count}
            claims.append(
                {
                    "claim_id": claim_id,
                    "claim_type": claim_type,
                    "headline": f"输入代谢物指向 {label} 代谢主题。",
                    "basis": f"{coverage.get('related_input_count', row.get('input_support_count', 0))} 条输入相关，其中 {coverage.get('matched_support_count', 0)} 条为严格匹配支持。",
                    "relation_chain": ["输入代谢物", label, "代谢主题解释"],
                    "possible_mechanism": f"提示可能存在 {label} 相关研究信号，但在写成机制前需要正交验证。",
                    "supporting_inputs": support_inputs,
                    "supporting_input_count": support_count,
                    "supporting_nodes": [{"node_uid": row.get("theme_id", ""), "display_name": label, "node_type": "metabolic_theme"}],
                    "supporting_edges": [],
                    "evidence_refs": row.get("evidence_refs", []),
                    "confidence_tier": confidence_tier,
                    "calibrated_confidence": row.get("calibrated_confidence", row.get("confidence_score", 0.0)),
                    "confidence_reasons": self.interpretation_confidence_reasons(row_for_confidence, support_inputs, input_summary, claim_type),
                    "boundary": self.interpretation_boundary(row, claim_type),
                    "related_candidates": [],
                    "next_validation": self.interpretation_next_validation(claim_type),
                    "review_required": confidence_tier in {"exploratory", "low"} or bool(row.get("needs_validation")),
                    "source": "prediction_model.metabolic_themes",
                }
            )
        return claims

    def interpretation_pathway_claims(
        self,
        analysis_pack: dict[str, Any],
        support_index: dict[str, list[dict[str, Any]]],
        input_summary: dict[str, Any],
        existing_count: int,
    ) -> list[dict[str, Any]]:
        claims = []
        for row in (analysis_pack.get("pathway_rankings") or [])[:5]:
            support_uids = list(row.get("matched_metabolite_uids") or [])
            support_inputs = self.interpretation_supporting_inputs(support_uids, support_index)
            if existing_count and len(support_inputs) == 0 and int(row.get("overlap_count") or 0) == 0:
                continue
            label = str(row.get("display_name") or row.get("pathway_uid") or "pathway")
            claim_type = "pathway"
            claim_id = f"claim_pathway_{len(claims) + 1:03d}"
            support_count = int(row.get("input_support_count") or row.get("overlap_count") or len(support_inputs) or 0)
            row_for_confidence = {**row, "supporting_input_count": support_count}
            related_targets = [
                item.get("display_name") or item.get("target_uid")
                for item in (analysis_pack.get("target_rankings") or [])[:3]
                if item.get("display_name") or item.get("target_uid")
            ]
            claims.append(
                {
                    "claim_id": claim_id,
                    "claim_type": claim_type,
                    "headline": f"输入代谢物连接到 {label} 通路。",
                    "basis": f"{len(support_inputs) or row.get('overlap_count', 0)} 个已匹配输入连接到该通路；当前排名为 {row.get('rank', '')}。",
                    "relation_chain": ["输入代谢物", label, "通路排序", *(related_targets[:2] or [])],
                    "possible_mechanism": f"{label} 是对上传代谢物集合的通路层解释，不是直接通路活性测定。",
                    "supporting_inputs": support_inputs,
                    "supporting_input_count": support_count,
                    "supporting_nodes": [{"node_uid": row.get("pathway_uid", ""), "display_name": label, "node_type": "pathway"}],
                    "supporting_edges": (row.get("claim_refs") or {}).get("edge_uids", []),
                    "evidence_refs": row.get("evidence_refs", []),
                    "confidence_tier": row.get("confidence_tier", "exploratory"),
                    "calibrated_confidence": row.get("calibrated_confidence", row.get("score", 0.0)),
                    "confidence_reasons": self.interpretation_confidence_reasons(row_for_confidence, support_inputs, input_summary, claim_type),
                    "boundary": self.interpretation_boundary(row, claim_type),
                    "related_candidates": related_targets[:5],
                    "next_validation": self.interpretation_next_validation(claim_type),
                    "review_required": row.get("confidence_tier") in {"exploratory", "low"} or bool(row.get("needs_validation")),
                    "source": "analysis_pack.pathway_rankings",
                }
            )
        return claims

    def interpretation_node_definition(self, node_type: str, name: str) -> str:
        lowered = name.lower()
        if node_type == "metabolite":
            return "来自输入表或 release 解析器的化学实体。"
        if node_type == "pathway":
            if "glycolysis" in lowered:
                return "葡萄糖分解相关的核心能量代谢通路，连接葡萄糖利用、丙酮酸和乳酸相关代谢。"
            if "tca" in lowered or "citric" in lowered:
                return "线粒体碳代谢通路，连接有机酸循环和能量产生。"
            if "glutamine" in lowered or "glutamate" in lowered:
                return "谷氨酰胺/谷氨酸相关的氮代谢和碳代谢主题，常与补充反应和氨基酸交换有关。"
            return "用于组织相关代谢物、基因和反应的 curated 通路或通路样条目。"
        if node_type in {"target", "gene"}:
            return "由代谢物/通路图谱或文献 overlay 连接到的基因或靶点节点。"
        if node_type == "disease":
            return "疾病或表型背景节点，应解读为背景关联，而不是诊断。"
        if node_type == "drug":
            return "来自 overlay 证据的研究用扰动或药物候选。"
        return "用于连接输入、证据和候选解释的图谱节点。"

    def interpretation_node_card(
        self,
        node_uid: str,
        display_name: str,
        node_type: str,
        why: str,
        relation: str,
        confidence_tier: str,
        confidence_reason: str,
        evidence_refs: list[dict[str, Any]] | None = None,
        boundary: str = "",
    ) -> dict[str, Any]:
        node = self.get_node(node_uid) if node_uid else None
        xrefs = parse_list((node or {}).get("external_xrefs"))[:10]
        primary_id = (node or {}).get("primary_external_id", "")
        evidence_count = len(evidence_refs or [])
        stable_ids = [value for value in [primary_id, *xrefs] if value]
        entity_label = display_name or (node or {}).get("display_name", "") or node_uid
        return {
            "node_uid": node_uid,
            "display_name": entity_label,
            "node_type": node_type,
            "entity_summary": {
                "label": entity_label,
                "type": node_type,
                "primary_id": primary_id,
                "node_uid": node_uid,
                "stable_ids": stable_ids[:6],
                "evidence_ref_count": evidence_count,
                "confidence_tier": confidence_tier,
                "is_precise_entity": bool(node_uid and node_type not in {"metabolic_theme", "disease", "drug"}),
                "display_hint": "精确图谱实体" if node_uid else "主题或候选实体",
            },
            "plain_language_definition": self.interpretation_node_definition(node_type, entity_label),
            "role_in_network": {
                "metabolite": "输入或解析种子",
                "pathway": "通路/主题组织节点",
                "target": "机制候选",
                "gene": "机制候选",
                "disease": "背景关联",
                "drug": "研究用 overlay 假设",
            }.get(node_type, "图谱连接节点"),
            "why_in_this_analysis": why,
            "relation_to_input": relation,
            "evidence_source_summary": "curated 图谱和关联文献引用" if evidence_refs else "curated 图谱或解析器元数据",
            "confidence_tier": confidence_tier,
            "confidence_reason": confidence_reason,
            "boundary": boundary or self.interpretation_boundary({}, node_type),
            "open_details": {
                "primary_external_id": primary_id,
                "xrefs": xrefs,
                "evidence_refs": evidence_refs or [],
            },
        }

    def interpretation_node_cards(
        self,
        analysis_pack: dict[str, Any],
        prediction_pack: dict[str, Any],
        claims: list[dict[str, Any]],
    ) -> list[dict[str, Any]]:
        cards: list[dict[str, Any]] = []
        seen: set[tuple[str, str]] = set()

        def add(card: dict[str, Any]) -> None:
            key = (str(card.get("node_type") or ""), str(card.get("node_uid") or card.get("display_name") or ""))
            if key not in seen:
                seen.add(key)
                cards.append(card)

        for row in (analysis_pack.get("matched") or [])[:6]:
            uid = str(row.get("metabolite_uid") or "")
            label = self.interpretation_candidate_label(row) or self.interpretation_input_label(row)
            add(
                self.interpretation_node_card(
                    uid,
                    label,
                    "metabolite",
                    "由用户输入解析得到，并作为分析种子使用。",
                    f"输入行：{self.interpretation_input_label(row)}",
                    "high" if row.get("match_status") == "matched" else "exploratory",
                    f"match score {row.get('match_score', row.get('top_score', 0.0))}",
                    [],
                    "只有 strict matched 的代谢物才支持精确化学身份层面的解释。",
                )
            )
        for claim in claims[:5]:
            for node in claim.get("supporting_nodes", []) or []:
                node_type = str(node.get("node_type") or "")
                if node_type == "metabolic_theme":
                    continue
                add(
                    self.interpretation_node_card(
                        str(node.get("node_uid") or ""),
                        str(node.get("display_name") or ""),
                        node_type,
                        f"属于核心结论 {claim.get('claim_id')} 的关键节点。",
                        "通过输入支持和证据引用与本次分析相连。",
                        str(claim.get("confidence_tier") or "exploratory"),
                        "; ".join((claim.get("confidence_reasons") or {}).get("positive_factors", [])[:2]) or "ranked in analysis pack",
                        claim.get("evidence_refs", []),
                        str(claim.get("boundary") or ""),
                    )
                )
        for row in (analysis_pack.get("target_rankings") or [])[:3]:
            add(
                self.interpretation_node_card(
                    str(row.get("target_uid") or ""),
                    str(row.get("display_name") or row.get("target_uid") or ""),
                    "target",
                    "由图传播或文献支持排序得到。",
                    "位于代谢物/通路信号下游的候选机制连接。",
                    str(row.get("confidence_tier") or "exploratory"),
                    str(row.get("downgrade_reason") or "target ranking is hypothesis-generating"),
                    row.get("evidence_refs", []),
                    self.interpretation_boundary(row, "target"),
                )
            )
        for row in (analysis_pack.get("disease_rankings") or [])[:2]:
            add(
                self.interpretation_node_card(
                    str(row.get("disease_uid") or ""),
                    str(row.get("display_name") or row.get("disease_uid") or ""),
                    "disease",
                    "作为疾病或表型关联被连接到。",
                    "仅表示背景关联，不是样本诊断。",
                    str(row.get("confidence_tier") or "exploratory"),
                    str(row.get("downgrade_reason") or "disease rows are shown as context hypotheses"),
                    row.get("evidence_refs", []),
                    self.interpretation_boundary(row, "disease"),
                )
            )
        for row in (prediction_pack.get("drug_rankings") or [])[:2]:
            add(
                self.interpretation_node_card(
                    str(row.get("drug_id") or ""),
                    str(row.get("drug_name") or row.get("drug_id") or ""),
                    "drug",
                    "通过可选 prediction overlay 匹配得到。",
                    "与排序靶点相关的研究优先级候选。",
                    str(row.get("display_confidence") or row.get("confidence_tier") or "low"),
                    str(row.get("downgrade_reason") or "drug hypotheses inherit upstream target confidence"),
                    row.get("evidence_refs", []),
                    self.interpretation_boundary(row, "drug"),
                )
            )
        return cards[:20]

    def build_interpretation_report(
        self,
        analysis_pack: dict[str, Any],
        prediction_pack: dict[str, Any] | None = None,
        question: str = "",
        context: Any = None,
    ) -> dict[str, Any]:
        prediction_pack = prediction_pack or {}
        input_summary = analysis_pack.get("input_summary") or {}
        support_index = self.interpretation_input_support_index(analysis_pack)
        theme_claims = self.interpretation_theme_claims(analysis_pack, support_index, input_summary)
        pathway_claims = self.interpretation_pathway_claims(analysis_pack, support_index, input_summary, len(theme_claims))
        v2_claims = []
        for index, fact in enumerate(((analysis_pack.get("database_accuracy") or {}).get("mechanism_ready_facts") or [])[:3], start=1):
            subject = str(fact.get("subject_name") or fact.get("subject_uid") or "input chemical")
            obj = str(fact.get("object_name") or fact.get("object_uid") or "mechanism fact")
            predicate = str(fact.get("predicate") or fact.get("fact_type") or "supports")
            v2_claims.append(
                {
                    "claim_id": f"v2_fact_{index:03d}",
                    "headline": f"{subject} has traceable mechanism fact: {obj}",
                    "confidence_tier": "high" if fact.get("source_record_uid") else "medium",
                    "calibrated_confidence": 0.85 if fact.get("source_record_uid") else 0.7,
                    "supporting_input_count": 1,
                    "supporting_inputs": [],
                    "basis": f"database_accuracy_store.v2 fact {fact.get('fact_uid', '')} links {subject} via {predicate} to {obj}.",
                    "relation_chain": [subject, predicate, obj],
                    "possible_mechanism": "This is a mechanism-ready database fact anchored to identity decisions and reaction/module/evidence source records.",
                    "boundary": "A mechanism-ready fact is a curated factual support unit; sample activation or causality still requires experimental effect direction and context validation.",
                    "confidence_reasons": {
                        "positive_factors": ["accepted exact identity", "v2 mechanism-ready source trace"],
                        "downgrade_factors": [] if fact.get("source_record_uid") else ["source record trace is incomplete"],
                    },
                    "evidence_refs": [
                        {
                            "ref_type": "database_accuracy_fact",
                            "source_record_id": fact.get("source_record_uid", ""),
                            "edge_uid": fact.get("evidence_assertion_uid", ""),
                            "source_release": self.release_id,
                        }
                    ],
                    "next_validation": [
                        "Check whether the reaction/module fact has consistent direction in the original sample statistics.",
                        "Validate representative metabolites and enzymes with orthogonal assays.",
                    ],
                }
            )
        claims = [*v2_claims, *theme_claims, *pathway_claims]

        def claim_sort_key(claim: dict[str, Any]) -> tuple[int, float, int, str]:
            tier_rank = {"high": 0, "medium": 1, "exploratory": 2, "low": 3}.get(str(claim.get("confidence_tier")), 4)
            confidence = float(claim.get("calibrated_confidence") or 0.0)
            support_count = int(claim.get("supporting_input_count") or len(claim.get("supporting_inputs") or []) or 0)
            return (tier_rank, -confidence, -support_count, str(claim.get("claim_id") or ""))

        claims = sorted(claims, key=claim_sort_key)[:5]
        for index, claim in enumerate(claims, start=1):
            claim["rank"] = index
            claim["claim_id"] = f"claim_{index:03d}"

        review_queue = []
        for bucket in ("ambiguous", "unmatched"):
            for row in analysis_pack.get(bucket, []) or []:
                review_queue.append(
                    {
                        "queue_type": bucket,
                        "input_id": row.get("input_id", ""),
                        "input_name": self.interpretation_input_label(row),
                        "reason": ";".join(row.get("identity_review_reasons") or []) or row.get("resolution_status", bucket),
                        "suggested_action": "补充稳定的 HMDB/ChEBI/PubChem/InChIKey 标识符，或提供 RT/MS2 等正交 LC-MS 证据。",
                    }
                )

        evidence_refs = self.analysis_pack_dedupe_refs(
            [ref for claim in claims for ref in (claim.get("evidence_refs") or [])],
            limit=50,
        )
        report = {
            "report_version": INTERPRETATION_REPORT_CONTRACT_VERSION,
            "question": question,
            "mode": "context_aware" if context else "generalized",
            "database_accuracy": {
                "contract_version": (analysis_pack.get("database_accuracy") or {}).get("contract_version", DATABASE_ACCURACY_CONTRACT_VERSION),
                "store_available": bool((analysis_pack.get("database_accuracy") or {}).get("store_available")),
                "identity_decision_summary": (analysis_pack.get("database_accuracy") or {}).get("identity_decision_summary", {}),
                "mechanism_ready_facts": (analysis_pack.get("database_accuracy") or {}).get("mechanism_ready_facts", [])[:10],
                "blocked_input_count": len((analysis_pack.get("database_accuracy") or {}).get("blocked_inputs") or []),
            },
            "conclusion_evaluation": analysis_pack.get("conclusion_evaluation") or {
                "contract_version": CONCLUSION_EVALUATION_CONTRACT_VERSION,
                "gold_standard": {"entry_count": 0, "files": []},
                "metrics": {},
                "bug_queue": [],
            },
            "executive_summary": [
                {
                    "claim_id": claim.get("claim_id"),
                    "headline": claim.get("headline"),
                    "confidence_tier": claim.get("confidence_tier"),
                    "supporting_input_count": int(claim.get("supporting_input_count") or len(claim.get("supporting_inputs") or []) or 0),
                    "evidence_ref_count": len(claim.get("evidence_refs") or []),
                    "review_required": bool(claim.get("review_required")),
                    "boundary": claim.get("boundary"),
                }
                for claim in claims
            ],
            "conclusion_chains": claims,
            "node_cards": self.interpretation_node_cards(analysis_pack, prediction_pack, claims),
            "confidence_explanations": [
                {
                    "claim_id": claim.get("claim_id"),
                    "confidence_tier": claim.get("confidence_tier"),
                    "calibrated_confidence": claim.get("calibrated_confidence"),
                    "positive_factors": (claim.get("confidence_reasons") or {}).get("positive_factors", []),
                    "downgrade_factors": (claim.get("confidence_reasons") or {}).get("downgrade_factors", []),
                }
                for claim in claims
            ],
            "review_queue": review_queue[:50],
            "evidence_sections": [
                {
                    "section": "core_conclusion_evidence",
                    "description": "Evidence refs attached to the core conclusion chains.",
                    "evidence_refs": evidence_refs,
                }
            ],
            "next_validation": [
                {
                    "claim_id": claim.get("claim_id"),
                    "headline": claim.get("headline"),
                    "suggestions": claim.get("next_validation", []),
                }
                for claim in claims
            ],
            "appendix": {
                "legacy_rankings_are_supporting_layer": True,
                "raw_ranking_counts": {
                    "pathway_rankings": len(analysis_pack.get("pathway_rankings") or []),
                    "target_rankings": len(analysis_pack.get("target_rankings") or []),
                    "disease_rankings": len(analysis_pack.get("disease_rankings") or []),
                },
                "interpretation_boundary": "This report organizes frozen analysis outputs for research interpretation. It does not create graph facts or clinical recommendations.",
            },
        }
        report["determinism"] = {
            "release_id": self.release_id,
            "analysis_pack_hash": (analysis_pack.get("determinism") or {}).get("analysis_pack_hash", ""),
            "report_hash": content_hash({key: value for key, value in report.items() if key != "determinism"}),
        }
        return report

    def build_analysis_pack(
        self,
        request: dict[str, Any],
        precheck: dict[str, Any],
        analysis_features: list[dict[str, Any]],
        features_by_uid: dict[str, list[dict[str, Any]]],
        seed_weights: dict[str, float],
        pathways: list[dict[str, Any]],
        targets: list[dict[str, Any]],
        diseases: list[dict[str, Any]],
        paths: list[dict[str, Any]],
        propagation: dict[str, Any],
        context: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        paths_by_id = {path.get("path_id", ""): path for path in paths if path.get("path_id")}
        expanded_summary = precheck.get("expanded_summary", {}) or {}
        input_summary = {
            "analysis_mode": precheck.get("analysis_mode", "metabolite_table"),
            "input_count": precheck.get("input_count", len(request.get("records", []))),
            "matched_count": precheck.get("summary", {}).get("matched", 0),
            "ambiguous_count": precheck.get("summary", {}).get("ambiguous", 0),
            "unmatched_count": precheck.get("summary", {}).get("unmatched", 0),
            "invalid_count": precheck.get("summary", {}).get("invalid", 0),
            "strict_matched_count": expanded_summary.get("strict_matched_count", precheck.get("summary", {}).get("matched", 0)),
            "expanded_candidate_count": expanded_summary.get("expanded_candidate_count", len(precheck.get("expanded_candidates", []) or [])),
            "expanded_candidate_input_count": expanded_summary.get("expanded_candidate_input_count", expanded_summary.get("expanded_input_count", 0)),
            "expanded_input_count": expanded_summary.get("expanded_input_count", 0),
            "expanded_by_class": expanded_summary.get("expanded_by_class", {}),
            "unresolved_input_count": expanded_summary.get("unresolved_input_count", 0),
            "analysis_seed_input_count": expanded_summary.get("analysis_seed_input_count", precheck.get("summary", {}).get("matched", 0)),
            "analysis_seed_feature_count": expanded_summary.get(
                "analysis_seed_feature_count",
                precheck.get("summary", {}).get("matched", 0) + len(precheck.get("expanded_candidates", []) or []),
            ),
            "weighted_seed_mass": expanded_summary.get("weighted_seed_mass", 0.0),
            "seed_weight_policy": expanded_summary.get("weight_policy", EXPANDED_SEED_WEIGHTS),
        }
        genetic_exposures = self.analysis_pack_genetic_exposures(precheck)
        input_summary.update(
            {
                "genetic_exposure_count": len(genetic_exposures),
                "genetic_exposure_available_count": sum(
                    1 for row in genetic_exposures if row.get("genetic_exposure_status") == "available"
                ),
            }
        )
        four_track_summary = self.analysis_pack_four_track_summary(input_summary, genetic_exposures, expanded_summary)
        pathway_rankings = self.analysis_pack_ranking_rows(
            pathways, "pathway_uid", paths_by_id, features_by_uid, input_summary, context
        )
        target_rankings = self.analysis_pack_ranking_rows(
            targets, "target_uid", paths_by_id, features_by_uid, input_summary, context
        )
        disease_rankings = self.analysis_pack_ranking_rows(
            diseases, "disease_uid", paths_by_id, features_by_uid, input_summary, context
        )
        selected_paths = self.select_analysis_pack_top_paths(paths, [pathway_rankings, target_rankings, disease_rankings])
        top_paths = [self.analysis_pack_path_row(path) for path in selected_paths]
        fallback_chains = []
        if not top_paths:
            fallback_chains = self.analysis_pack_fallback_explanation_chains(
                pathway_rankings,
                target_rankings,
                disease_rankings,
                input_summary,
            )
        evidence_refs = self.analysis_pack_dedupe_refs(
            [
                ref
                for section in (pathway_rankings, target_rankings, disease_rankings, top_paths, fallback_chains)
                for row in section
                for ref in row.get("evidence_refs", [])
            ]
        )
        release = self.release_meta()
        feature_summary = self.feature_summary(analysis_features)
        duplicate_uids = sorted(uid for uid, features in features_by_uid.items() if len(features) > 1)
        quality_warnings = self.analysis_pack_quality_warnings(precheck, analysis_features, features_by_uid, propagation, paths)
        database_accuracy = self.database_accuracy_pack(precheck, features_by_uid, context=context)
        pack = {
            "contract_version": ANALYSIS_PACK_CONTRACT_VERSION,
            "input_summary": {
                "analysis_mode": precheck.get("analysis_mode", "metabolite_table"),
                **input_summary,
                "unique_matched_metabolite_count": len(features_by_uid),
                "duplicate_matched_metabolite_uids": duplicate_uids,
                "feature_summary": feature_summary,
                "seed_weights": seed_weights,
            },
            "matched": [self.analysis_pack_resolution_row(row, "matched") for row in precheck.get("matched", [])],
            "expanded_candidates": [
                self.analysis_pack_resolution_row(row, "expanded_candidate") for row in precheck.get("expanded_candidates", [])
            ],
            "expanded_summary": expanded_summary,
            "genetic_exposures": genetic_exposures,
            "four_track_summary": four_track_summary,
            "ambiguous": [self.analysis_pack_resolution_row(row, "ambiguous") for row in precheck.get("ambiguous", [])],
            "unmatched": [self.analysis_pack_resolution_row(row, "unmatched") for row in precheck.get("unmatched", [])]
            + [self.analysis_pack_resolution_row(row, "invalid") for row in precheck.get("invalid", [])],
            "biological_entity_pools": self.build_biological_entity_pool_summary(precheck),
            "database_accuracy": database_accuracy,
            "quality_warnings": quality_warnings,
            "pathway_rankings": pathway_rankings,
            "target_rankings": target_rankings,
            "disease_rankings": disease_rankings,
            "top_explanation_paths": top_paths,
            "fallback_explanation_chains": fallback_chains,
            "literature_evidence_pack": self.analysis_pack_literature_summary(evidence_refs),
            "evidence_refs": evidence_refs,
            "release": release,
            "blocked_reasons": self.analysis_pack_blocked_reasons(precheck, pathways, targets, diseases),
        }
        pack["conclusion_evaluation"] = self.build_conclusion_evaluation(pack, context)
        pack["prediction_model"] = self.build_prediction_model_pack(
            precheck,
            features_by_uid,
            request.get("records", []),
            input_summary,
            pathway_rankings,
            target_rankings,
            disease_rankings,
            quality_warnings,
            context,
        )
        pack["determinism"] = {
            "input_hash": content_hash(request),
            "release_id": self.release_id,
            "config_hash": release.get("config_hash", ""),
            "contract_hash": content_hash({"contract_version": ANALYSIS_PACK_CONTRACT_VERSION, "fields": sorted(pack.keys())})[:16],
        }
        pack["determinism"]["analysis_pack_hash"] = content_hash({key: value for key, value in pack.items() if key != "determinism"})
        return pack

    def analyze_metabolites(
        self,
        records: list[Any],
        max_paths: int | None = None,
        max_hops: int | None = None,
        context: Any = None,
        context_mode: str = "soft",
    ) -> dict[str, Any]:
        records, input_normalization = prepare_analysis_records(records)
        request = {
            "records": records,
            "max_paths": max_paths,
            "max_hops": max_hops,
            "context": context,
            "context_mode": context_mode,
            "input_normalization": input_normalization,
        }
        precheck = self.precheck_metabolites(records, input_normalization=input_normalization)
        seed_rows = [*precheck.get("matched", []), *precheck.get("expanded_candidates", [])]
        analysis_features, features_by_uid = self.matched_analysis_features(seed_rows)
        seed_weights = self.seed_weights_from_features(features_by_uid)
        seed_presence_weights = self.seed_presence_weights_from_features(features_by_uid)
        matched_uids = sorted(features_by_uid)
        enrichment_rows = self.pathway_enrichment(matched_uids, seed_presence_weights)
        paths = self.enumerate_paths(matched_uids, max_hops=max_hops, max_paths=max_paths)
        propagation = self.propagate_scores(seed_weights, max_hops=max_hops)
        best_paths = self.best_paths_by_terminal(paths)
        pathways = self.merge_pathway_rankings(enrichment_rows, propagation, best_paths, features_by_uid, input_records=records)
        targets = self.propagation_ranking(propagation, "target", "target_uid", best_paths)
        diseases = self.propagation_ranking(propagation, "disease", "disease_uid", best_paths)
        directional_enrichment = self.directional_pathway_enrichment(features_by_uid)
        normalized_context = self.normalize_prediction_context(context, records)
        analysis_pack = self.build_analysis_pack(
            request,
            precheck,
            analysis_features,
            features_by_uid,
            seed_weights,
            pathways[:50],
            targets[:50],
            diseases[:50],
            paths,
            propagation,
            context=normalized_context,
        )
        prediction_pack = self.build_prediction_pack(
            records,
            analysis_features,
            analysis_pack.get("pathway_rankings", pathways[:50]),
            analysis_pack.get("target_rankings", targets[:50]),
            analysis_pack.get("disease_rankings", diseases[:50]),
            context=normalized_context,
            context_mode=context_mode,
        )
        structured_prediction = self.structured_prediction_json(analysis_pack, prediction_pack, normalized_context)
        interpretation_report = self.build_interpretation_report(
            analysis_pack,
            prediction_pack,
            question="",
            context=normalized_context,
        )
        payload = {
            "analysis_pack": analysis_pack,
            "structured_prediction": structured_prediction,
            "interpretation_report": interpretation_report,
            "precheck": {
                "analysis_mode": precheck.get("analysis_mode", "metabolite_table"),
                "input_normalization": precheck.get("input_normalization", input_normalization),
                "summary": precheck["summary"],
                "matched": precheck["matched"],
                "expanded_candidates": precheck.get("expanded_candidates", []),
                "expanded_summary": precheck.get("expanded_summary", {}),
                "identity_resolution_v2_summary": precheck.get("identity_resolution_v2_summary", {}),
                "ambiguous": precheck["ambiguous"],
                "unmatched": precheck["unmatched"],
                "invalid": precheck["invalid"],
            },
            "input_normalization": input_normalization,
            "analysis_features": {
                "summary": self.feature_summary(analysis_features),
                "matched": analysis_features,
                "seed_weights": seed_weights,
            },
            "rankings": {
                "pathways": pathways[:50],
                "targets": targets[:50],
                "diseases": diseases[:50],
            },
            "directional_enrichment": directional_enrichment,
            "predictions": prediction_pack,
            "propagation": {
                "formula": propagation["formula"],
                "alpha": propagation["alpha"],
                "iterations": propagation["iterations"],
                "converged": propagation["converged"],
                "graph": propagation["graph"],
                "top_nodes": propagation["nodes"][:50],
            },
            "explanation_paths": paths,
            "compressed_explanation_paths": self.compress_explanation_paths(paths),
            "score_notes": {
                "entity_resolution": "matched only when top1 >= 60 and top1 - top2 >= 15",
                "expanded_resolution": "ambiguous/unmatched rows may contribute low-weight expanded seeds as soft_identity, analog_candidate, class_or_pool, or ratio_component; ratio components are exploratory and require validation before abundance interpretation",
                "relation_probability": "p_final(e)=1-(1-p_curated)(1-p_literature)(1-p_topology)(1-p_user)",
                "path_cost": "C(path)=sum(-log(p_final(edge)))",
                "propagation": "pi=(1-alpha)y+alpha*W*pi over a compressed typed beam graph induced by strict and low-weight expanded metabolite seeds",
                "ranking_score": "pathway score=enrichment+propagation+literature_boost; target/disease score=propagation+literature_boost",
                "compression": "typed beam keeps high-mass transitions; path signatures group explanation paths by node and edge type",
                "user_probability": "p_user for uploaded differential features is reported on seeds and used to weight y",
                "literature_probability": "confirm/support_direction support can boost existing graph-edge p_final; novel_candidate rows stay evidence-labeled and can support explanations, but do not create graph-scoring edges",
            },
        }
        return self.response_envelope("/analyze/metabolites", request, payload)

    def analyze_trait_score(
        self,
        records: list[Any],
        max_paths: int | None = None,
        max_hops: int | None = None,
        context: Any = None,
        context_mode: str = "soft",
        q_threshold: float | None = 0.05,
        min_abs_effect: float = 0.0,
        top_per_celltype: int = 80,
        max_records: int = 120,
    ) -> dict[str, Any]:
        selected_records, selection = select_trait_score_records(
            records,
            q_threshold=q_threshold,
            min_abs_effect=min_abs_effect,
            top_per_celltype=top_per_celltype,
            max_records=max_records,
        )
        effective_context = context or infer_trait_score_context(selected_records)
        analyzed = self.analyze_metabolites(
            selected_records,
            max_paths=max_paths,
            max_hops=max_hops,
            context=effective_context,
            context_mode=context_mode,
        )
        request = {
            "input_count": len(records),
            "selected_count": len(selected_records),
            "max_paths": max_paths,
            "max_hops": max_hops,
            "context": effective_context,
            "context_mode": context_mode,
            "q_threshold": q_threshold,
            "min_abs_effect": min_abs_effect,
            "top_per_celltype": top_per_celltype,
            "max_records": max_records,
        }
        payload = {
            "trait_score_selection": selection,
            "selected_records": selected_records,
            "analysis_pack": analyzed.get("analysis_pack", {}),
            "structured_prediction": analyzed.get("structured_prediction", {}),
            "interpretation_report": analyzed.get("interpretation_report", {}),
            "precheck": analyzed.get("precheck", {}),
            "analysis_features": analyzed.get("analysis_features", {}),
            "rankings": analyzed.get("rankings", {}),
            "directional_enrichment": analyzed.get("directional_enrichment", {}),
            "predictions": analyzed.get("predictions", {}),
            "propagation": analyzed.get("propagation", {}),
            "explanation_paths": analyzed.get("explanation_paths", []),
            "compressed_explanation_paths": analyzed.get("compressed_explanation_paths", {}),
            "score_notes": {
                **(analyzed.get("score_notes", {}) or {}),
                "trait_score_selection": "TraitScore rows are selected by q value and effect size before GCST-to-chemical mapping; outputs are research-prioritization signals, not direct abundance claims.",
            },
            "source_analysis_response_hash": (analyzed.get("determinism") or {}).get("response_hash", ""),
        }
        return self.response_envelope("/analyze/trait-score", request, payload)

    def analyze_differential_table(
        self,
        records: list[Any],
        max_paths: int | None = None,
        max_hops: int | None = None,
        context: Any = None,
        context_mode: str = "soft",
        q_threshold: float | None = 0.05,
        p_threshold: float | None = None,
        min_abs_effect: float = 0.0,
        top_per_group: int = 120,
        max_records: int = 300,
    ) -> dict[str, Any]:
        selected_records, selection = select_differential_table_records(
            records,
            q_threshold=q_threshold,
            p_threshold=p_threshold,
            min_abs_effect=min_abs_effect,
            top_per_group=top_per_group,
            max_records=max_records,
        )
        effective_context = context or infer_differential_table_context(selected_records)
        analyzed = self.analyze_metabolites(
            selected_records,
            max_paths=max_paths,
            max_hops=max_hops,
            context=effective_context,
            context_mode=context_mode,
        )
        request = {
            "input_count": len(records),
            "selected_count": len(selected_records),
            "max_paths": max_paths,
            "max_hops": max_hops,
            "context": effective_context,
            "context_mode": context_mode,
            "q_threshold": q_threshold,
            "p_threshold": p_threshold,
            "min_abs_effect": min_abs_effect,
            "top_per_group": top_per_group,
            "max_records": max_records,
        }
        payload = {
            "differential_table_selection": selection,
            "selected_records": selected_records,
            "analysis_pack": analyzed.get("analysis_pack", {}),
            "structured_prediction": analyzed.get("structured_prediction", {}),
            "interpretation_report": analyzed.get("interpretation_report", {}),
            "precheck": analyzed.get("precheck", {}),
            "analysis_features": analyzed.get("analysis_features", {}),
            "rankings": analyzed.get("rankings", {}),
            "directional_enrichment": analyzed.get("directional_enrichment", {}),
            "predictions": analyzed.get("predictions", {}),
            "propagation": analyzed.get("propagation", {}),
            "explanation_paths": analyzed.get("explanation_paths", []),
            "compressed_explanation_paths": analyzed.get("compressed_explanation_paths", {}),
            "score_notes": {
                **(analyzed.get("score_notes", {}) or {}),
                "differential_table_selection": (
                    "Input rows are treated as precomputed differential results. The adapter selects rows by "
                    "adjusted p/q or p value and effect size; it does not require or read raw abundance matrices."
                ),
                "effect_semantics": (
                    "When log2FC is absent, mean_diff, cohen_d, z statistics, or other directional effects are "
                    "used only as seed-weight surrogates and remain labeled as non-abundance effect fields."
                ),
            },
            "source_analysis_response_hash": (analyzed.get("determinism") or {}).get("response_hash", ""),
        }
        return self.response_envelope("/analyze/differential-table", request, payload)

    def llm_adapter_evidence_request(self, analyzed: dict[str, Any], evidence_limit: int | None = None) -> dict[str, Any]:
        limit = min(100, max(1, int(evidence_limit or self.config.max_evidence_items)))
        pack = analyzed.get("analysis_pack", {})
        for ref in pack.get("evidence_refs", []) or []:
            edge_uid = str(ref.get("edge_uid") or "")
            if edge_uid:
                return {"edge_uid": edge_uid, "limit": limit}
        for row in pack.get("top_explanation_paths", []) or []:
            claim_refs = row.get("claim_refs", {})
            edge_uids = claim_refs.get("edge_uids") or []
            if edge_uids:
                return {"edge_uid": str(edge_uids[0]), "limit": limit}
        for ranking_key, uid_key in (
            ("target_rankings", "target_uid"),
            ("disease_rankings", "disease_uid"),
            ("pathway_rankings", "pathway_uid"),
        ):
            rows = pack.get(ranking_key) or []
            if rows and rows[0].get(uid_key):
                return {"entity_uid": str(rows[0][uid_key]), "limit": limit}
        matched = pack.get("matched") or []
        if matched and matched[0].get("metabolite_uid"):
            return {"entity_uid": str(matched[0]["metabolite_uid"]), "limit": limit}
        return {"limit": limit}

    def llm_adapter_subgraph_seed_ids(self, analyzed: dict[str, Any], seed_limit: int = 6) -> list[str]:
        pack = analyzed.get("analysis_pack", {})
        seeds: list[str] = []

        def add_seed(value: Any) -> None:
            text = str(value or "")
            if text and text not in seeds and len(seeds) < seed_limit:
                seeds.append(text)

        for row in pack.get("matched", []) or []:
            add_seed(row.get("metabolite_uid"))
        for row in pack.get("top_explanation_paths", []) or []:
            add_seed(row.get("terminal_node_uid"))
        for ranking_key, uid_key in (
            ("pathway_rankings", "pathway_uid"),
            ("target_rankings", "target_uid"),
            ("disease_rankings", "disease_uid"),
        ):
            for row in pack.get(ranking_key, []) or []:
                add_seed(row.get(uid_key))
                if len(seeds) >= seed_limit:
                    break
        return seeds

    def active_release_payload(self, releases: dict[str, Any]) -> dict[str, Any]:
        active = None
        for row in releases.get("releases", []) or []:
            if row.get("release_id") == self.release_id:
                active = row
                break
        return {
            "release_id": self.release_id,
            "service_release": self.release_meta(),
            "active_release": active or {},
            "releases": releases.get("releases", []),
            "source_endpoint": "/releases",
            "source_response_hash": releases.get("determinism", {}).get("response_hash", ""),
        }

    def build_llm_adapter_input_pack(
        self,
        records: list[Any],
        question: str = "",
        max_paths: int | None = None,
        max_hops: int | None = None,
        evidence_limit: int | None = None,
        subgraph_max_hops: int = 1,
        context: Any = None,
        context_mode: str = "soft",
        analyzed_response: dict[str, Any] | None = None,
        analysis_endpoint: str = "/analyze/metabolites",
    ) -> dict[str, Any]:
        analyzed = analyzed_response or self.analyze_metabolites(
            records, max_paths=max_paths, max_hops=max_hops, context=context, context_mode=context_mode
        )
        evidence_request = self.llm_adapter_evidence_request(analyzed, evidence_limit=evidence_limit)
        evidence = self.evidence(**evidence_request)
        seed_ids = self.llm_adapter_subgraph_seed_ids(analyzed)
        subgraph = self.subgraph(seed_ids, max_hops=subgraph_max_hops) if seed_ids else self.subgraph([], max_hops=subgraph_max_hops)
        releases = self.releases()
        analysis_payload = dict(analyzed.get("analysis_pack", {}) or {})
        structured_prediction = analyzed.get("structured_prediction", {})
        if structured_prediction:
            analysis_payload["structured_prediction"] = structured_prediction
        interpretation_report = analyzed.get("interpretation_report", {})
        if interpretation_report:
            analysis_payload["interpretation_report"] = interpretation_report
        for key in (
            "trait_score_selection",
            "differential_table_selection",
            "analysis_features",
            "directional_enrichment",
            "score_notes",
        ):
            value = analyzed.get(key)
            if value:
                analysis_payload[key] = value
        source_endpoints = [
            {
                "endpoint": analysis_endpoint,
                "payload_key": "analysis_pack",
                "request_hash": analyzed.get("determinism", {}).get("request_hash", ""),
                "response_hash": analyzed.get("determinism", {}).get("response_hash", ""),
            },
            {
                "endpoint": "/evidence",
                "payload_key": "evidence",
                "request_hash": evidence.get("determinism", {}).get("request_hash", ""),
                "response_hash": evidence.get("determinism", {}).get("response_hash", ""),
            },
            {
                "endpoint": "/subgraph",
                "payload_key": "subgraph",
                "request_hash": subgraph.get("determinism", {}).get("request_hash", ""),
                "response_hash": subgraph.get("determinism", {}).get("response_hash", ""),
            },
            {
                "endpoint": "/releases",
                "payload_key": "release",
                "request_hash": releases.get("determinism", {}).get("request_hash", ""),
                "response_hash": releases.get("determinism", {}).get("response_hash", ""),
            },
        ]
        return {
            "contract_version": LLM_ADAPTER_INPUT_CONTRACT_VERSION,
            "adapter_intent": "analysis_explanation",
            "source_endpoints": source_endpoints,
            "request": {
                "question": str(question or ""),
                "max_paths": max_paths,
                "max_hops": max_hops,
                "context": context,
                "context_mode": context_mode,
                "evidence_request": evidence_request,
                "subgraph_seed_ids": seed_ids,
                "subgraph_max_hops": subgraph_max_hops,
                "analysis_endpoint": analysis_endpoint,
            },
            "analysis_pack": analysis_payload,
            "prediction_pack": analyzed.get("predictions", {}),
            "evidence": {
                key: value
                for key, value in evidence.items()
                if key not in {"api_version", "endpoint", "release", "determinism"}
            },
            "subgraph": {
                key: value
                for key, value in subgraph.items()
                if key not in {"api_version", "endpoint", "release", "determinism"}
            },
            "release": self.active_release_payload(releases),
        }

    def blocked_llm_adapter_output(
        self,
        adapter_output: dict[str, Any],
        guard: dict[str, Any],
        input_pack: dict[str, Any],
        backend_audit: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        output = {
            "contract_version": LLM_ADAPTER_OUTPUT_CONTRACT_VERSION,
            "adapter_version": str(adapter_output.get("adapter_version") or "external"),
            "status": "blocked_by_guard",
            "guard": guard,
            "insufficient_evidence": {
                "reason": "Adapter output failed the deterministic guard and was withheld.",
                "source_refs": [
                    {
                        "source_type": "release",
                        "ref_id": f"release:{self.release_id}",
                        "path": "$.release.release_id",
                    }
                ],
            },
        }
        output["determinism"] = {
            "adapter_input_hash": content_hash(input_pack),
            "rejected_adapter_output_hash": content_hash(adapter_output),
            "guard_policy_hash": guard_policy_hash(),
        }
        if backend_audit:
            output["determinism"]["backend_audit"] = backend_audit
        output["determinism"]["response_hash"] = content_hash({key: value for key, value in output.items() if key != "determinism"})
        return output

    def force_guard_issue(self, guard: dict[str, Any], issue: dict[str, Any]) -> dict[str, Any]:
        return {
            **guard,
            "passed": False,
            "issue_count": int(guard.get("issue_count", 0) or 0) + 1,
            "issues": [*guard.get("issues", []), issue],
        }

    def external_llm_fallback_output(
        self,
        input_pack: dict[str, Any],
        error: LLMAdapterError | Exception,
        backend_name: str,
        config: ExternalLLMConfig,
    ) -> tuple[dict[str, Any], dict[str, Any]]:
        issue = (
            error.as_issue()
            if isinstance(error, LLMAdapterError)
            else {
                "code": "external_llm_unhandled_error",
                "path": "$.adapter_backend",
                "detail": {"message": str(error), "type": type(error).__name__},
            }
        )
        explanation = build_local_adapter_output(input_pack)
        audit = {
            "backend": backend_name,
            "config": config.sanitized(),
            "fallback_backend": "local",
            "external_issue": issue,
        }
        explanation.setdefault("determinism", {})["backend_audit"] = audit
        explanation["determinism"]["response_hash"] = content_hash({key: value for key, value in explanation.items() if key != "determinism"})
        return explanation, issue

    def finalize_guarded_adapter_output(
        self,
        input_pack: dict[str, Any],
        raw_output: dict[str, Any],
        guard: dict[str, Any],
        backend_audit: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        if not guard.get("passed"):
            return self.blocked_llm_adapter_output(raw_output, guard, input_pack, backend_audit=backend_audit)
        final = copy.deepcopy(raw_output)
        final["guard"] = guard
        final["determinism"] = {
            "adapter_input_hash": content_hash(input_pack),
            "adapter_output_hash": content_hash(raw_output),
            "guard_policy_hash": guard_policy_hash(),
        }
        if backend_audit:
            final["determinism"]["backend_audit"] = backend_audit
        final["determinism"]["response_hash"] = content_hash({key: value for key, value in final.items() if key != "determinism"})
        return final

    def external_llm_error_output(self, error: LLMAdapterError | Exception) -> dict[str, Any]:
        return {
            "contract_version": LLM_ADAPTER_OUTPUT_CONTRACT_VERSION,
            "adapter_version": EXTERNAL_ADAPTER_VERSION,
            "status": "blocked_by_guard",
            "insufficient_evidence": {
                "reason": "External LLM backend did not produce a usable adapter output. Check the guard issue code and service configuration.",
                "source_refs": [
                    {
                        "source_type": "release",
                        "ref_id": f"release:{self.release_id}",
                        "path": "$.release.release_id",
                    }
                ],
            },
        }

    def llm_status(self) -> dict[str, Any]:
        config = self.llm_config
        issues: list[dict[str, Any]] = []
        ready = True
        try:
            config.validate_ready()
        except LLMAdapterError as exc:
            ready = False
            issues.append(exc.as_issue())
        return self.response_envelope(
            "/llm/status",
            {},
            {
                "llm": {
                    "ready": ready,
                    "config": config.sanitized(),
                    "has_api_key": bool(config.api_key),
                    "issues": issues,
                    "request_scoped_config_supported": True,
                    "request_scoped_config_note": "POST /chat, /explain, and /llm/test may include llm_config with endpoint, model, api_key, timeout_seconds, max_output_tokens, and proxy_url.",
                }
            },
        )

    def llm_test_diagnostic(self, code: str, detail: dict[str, Any]) -> str:
        status = detail.get("status")
        try:
            status_int = int(status) if status is not None else None
        except (TypeError, ValueError):
            status_int = None
        reason = str(detail.get("reason") or detail.get("message") or "").casefold()
        if code in EXTERNAL_LLM_CONFIG_ERROR_CODES:
            return "配置未完整：请检查 endpoint、model、API key，并确认已选择外部 LLM 后端。"
        if code == "external_llm_http_error":
            if status_int in {401, 403}:
                return "服务已连通，但鉴权失败：请检查 API key 是否正确、是否属于该服务商。"
            if status_int == 404:
                return "服务已连通，但 endpoint 或 model 可能不存在：请核对接口地址和模型名。"
            if status_int == 429:
                return "服务已连通，但触发限流或额度不足：稍后重试或检查账户额度。"
            if status_int and 500 <= status_int < 600:
                return "服务端返回 5xx：通常是上游服务暂时不可用，建议稍后重试。"
            return "服务返回 HTTP 错误：请查看状态码和响应预览。"
        if code == "external_llm_transport_error":
            if "timed out" in reason or "timeout" in reason:
                return "网络已开始连接但等待响应超时：可增大超时、换更快模型，或检查代理稳定性。"
            if "getaddrinfo" in reason or "name resolution" in reason or "nodename" in reason:
                return "DNS 解析失败：本机可能无法直接解析服务域名，建议填写代理 URL。"
            if "proxy" in reason or "tunnel" in reason:
                return "代理连接失败：请确认代理软件正在监听该端口，并且协议是 http://。"
            if "certificate" in reason or "ssl" in reason or "tls" in reason:
                return "TLS/证书握手失败：常见原因是代理证书、系统时间或 HTTPS 拦截配置异常。"
            return "传输层失败：请检查网络、代理地址和服务商域名是否可达。"
        if code in {"external_llm_invalid_response", "external_llm_empty_content"}:
            return "模型返回格式不符合预期：外部文字模式通常更宽容，JSON 模式需要模型严格输出结构化内容。"
        return "外部 LLM 自检失败：请查看错误代码、状态码和响应预览。"

    def llm_test(self, body: dict[str, Any]) -> dict[str, Any]:
        llm_config_override = self.request_scoped_llm_config(body.get("llm_config"))
        config = llm_config_override or self.llm_config
        request = {
            "llm_config_override": config.sanitized() if llm_config_override else None,
            "service_config": None if llm_config_override else config.sanitized(),
        }
        started = time.perf_counter()

        def elapsed_ms() -> int:
            return int((time.perf_counter() - started) * 1000)

        try:
            config.validate_ready()
        except LLMAdapterError as exc:
            issue = exc.as_issue()
            detail = issue.get("detail", {})
            return self.response_envelope(
                "/llm/test",
                request,
                {
                    "status": "blocked",
                    "llm_test": {
                        "ready": False,
                        "phase": "config",
                        "code": exc.code,
                        "message": str(exc),
                        "diagnostic": self.llm_test_diagnostic(exc.code, detail),
                        "elapsed_ms": elapsed_ms(),
                        "config": config.sanitized(),
                    },
                },
            )

        max_tokens = max(1, min(int(config.max_output_tokens or 16), 32))
        test_config = ExternalLLMConfig(
            enabled=True,
            provider=config.provider,
            endpoint=config.endpoint,
            model=config.model,
            api_key=config.api_key,
            timeout_seconds=config.timeout_seconds,
            max_output_tokens=max_tokens,
            proxy_url=config.proxy_url,
        )
        request_payload = {
            "model": test_config.model,
            "messages": [
                {
                    "role": "system",
                    "content": "You are a connection test. Reply with exactly OK.",
                },
                {"role": "user", "content": "Return exactly OK."},
            ],
            "temperature": 0,
            "max_tokens": max_tokens,
            "stream": False,
        }
        try:
            transport = self.llm_transport or default_openai_compatible_transport
            response = transport(test_config, request_payload)
            content = extract_chat_completion_content(response).strip()
        except LLMAdapterError as exc:
            issue = exc.as_issue()
            detail = issue.get("detail", {})
            status = detail.get("status")
            body_preview = str(detail.get("body") or "")[:500]
            return self.response_envelope(
                "/llm/test",
                request,
                {
                    "status": "failed",
                    "llm_test": {
                        "ready": False,
                        "phase": "transport" if exc.code in {"external_llm_http_error", "external_llm_transport_error"} else "response",
                        "code": exc.code,
                        "message": str(exc),
                        "http_status": status,
                        "body_preview": body_preview,
                        "diagnostic": self.llm_test_diagnostic(exc.code, detail),
                        "elapsed_ms": elapsed_ms(),
                        "config": test_config.sanitized(),
                    },
                },
            )
        except Exception as exc:  # pragma: no cover - defensive boundary for custom transports.
            code = "external_llm_unhandled_error"
            detail = {"message": str(exc), "type": type(exc).__name__}
            return self.response_envelope(
                "/llm/test",
                request,
                {
                    "status": "failed",
                    "llm_test": {
                        "ready": False,
                        "phase": "transport",
                        "code": code,
                        "message": str(exc),
                        "diagnostic": self.llm_test_diagnostic(code, detail),
                        "elapsed_ms": elapsed_ms(),
                        "config": test_config.sanitized(),
                    },
                },
            )

        return self.response_envelope(
            "/llm/test",
            request,
            {
                "status": "ok",
                "llm_test": {
                    "ready": True,
                    "phase": "response",
                    "code": "ok",
                    "message": "External LLM connection test succeeded.",
                    "diagnostic": "连接、鉴权和模型响应均正常；可以切回正式分析。",
                    "elapsed_ms": elapsed_ms(),
                    "response_preview": content[:160],
                    "usage": response.get("usage", {}) if isinstance(response, dict) else {},
                    "config": test_config.sanitized(),
                },
            },
        )

    def request_scoped_llm_config(self, raw_config: Any) -> ExternalLLMConfig | None:
        if not isinstance(raw_config, dict):
            return None
        has_value = any(str(raw_config.get(key) or "").strip() for key in ("endpoint", "model", "api_key", "provider"))
        if not has_value and raw_config.get("enabled") is None:
            return None

        def text_value(key: str, fallback: str = "") -> str:
            return str(raw_config.get(key) or fallback or "").strip()

        def float_value(key: str, fallback: float) -> float:
            try:
                return float(raw_config.get(key) or fallback)
            except (TypeError, ValueError):
                return fallback

        def int_value(key: str, fallback: int) -> int:
            try:
                return int(raw_config.get(key) or fallback)
            except (TypeError, ValueError):
                return fallback

        base = self.llm_config
        return ExternalLLMConfig(
            enabled=True,
            provider=text_value("provider", base.provider or "openai_compatible"),
            endpoint=text_value("endpoint", base.endpoint),
            model=text_value("model", base.model),
            api_key=str(raw_config.get("api_key") or base.api_key or ""),
            timeout_seconds=float_value("timeout_seconds", base.timeout_seconds),
            max_output_tokens=int_value("max_output_tokens", base.max_output_tokens),
            proxy_url=text_value("proxy_url", base.proxy_url),
        )

    def explain_analysis(
        self,
        records: list[Any],
        question: str = "",
        max_paths: int | None = None,
        max_hops: int | None = None,
        evidence_limit: int | None = None,
        subgraph_max_hops: int = 1,
        adapter_backend: str = "local",
        adapter_output: dict[str, Any] | None = None,
        include_input_pack: bool = False,
        context: Any = None,
        context_mode: str = "soft",
        llm_config_override: ExternalLLMConfig | None = None,
        analyzed_response: dict[str, Any] | None = None,
        analysis_endpoint: str = "/analyze/metabolites",
    ) -> dict[str, Any]:
        effective_llm_config = llm_config_override or self.llm_config
        request = {
            "records": records,
            "question": str(question or ""),
            "max_paths": max_paths,
            "max_hops": max_hops,
            "context": context,
            "context_mode": context_mode,
            "evidence_limit": evidence_limit,
            "subgraph_max_hops": subgraph_max_hops,
            "adapter_backend": adapter_backend,
            "include_input_pack": bool(include_input_pack),
            "adapter_output_supplied": adapter_output is not None,
            "llm_config_override": effective_llm_config.sanitized() if llm_config_override else None,
            "analyzed_response_supplied": analyzed_response is not None,
            "analysis_endpoint": analysis_endpoint,
        }
        input_pack = self.build_llm_adapter_input_pack(
            records,
            question=question,
            max_paths=max_paths,
            max_hops=max_hops,
            evidence_limit=evidence_limit,
            subgraph_max_hops=subgraph_max_hops,
            context=context,
            context_mode=context_mode,
            analyzed_response=analyzed_response,
            analysis_endpoint=analysis_endpoint,
        )
        if adapter_output is not None:
            raw_output = adapter_output
            guard = guard_adapter_output(input_pack, raw_output)
            explanation = self.finalize_guarded_adapter_output(input_pack, raw_output, guard)
            external_fallback_issue = None
        elif adapter_backend == "local":
            explanation = build_local_adapter_output(input_pack)
            guard = explanation.get("guard", {})
            if not guard.get("passed"):
                explanation = self.blocked_llm_adapter_output(explanation, guard, input_pack)
            external_fallback_issue = None
        elif adapter_backend in {"external_text", "external_llm_text", "text_llm"}:
            backend = ExternalLLMTextNarratorBackend(effective_llm_config, transport=self.llm_transport)
            backend_audit: dict[str, Any] = {"backend": "external_text", "config": effective_llm_config.sanitized()}
            external_fallback_issue = None
            try:
                backend_result = backend.generate(input_pack)
                raw_output = backend_result.output
                backend_audit = backend_result.audit
                guard = guard_adapter_output(input_pack, raw_output)
            except LLMAdapterError as exc:
                if exc.code not in EXTERNAL_LLM_CONFIG_ERROR_CODES:
                    explanation, external_fallback_issue = self.external_llm_fallback_output(
                        input_pack, exc, "external_text", effective_llm_config
                    )
                else:
                    raw_output = self.external_llm_error_output(exc)
                    guard = guard_adapter_output(input_pack, raw_output)
                    guard = self.force_guard_issue(guard, exc.as_issue())
                    explanation = self.finalize_guarded_adapter_output(input_pack, raw_output, guard, backend_audit=backend_audit)
            except Exception as exc:  # pragma: no cover - defensive boundary for external client code.
                explanation, external_fallback_issue = self.external_llm_fallback_output(
                    input_pack, exc, "external_text", effective_llm_config
                )
            else:
                explanation = self.finalize_guarded_adapter_output(input_pack, raw_output, guard, backend_audit=backend_audit)
        elif adapter_backend in {"external_text_chunked", "chunked_text", "hierarchical_text", "hierarchical_llm"}:
            backend = ExternalLLMChunkedTextNarratorBackend(effective_llm_config, transport=self.llm_transport)
            backend_audit = {"backend": "external_text_chunked", "config": effective_llm_config.sanitized()}
            external_fallback_issue = None
            try:
                backend_result = backend.generate(input_pack)
                raw_output = backend_result.output
                backend_audit = backend_result.audit
                guard = guard_adapter_output(input_pack, raw_output)
            except LLMAdapterError as exc:
                if exc.code not in EXTERNAL_LLM_CONFIG_ERROR_CODES:
                    explanation, external_fallback_issue = self.external_llm_fallback_output(
                        input_pack, exc, "external_text_chunked", effective_llm_config
                    )
                else:
                    raw_output = self.external_llm_error_output(exc)
                    guard = guard_adapter_output(input_pack, raw_output)
                    guard = self.force_guard_issue(guard, exc.as_issue())
                    explanation = self.finalize_guarded_adapter_output(input_pack, raw_output, guard, backend_audit=backend_audit)
            except Exception as exc:  # pragma: no cover - defensive boundary for external client code.
                explanation, external_fallback_issue = self.external_llm_fallback_output(
                    input_pack, exc, "external_text_chunked", effective_llm_config
                )
            else:
                explanation = self.finalize_guarded_adapter_output(input_pack, raw_output, guard, backend_audit=backend_audit)
        elif adapter_backend in {"external_llm", "llm", "openai_compatible"}:
            backend = ExternalLLMNarratorBackend(effective_llm_config, transport=self.llm_transport)
            backend_audit: dict[str, Any] = {"backend": "external_llm", "config": effective_llm_config.sanitized()}
            external_fallback_issue = None
            try:
                backend_result = backend.generate(input_pack)
                raw_output = backend_result.output
                backend_audit = backend_result.audit
                guard = guard_adapter_output(input_pack, raw_output)
            except LLMAdapterError as exc:
                if exc.code not in EXTERNAL_LLM_CONFIG_ERROR_CODES:
                    explanation, external_fallback_issue = self.external_llm_fallback_output(
                        input_pack, exc, "external_llm", effective_llm_config
                    )
                else:
                    raw_output = self.external_llm_error_output(exc)
                    guard = guard_adapter_output(input_pack, raw_output)
                    guard = self.force_guard_issue(guard, exc.as_issue())
                    explanation = self.finalize_guarded_adapter_output(input_pack, raw_output, guard, backend_audit=backend_audit)
            except Exception as exc:  # pragma: no cover - defensive boundary for external client code.
                explanation, external_fallback_issue = self.external_llm_fallback_output(
                    input_pack, exc, "external_llm", effective_llm_config
                )
            else:
                explanation = self.finalize_guarded_adapter_output(input_pack, raw_output, guard, backend_audit=backend_audit)
        else:
            external_fallback_issue = None
            raw_output = {
                "contract_version": LLM_ADAPTER_OUTPUT_CONTRACT_VERSION,
                "adapter_version": str(adapter_backend or "external"),
                "status": "blocked_by_guard",
                "insufficient_evidence": {
                    "reason": f"Adapter backend {adapter_backend!r} is not enabled in the local safety layer.",
                    "source_refs": [
                        {
                            "source_type": "release",
                            "ref_id": f"release:{self.release_id}",
                            "path": "$.release.release_id",
                        }
                    ],
                },
            }
            guard = guard_adapter_output(input_pack, raw_output)
            guard = self.force_guard_issue(
                guard,
                {
                    "code": "adapter_backend_not_enabled",
                    "path": "$.adapter_backend",
                    "detail": {"adapter_backend": adapter_backend},
                },
            )
            explanation = self.blocked_llm_adapter_output(raw_output, guard, input_pack)

        payload = {
            "contract_version": EXPLAIN_CONTRACT_VERSION,
            "status": explanation.get("status", "blocked_by_guard"),
            "adapter_backend": str(adapter_backend or "local") if adapter_output is None else "supplied",
            "adapter_version": explanation.get("adapter_version", LOCAL_LLM_ADAPTER_VERSION),
            "adapter_input_contract_version": input_pack.get("contract_version", ""),
            "adapter_input_hash": content_hash(input_pack),
            "guard_passed": bool(explanation.get("guard", {}).get("passed")),
            "guard_issue_count": int(explanation.get("guard", {}).get("issue_count", 0) or 0),
            "explanation": explanation,
        }
        if external_fallback_issue:
            payload["external_fallback_issue"] = external_fallback_issue
        if include_input_pack:
            payload["adapter_input_pack"] = input_pack
        return self.response_envelope("/explain", request, payload)

    def chat(self, body: dict[str, Any]) -> dict[str, Any]:
        records, normalization = chat_records_from_body(body)
        question = chat_question_from_body(body)
        adapter_backend = str(body.get("adapter_backend") or "local")
        input_mode = str(body.get("input_mode") or "metabolite_table").strip() or "metabolite_table"
        max_paths = body.get("max_paths")
        max_hops = body.get("max_hops")
        evidence_limit = body.get("evidence_limit")
        subgraph_max_hops = int(body.get("subgraph_max_hops") or 1)
        context = body.get("context")
        context_mode = str(body.get("context_mode") or "soft")
        llm_config_override = self.request_scoped_llm_config(body.get("llm_config"))
        request = {
            "question": question,
            "table_text_sha256": hashlib.sha256(str(body.get("table_text") or "").encode("utf-8")).hexdigest()
            if body.get("table_text")
            else "",
            "has_records": bool(body.get("records")),
            "record_count": len(records),
            "input_mode": input_mode,
            "adapter_backend": adapter_backend,
            "max_paths": max_paths,
            "max_hops": max_hops,
            "context": context,
            "context_mode": context_mode,
            "evidence_limit": evidence_limit,
            "subgraph_max_hops": subgraph_max_hops,
            "llm_config_override": llm_config_override.sanitized() if llm_config_override else None,
        }
        if not records:
            payload = {
                "contract_version": CHAT_CONTRACT_VERSION,
                "status": "needs_records",
                "question": question,
                "records": [],
                "input_normalization": normalization,
                "message": "Upload a metabolite table, paste a table, or provide a metabolite list before asking for an explanation.",
            }
            return self.response_envelope("/chat", request, payload)

        if input_mode in {"differential_table", "trait_score_differential", "trait_score_diff", "traitscore_differential"}:
            q_threshold = body.get("q_threshold", 0.05)
            p_threshold = body.get("p_threshold", None)
            analyzed = self.analyze_differential_table(
                records,
                max_paths=max_paths,
                max_hops=max_hops,
                context=context,
                context_mode=context_mode,
                q_threshold=None if q_threshold in {None, ""} else float(q_threshold),
                p_threshold=None if p_threshold in {None, ""} else float(p_threshold),
                min_abs_effect=float(body.get("min_abs_effect") or 0.0),
                top_per_group=int(body.get("top_per_group") or 120),
                max_records=int(body.get("max_records") or 300),
            )
            selected_records = list(analyzed.get("selected_records") or [])
            effective_context = context or infer_differential_table_context(selected_records)
            if not selected_records:
                payload = {
                    "contract_version": CHAT_CONTRACT_VERSION,
                    "status": "needs_records",
                    "question": question,
                    "input_mode": input_mode,
                    "records": [],
                    "selected_records": [],
                    "input_normalization": normalization,
                    "differential_table_selection": analyzed.get("differential_table_selection", {}),
                    "analysis_pack": analyzed.get("analysis_pack", {}),
                    "precheck": analyzed.get("precheck", {}),
                    "message": "No differential table rows passed the current selection thresholds; relax q/p/effect filters or provide eligible rows.",
                }
                return self.response_envelope("/chat", request, payload)

            explained = self.explain_analysis(
                selected_records,
                question=question,
                max_paths=max_paths,
                max_hops=max_hops,
                evidence_limit=evidence_limit,
                subgraph_max_hops=subgraph_max_hops,
                adapter_backend=adapter_backend,
                include_input_pack=True,
                context=effective_context,
                context_mode=context_mode,
                llm_config_override=llm_config_override,
                analyzed_response=analyzed,
                analysis_endpoint="/analyze/differential-table",
            )
            input_pack = explained.pop("adapter_input_pack", {})
            analysis_pack = input_pack.get("analysis_pack", {}) if isinstance(input_pack, dict) else {}
            prediction_pack = input_pack.get("prediction_pack", {}) if isinstance(input_pack, dict) else {}
            interpretation_report = (
                analysis_pack.get("interpretation_report")
                if isinstance(analysis_pack, dict)
                else {}
            ) or analyzed.get("interpretation_report", {}) or self.build_interpretation_report(
                analysis_pack if isinstance(analysis_pack, dict) else {},
                prediction_pack if isinstance(prediction_pack, dict) else {},
                question=question,
                context=effective_context,
            )
            payload = {
                "contract_version": CHAT_CONTRACT_VERSION,
                "status": explained.get("status", "blocked_by_guard"),
                "question": question,
                "input_mode": input_mode,
                "records": selected_records,
                "selected_records": selected_records,
                "input_normalization": normalization,
                "differential_table_selection": analyzed.get("differential_table_selection", {}),
                "analysis_pack": analysis_pack,
                "structured_prediction": analyzed.get("structured_prediction", {}),
                "precheck": analyzed.get("precheck", {}),
                "analysis_features": analyzed.get("analysis_features", {}),
                "rankings": analyzed.get("rankings", {}),
                "directional_enrichment": analyzed.get("directional_enrichment", {}),
                "propagation": analyzed.get("propagation", {}),
                "score_notes": analyzed.get("score_notes", {}),
                "predictions": prediction_pack,
                "interpretation_report": interpretation_report,
                "source_analysis_response_hash": (analyzed.get("determinism") or {}).get("response_hash", ""),
                "explanation": explained.get("explanation", {}),
                "explain_response": {
                    key: value
                    for key, value in explained.items()
                    if key not in {"api_version", "endpoint", "release", "determinism", "explanation"}
                },
            }
            payload["explain_response"]["adapter_input_source_endpoints"] = input_pack.get("source_endpoints", [])
            return self.response_envelope("/chat", request, payload)

        explained = self.explain_analysis(
            records,
            question=question,
            max_paths=max_paths,
            max_hops=max_hops,
            evidence_limit=evidence_limit,
            subgraph_max_hops=subgraph_max_hops,
            adapter_backend=adapter_backend,
            include_input_pack=True,
            context=context,
            context_mode=context_mode,
            llm_config_override=llm_config_override,
        )
        input_pack = explained.pop("adapter_input_pack", {})
        analysis_pack = input_pack.get("analysis_pack", {}) if isinstance(input_pack, dict) else {}
        prediction_pack = input_pack.get("prediction_pack", {}) if isinstance(input_pack, dict) else {}
        interpretation_report = (
            analysis_pack.get("interpretation_report")
            if isinstance(analysis_pack, dict)
            else {}
        ) or self.build_interpretation_report(
            analysis_pack if isinstance(analysis_pack, dict) else {},
            prediction_pack if isinstance(prediction_pack, dict) else {},
            question=question,
            context=context,
        )
        payload = {
            "contract_version": CHAT_CONTRACT_VERSION,
            "status": explained.get("status", "blocked_by_guard"),
            "question": question,
            "records": records,
            "input_normalization": normalization,
            "analysis_pack": analysis_pack,
            "predictions": prediction_pack,
            "interpretation_report": interpretation_report,
            "explanation": explained.get("explanation", {}),
            "explain_response": {
                key: value
                for key, value in explained.items()
                if key not in {"api_version", "endpoint", "release", "determinism", "explanation"}
            },
        }
        payload["explain_response"]["adapter_input_source_endpoints"] = input_pack.get("source_endpoints", [])
        return self.response_envelope("/chat", request, payload)


class MetaboRequestHandler(BaseHTTPRequestHandler):
    service: MetaboService

    def log_message(self, format: str, *args: Any) -> None:  # noqa: A002
        sys.stderr.write("%s - - [%s] %s\n" % (self.address_string(), self.log_date_time_string(), format % args))

    def read_json_body(self) -> dict[str, Any]:
        length = int(self.headers.get("Content-Length", "0") or 0)
        raw = self.rfile.read(length).decode("utf-8") if length else "{}"
        if not raw.strip():
            return {}
        return json.loads(raw)

    def send_json(self, status: int, body: dict[str, Any]) -> None:
        payload = json.dumps(body, ensure_ascii=False, sort_keys=True).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(payload)))
        self.end_headers()
        self.wfile.write(payload)

    def send_html(self, status: int, text: str) -> None:
        payload = text.encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Content-Length", str(len(payload)))
        self.end_headers()
        self.wfile.write(payload)

    def chat_ui_html(self) -> str:
        path = self.service.workspace / "web" / "chat.html"
        if path.exists():
            return path.read_text(encoding="utf-8")
        return "<!doctype html><html><body><h1>Metabo Chat</h1><p>web/chat.html is missing.</p></body></html>"

    def do_GET(self) -> None:  # noqa: N802
        try:
            parsed = urlparse(self.path)
            path = parsed.path.rstrip("/") or "/"
            query = parse_qs(parsed.query)
            if path in {"/", "/chat", "/flow-test"}:
                self.send_html(200, self.chat_ui_html())
            elif path == "/releases":
                self.send_json(200, self.service.releases())
            elif path == "/llm/status":
                self.send_json(200, self.service.llm_status())
            elif path.startswith("/entity/"):
                self.send_json(200, self.service.entity_detail(unquote(path.split("/entity/", 1)[1])))
            elif path.startswith("/pubchem/"):
                self.send_json(200, self.service.pubchem(unquote(path.split("/pubchem/", 1)[1])))
            elif path == "/subgraph":
                seed_ids = query.get("seed_id", []) + query.get("seed_ids", [])
                if len(seed_ids) == 1 and "," in seed_ids[0]:
                    seed_ids = [item for item in seed_ids[0].split(",") if item]
                max_hops = int(query.get("max_hops", [1])[0])
                self.send_json(200, self.service.subgraph(seed_ids, max_hops=max_hops))
            elif path == "/evidence":
                self.send_json(
                    200,
                    self.service.evidence(
                        edge_uid=query.get("edge_uid", [""])[0],
                        subject_uid=query.get("subject_uid", [""])[0],
                        predicate=query.get("predicate", [""])[0],
                        object_uid=query.get("object_uid", [""])[0],
                        entity_uid=query.get("entity_uid", [""])[0],
                        sentence_uid=query.get("sentence_uid", [""])[0],
                        relation_uid=query.get("relation_uid", [""])[0],
                        limit=int(query.get("limit", [self.service.config.max_evidence_items])[0]),
                    ),
                )
            else:
                if path == "/llm/status":
                    self.send_json(200, self.service.llm_status())
                else:
                    self.send_json(404, {"error": "not_found", "path": path})
        except Exception as exc:  # pragma: no cover - exercised manually.
            self.send_json(500, {"error": type(exc).__name__, "message": str(exc)})

    def do_POST(self) -> None:  # noqa: N802
        try:
            parsed = urlparse(self.path)
            path = parsed.path.rstrip("/") or "/"
            body = self.read_json_body()
            if path == "/resolve":
                self.send_json(
                    200,
                    self.service.resolve(
                        str(body.get("query", "")),
                        entity_type=body.get("entity_type"),
                        namespace=body.get("namespace"),
                        limit=int(body.get("limit") or self.service.config.max_candidates),
                    ),
                )
            elif path == "/precheck/metabolites":
                if body.get("table_text") or body.get("table") or body.get("csv") or body.get("tsv"):
                    records, normalization = prepare_analysis_records(
                        *parse_chat_table_text(body.get("table_text") or body.get("table") or body.get("csv") or body.get("tsv"))
                    )
                    self.send_json(200, self.service.precheck_metabolites(records, input_normalization=normalization))
                else:
                    self.send_json(200, self.service.precheck_metabolites(list(body.get("records") or [])))
            elif path == "/analyze/metabolites":
                records = list(body.get("records") or [])
                if not records and (body.get("table_text") or body.get("table") or body.get("csv") or body.get("tsv")):
                    records, _normalization = prepare_analysis_records(
                        *parse_chat_table_text(body.get("table_text") or body.get("table") or body.get("csv") or body.get("tsv"))
                    )
                self.send_json(
                    200,
                    self.service.analyze_metabolites(
                        records,
                        max_paths=body.get("max_paths"),
                        max_hops=body.get("max_hops"),
                        context=body.get("context"),
                        context_mode=str(body.get("context_mode") or "soft"),
                    ),
                )
            elif path == "/analyze/trait-score":
                records = list(body.get("records") or [])
                if not records and (body.get("table_text") or body.get("table") or body.get("csv") or body.get("tsv")):
                    records, _normalization = parse_chat_table_text(
                        body.get("table_text") or body.get("table") or body.get("csv") or body.get("tsv")
                    )
                q_threshold = body.get("q_threshold", 0.05)
                self.send_json(
                    200,
                    self.service.analyze_trait_score(
                        records,
                        max_paths=body.get("max_paths"),
                        max_hops=body.get("max_hops"),
                        context=body.get("context"),
                        context_mode=str(body.get("context_mode") or "soft"),
                        q_threshold=None if q_threshold in {None, ""} else float(q_threshold),
                        min_abs_effect=float(body.get("min_abs_effect") or 0.0),
                        top_per_celltype=int(body.get("top_per_celltype") or 80),
                        max_records=int(body.get("max_records") or 120),
                    ),
                )
            elif path == "/analyze/differential-table":
                records = list(body.get("records") or [])
                if not records and (body.get("table_text") or body.get("table") or body.get("csv") or body.get("tsv")):
                    records, _normalization = parse_chat_table_text(
                        body.get("table_text") or body.get("table") or body.get("csv") or body.get("tsv")
                    )
                q_threshold = body.get("q_threshold", 0.05)
                p_threshold = body.get("p_threshold", None)
                self.send_json(
                    200,
                    self.service.analyze_differential_table(
                        records,
                        max_paths=body.get("max_paths"),
                        max_hops=body.get("max_hops"),
                        context=body.get("context"),
                        context_mode=str(body.get("context_mode") or "soft"),
                        q_threshold=None if q_threshold in {None, ""} else float(q_threshold),
                        p_threshold=None if p_threshold in {None, ""} else float(p_threshold),
                        min_abs_effect=float(body.get("min_abs_effect") or 0.0),
                        top_per_group=int(body.get("top_per_group") or 120),
                        max_records=int(body.get("max_records") or 300),
                    ),
                )
            elif path == "/explain":
                records = list(body.get("records") or [])
                if not records and (body.get("table_text") or body.get("table") or body.get("csv") or body.get("tsv")):
                    records, _normalization = prepare_analysis_records(
                        *parse_chat_table_text(body.get("table_text") or body.get("table") or body.get("csv") or body.get("tsv"))
                    )
                self.send_json(
                    200,
                    self.service.explain_analysis(
                        records,
                        question=str(body.get("question", "") or ""),
                        max_paths=body.get("max_paths"),
                        max_hops=body.get("max_hops"),
                        evidence_limit=body.get("evidence_limit"),
                        subgraph_max_hops=int(body.get("subgraph_max_hops") or 1),
                        adapter_backend=str(body.get("adapter_backend", "local") or "local"),
                        include_input_pack=bool(body.get("include_input_pack", False)),
                        context=body.get("context"),
                        context_mode=str(body.get("context_mode") or "soft"),
                        llm_config_override=self.service.request_scoped_llm_config(body.get("llm_config")),
                    ),
                )
            elif path == "/chat":
                self.send_json(200, self.service.chat(body))
            elif path == "/llm/test":
                self.send_json(200, self.service.llm_test(body))
            elif path == "/subgraph":
                self.send_json(
                    200,
                    self.service.subgraph(
                        list(body.get("seed_ids") or []),
                        max_hops=int(body.get("max_hops") or 1),
                        edge_types=body.get("edge_types"),
                    ),
                )
            elif path == "/evidence":
                self.send_json(
                    200,
                    self.service.evidence(
                        edge_uid=str(body.get("edge_uid", "") or ""),
                        subject_uid=str(body.get("subject_uid", "") or ""),
                        predicate=str(body.get("predicate", "") or ""),
                        object_uid=str(body.get("object_uid", "") or ""),
                        entity_uid=str(body.get("entity_uid", "") or ""),
                        sentence_uid=str(body.get("sentence_uid", "") or ""),
                        relation_uid=str(body.get("relation_uid", "") or ""),
                        limit=int(body.get("limit") or self.service.config.max_evidence_items),
                    ),
                )
            else:
                self.send_json(404, {"error": "not_found", "path": path})
        except Exception as exc:  # pragma: no cover - exercised manually.
            self.send_json(500, {"error": type(exc).__name__, "message": str(exc)})


def load_records_arg(path_or_json: str) -> list[Any]:
    stripped = path_or_json.strip()
    if stripped.startswith(("[", "{")):
        text = stripped
        parsed = json.loads(text)
        if isinstance(parsed, list):
            return parsed
        if isinstance(parsed, dict) and isinstance(parsed.get("records"), list):
            return parsed["records"]
        raise ValueError("Expected a JSON list or an object with a records list.")
    else:
        path = Path(path_or_json)
        if not path.exists():
            raise FileNotFoundError(f"Records file does not exist: {path}")
        if path.suffix.casefold() in {".csv", ".tsv", ".txt", ".json", ".jsonl", ".gz"}:
            return load_records_file(path)
        text = path.read_text(encoding="utf-8-sig")
    parsed = json.loads(text)
    if isinstance(parsed, list):
        return parsed
    if isinstance(parsed, dict) and isinstance(parsed.get("records"), list):
        return parsed["records"]
    raise ValueError("Expected a JSON/JSONL/CSV/TSV records file, a JSON list, or an object with a records list.")


def load_context_arg(path_or_json: str) -> Any:
    stripped = str(path_or_json or "").strip()
    if not stripped:
        return None
    if stripped.startswith(("{", "[")):
        return json.loads(stripped)
    path = Path(stripped)
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8-sig"))
    return {"context_terms": stripped}


def compact_precheck_summary(precheck: dict[str, Any]) -> dict[str, Any]:
    input_normalization = precheck.get("input_normalization") or {}
    expanded_summary = precheck.get("expanded_summary") or {}
    return {
        "endpoint": precheck.get("endpoint", "/precheck/metabolites"),
        "analysis_mode": precheck.get("analysis_mode", input_normalization.get("analysis_mode", "metabolite_table")),
        "input_count": precheck.get("input_count", 0),
        "summary": precheck.get("summary", {}),
        "expanded_summary": {
            key: expanded_summary.get(key)
            for key in (
                "strict_matched_count",
                "expanded_candidate_count",
                "expanded_candidate_input_count",
                "unresolved_input_count",
                "analysis_seed_input_count",
                "weighted_seed_mass",
            )
            if key in expanded_summary
        },
        "input_normalization": {
            key: input_normalization.get(key)
            for key in (
                "source",
                "record_count",
                "analysis_mode",
                "input_format_detected",
                "converted_record_count",
                "comparison_primary_group_counts",
            )
            if key in input_normalization
        },
        "warning_codes": [row.get("code", "") for row in input_normalization.get("warnings", []) or []],
        "warnings": input_normalization.get("warnings", []) or [],
        "supported_identifier_columns": precheck.get("supported_identifier_columns", []),
    }


def build_service_from_args(args: argparse.Namespace) -> MetaboService:
    workspace = Path(args.workspace).resolve()
    return MetaboService(
        workspace=workspace,
        normalized_root=(workspace / args.normalized_root).resolve(),
        graph_root=(workspace / args.graph_root).resolve(),
        pubchem_root=(workspace / args.pubchem_root).resolve(),
        compound_root=(workspace / args.compound_root).resolve(),
        literature_root=(workspace / args.literature_root).resolve(),
        prediction_root=(workspace / args.prediction_root).resolve(),
        release_id=args.release_id or None,
        llm_config=ExternalLLMConfig.from_env(
            enabled=True if getattr(args, "enable_external_llm", False) else None,
            provider=getattr(args, "llm_provider", "") or "",
            endpoint=getattr(args, "llm_endpoint", "") or "",
            model=getattr(args, "llm_model", "") or "",
            api_key_env=getattr(args, "llm_api_key_env", "") or "",
            timeout_seconds=getattr(args, "llm_timeout_seconds", None) or None,
            max_output_tokens=getattr(args, "llm_max_output_tokens", None) or None,
            proxy_url=getattr(args, "llm_proxy", "") or "",
        ),
    )


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Read-only MVP API/service layer for frozen metabolism releases.")
    parser.add_argument("--workspace", default=".")
    parser.add_argument("--normalized-root", default=DEFAULT_NORMALIZED_ROOT)
    parser.add_argument("--graph-root", default=DEFAULT_GRAPH_ROOT)
    parser.add_argument("--pubchem-root", default=DEFAULT_PUBCHEM_ROOT)
    parser.add_argument("--compound-root", default=DEFAULT_COMPOUND_MATCH_ROOT)
    parser.add_argument("--literature-root", default=DEFAULT_LITERATURE_ROOT)
    parser.add_argument("--prediction-root", default=DEFAULT_PREDICTION_OVERLAY_ROOT)
    parser.add_argument("--release-id", default="")
    parser.add_argument("--enable-external-llm", action="store_true")
    parser.add_argument("--llm-provider", default="")
    parser.add_argument("--llm-endpoint", default="")
    parser.add_argument("--llm-model", default="")
    parser.add_argument("--llm-api-key-env", default="")
    parser.add_argument("--llm-timeout-seconds", type=float, default=None)
    parser.add_argument("--llm-max-output-tokens", type=int, default=None)
    parser.add_argument("--llm-proxy", default="", help="Optional HTTP proxy URL for external LLM requests, e.g. http://127.0.0.1:7890")
    subparsers = parser.add_subparsers(dest="command", required=True)

    resolve_parser = subparsers.add_parser("resolve")
    resolve_parser.add_argument("query")
    resolve_parser.add_argument("--entity-type", default="")
    resolve_parser.add_argument("--namespace", default="")

    precheck_parser = subparsers.add_parser("precheck-metabolites")
    precheck_parser.add_argument("records_json", metavar="records")
    precheck_parser.add_argument("--summary", action="store_true", help="Print compact mode/count/warning summary instead of full row-level JSON.")

    analyze_parser = subparsers.add_parser("analyze-metabolites")
    analyze_parser.add_argument("records_json", metavar="records")
    analyze_parser.add_argument("--max-paths", type=int, default=DEFAULT_MAX_PATHS)
    analyze_parser.add_argument("--max-hops", type=int, default=DEFAULT_MAX_HOPS)
    analyze_parser.add_argument("--context-json", default="")
    analyze_parser.add_argument("--context-mode", default="soft", choices=["soft", "hard"])

    trait_score_parser = subparsers.add_parser("analyze-trait-score")
    trait_score_parser.add_argument("records_json", metavar="records")
    trait_score_parser.add_argument("--max-paths", type=int, default=DEFAULT_MAX_PATHS)
    trait_score_parser.add_argument("--max-hops", type=int, default=DEFAULT_MAX_HOPS)
    trait_score_parser.add_argument("--context-json", default="")
    trait_score_parser.add_argument("--context-mode", default="soft", choices=["soft", "hard"])
    trait_score_parser.add_argument("--q-threshold", type=float, default=0.05)
    trait_score_parser.add_argument("--no-q-threshold", action="store_true")
    trait_score_parser.add_argument("--min-abs-effect", type=float, default=0.0)
    trait_score_parser.add_argument("--top-per-celltype", type=int, default=80)
    trait_score_parser.add_argument("--max-records", type=int, default=120)

    differential_parser = subparsers.add_parser("analyze-differential-table")
    differential_parser.add_argument("records_json", metavar="records")
    differential_parser.add_argument("--max-paths", type=int, default=DEFAULT_MAX_PATHS)
    differential_parser.add_argument("--max-hops", type=int, default=DEFAULT_MAX_HOPS)
    differential_parser.add_argument("--context-json", default="")
    differential_parser.add_argument("--context-mode", default="soft", choices=["soft", "hard"])
    differential_parser.add_argument("--q-threshold", type=float, default=0.05)
    differential_parser.add_argument("--no-q-threshold", action="store_true")
    differential_parser.add_argument("--p-threshold", type=float, default=None)
    differential_parser.add_argument("--min-abs-effect", type=float, default=0.0)
    differential_parser.add_argument("--top-per-group", type=int, default=120)
    differential_parser.add_argument("--max-records", type=int, default=300)

    predict_parser = subparsers.add_parser("predict-metabolites")
    predict_parser.add_argument("records_json", metavar="records")
    predict_parser.add_argument("--max-paths", type=int, default=DEFAULT_MAX_PATHS)
    predict_parser.add_argument("--max-hops", type=int, default=DEFAULT_MAX_HOPS)
    predict_parser.add_argument("--context-json", default="")
    predict_parser.add_argument("--context-mode", default="soft", choices=["soft", "hard"])

    explain_parser = subparsers.add_parser("explain")
    explain_parser.add_argument("records_json", metavar="records")
    explain_parser.add_argument("--question", default="")
    explain_parser.add_argument("--max-paths", type=int, default=DEFAULT_MAX_PATHS)
    explain_parser.add_argument("--max-hops", type=int, default=DEFAULT_MAX_HOPS)
    explain_parser.add_argument("--context-json", default="")
    explain_parser.add_argument("--context-mode", default="soft", choices=["soft", "hard"])
    explain_parser.add_argument("--evidence-limit", type=int, default=DEFAULT_MAX_EVIDENCE_ITEMS)
    explain_parser.add_argument("--subgraph-max-hops", type=int, default=1)
    explain_parser.add_argument(
        "--adapter-backend",
        default="local",
        choices=[
            "local",
            "external_text",
            "external_text_chunked",
            "chunked_text",
            "hierarchical_text",
            "hierarchical_llm",
            "external_llm_text",
            "text_llm",
            "external_llm",
            "llm",
            "openai_compatible",
        ],
    )
    explain_parser.add_argument("--include-input-pack", action="store_true")

    entity_parser = subparsers.add_parser("entity")
    entity_parser.add_argument("entity_id")

    pubchem_parser = subparsers.add_parser("pubchem")
    pubchem_parser.add_argument("cid")

    subgraph_parser = subparsers.add_parser("subgraph")
    subgraph_parser.add_argument("seed_ids", nargs="+")
    subgraph_parser.add_argument("--max-hops", type=int, default=1)

    evidence_parser = subparsers.add_parser("evidence")
    evidence_parser.add_argument("--edge-uid", default="")
    evidence_parser.add_argument("--subject-uid", default="")
    evidence_parser.add_argument("--predicate", default="")
    evidence_parser.add_argument("--object-uid", default="")
    evidence_parser.add_argument("--entity-uid", default="")
    evidence_parser.add_argument("--sentence-uid", default="")
    evidence_parser.add_argument("--relation-uid", default="")
    evidence_parser.add_argument("--limit", type=int, default=DEFAULT_MAX_EVIDENCE_ITEMS)

    subparsers.add_parser("releases")

    serve_parser = subparsers.add_parser("serve")
    serve_parser.add_argument("--host", default="127.0.0.1")
    serve_parser.add_argument("--port", type=int, default=8765)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    args = parse_args(argv or sys.argv[1:])
    service = build_service_from_args(args)
    if args.command == "resolve":
        output = service.resolve(args.query, entity_type=args.entity_type or None, namespace=args.namespace or None)
    elif args.command == "precheck-metabolites":
        output = service.precheck_metabolites(load_records_arg(args.records_json))
        if args.summary:
            output = compact_precheck_summary(output)
    elif args.command == "analyze-metabolites":
        output = service.analyze_metabolites(
            load_records_arg(args.records_json),
            max_paths=args.max_paths,
            max_hops=args.max_hops,
            context=load_context_arg(args.context_json),
            context_mode=args.context_mode,
        )
    elif args.command == "analyze-trait-score":
        output = service.analyze_trait_score(
            load_records_arg(args.records_json),
            max_paths=args.max_paths,
            max_hops=args.max_hops,
            context=load_context_arg(args.context_json),
            context_mode=args.context_mode,
            q_threshold=None if args.no_q_threshold else args.q_threshold,
            min_abs_effect=args.min_abs_effect,
            top_per_celltype=args.top_per_celltype,
            max_records=args.max_records,
        )
    elif args.command == "analyze-differential-table":
        output = service.analyze_differential_table(
            load_records_arg(args.records_json),
            max_paths=args.max_paths,
            max_hops=args.max_hops,
            context=load_context_arg(args.context_json),
            context_mode=args.context_mode,
            q_threshold=None if args.no_q_threshold else args.q_threshold,
            p_threshold=args.p_threshold,
            min_abs_effect=args.min_abs_effect,
            top_per_group=args.top_per_group,
            max_records=args.max_records,
        )
    elif args.command == "predict-metabolites":
        analyzed = service.analyze_metabolites(
            load_records_arg(args.records_json),
            max_paths=args.max_paths,
            max_hops=args.max_hops,
            context=load_context_arg(args.context_json),
            context_mode=args.context_mode,
        )
        output = service.response_envelope(
            "/predict/metabolites",
            {
                "records": load_records_arg(args.records_json),
                "max_paths": args.max_paths,
                "max_hops": args.max_hops,
                "context": load_context_arg(args.context_json),
                "context_mode": args.context_mode,
            },
            {"predictions": analyzed.get("predictions", {})},
        )
    elif args.command == "explain":
        output = service.explain_analysis(
            load_records_arg(args.records_json),
            question=args.question,
            max_paths=args.max_paths,
            max_hops=args.max_hops,
            evidence_limit=args.evidence_limit,
            subgraph_max_hops=args.subgraph_max_hops,
            adapter_backend=args.adapter_backend,
            include_input_pack=args.include_input_pack,
            context=load_context_arg(args.context_json),
            context_mode=args.context_mode,
        )
    elif args.command == "entity":
        output = service.entity_detail(args.entity_id)
    elif args.command == "pubchem":
        output = service.pubchem(args.cid)
    elif args.command == "subgraph":
        output = service.subgraph(args.seed_ids, max_hops=args.max_hops)
    elif args.command == "evidence":
        output = service.evidence(
            edge_uid=args.edge_uid,
            subject_uid=args.subject_uid,
            predicate=args.predicate,
            object_uid=args.object_uid,
            entity_uid=args.entity_uid,
            sentence_uid=args.sentence_uid,
            relation_uid=args.relation_uid,
            limit=args.limit,
        )
    elif args.command == "releases":
        output = service.releases()
    elif args.command == "serve":
        MetaboRequestHandler.service = service
        server = ThreadingHTTPServer((args.host, args.port), MetaboRequestHandler)
        print(f"Serving {service.release_id} at http://{args.host}:{args.port}", file=sys.stderr)
        server.serve_forever()
        return 0
    else:  # pragma: no cover
        raise ValueError(f"Unsupported command: {args.command}")
    print(json.dumps(output, ensure_ascii=False, sort_keys=True, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
