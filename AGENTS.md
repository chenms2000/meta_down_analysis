# Project Goal

When working in this repository, keep the default goal aligned with the Codex
`/goal` objective:

Build the tumor metabolism knowledge base and output layer into an
evidence-grounded research explanation system. Organize metabolite matching,
knowledge-graph relationships, literature evidence, ranking results, and
confidence tiers into conclusions that users can audit.

Every conclusion should:

- Trace back to input metabolites, graph edges, literature evidence, ranking
  fields, or explicit model outputs.
- State confidence, downgrade reasons, and interpretation boundaries.
- Treat ambiguous matches, weak evidence, graph propagation, and model-only
  rankings as lower-confidence signals.
- Provide research interpretation and validation priorities only; do not make
  clinical decisions or treatment recommendations.
