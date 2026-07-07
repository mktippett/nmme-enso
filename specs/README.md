# specs/ — Behavioral Specifications

## What these files are

Each Markdown file documents the **external behavior** of one analysis script:
what it reads, what it writes, and precisely how it computes its outputs.  The
goal is that a collaborator (or yourself six months later) can understand any
figure or metric without reading the source code, and that the specs are
accurate enough to serve as a draft Methods section in a paper.

| Spec file | Script |
|-----------|--------|
| `latest_forecast.md` | `scripts/latest_forecast.py` |

## Structure of each spec

Every spec follows the same 7-section template:

1. **Purpose** — 1–2 sentence scientific purpose
2. **Inputs** — table of files, relevant columns, and filters applied
3. **Outputs** — table of files/plots and what each contains
4. **Algorithm** — numbered prose + formulae; enough detail to re-implement
5. **Constants & Scientific Rationale** — table of name / value / why
6. **Edge Cases & Error Handling** — non-obvious decisions and failure modes
7. **Synchronization Log** — date / code change / spec updated

Each spec also ends with a **Verification snippet**: a short runnable Python
block that asserts key invariants (no NaN in required columns, row counts in
expected range). Run it after changes as a smoke test.

## How to keep specs in sync

> **Rule:** when you change a constant, threshold, formula, or algorithm in a
> script, update the corresponding spec the **same session** and add a row to
> the Synchronization Log.

This is discipline-based, not automated.  The `> Last reviewed` badge at the
top of each spec shows the last date the spec was manually verified against
the code.

No tooling required — the cost of divergence is a wrong Methods section, which
is much worse than spending five minutes updating a Markdown file.
