# <script_name>.py — Behavioral Specification

> Last reviewed against code: YYYY-MM-DD

## Purpose

<!-- 1–2 sentences: what scientific question does this script address,
     what does it produce, and why. -->

## Inputs

| File | Relevant columns | Filters applied |
|------|-----------------|-----------------|
| `<path/to/data.parquet>` | all columns | year ∈ [Y1, Y2]; `filter_data()` |
| `observations/<obs.csv>` | `col_a, col_b` | field == value; no NaN in col_a |

## Outputs

| File | Contents | Format |
|------|----------|--------|
| `plots/01_example.pdf` | Description of what the figure shows | PDF, dpi=150 |
| `output.parquet` | Description of the intermediate table | Parquet, snappy |

## Algorithm

### 1. Load and filter

Call `load_data()` from `config.py`.  <!-- describe any additional selection. -->

### 2. <Step name>

<!-- Numbered prose. Include formulae where relevant, e.g.:

    score = Σ(f_i − o_i)² / Σ(clim − o_i)²

Enough detail that a reader could re-implement without reading the code.
-->

### 3. <Step name>

<!-- ... -->

## Constants & Scientific Rationale

| Name | Value | Rationale |
|------|-------|-----------|
| `YEARS` | range(Y1, Y2+1) | ... |
| `<CONST>` | `<value>` | citation or explanation |

## Edge Cases & Error Handling

- **<Situation>**: <how it is handled and why>
- **<Situation>**: <how it is handled and why>

## Verification Snippet

```python
# Run after changes to confirm key invariants
import pandas as pd
df = pd.read_parquet("output.parquet")

assert df["required_col"].notna().all(), "NaN in required_col"
assert 1_000 < len(df) < 10_000_000, f"Unexpected row count: {len(df)}"
# add project-specific assertions here
print("Verification passed.")
```

## Synchronization Log

| Date | Code change | Spec updated |
|------|-------------|--------------|
| YYYY-MM-DD | Initial spec written | ✓ |
