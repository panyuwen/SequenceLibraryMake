# Index Balance Tool

## Overview
This tool is designed to ensure robust index balance across sequencing lanes, with ratio-aware optimization.
This script performs index balance checking and **lane assignment optimization** for Illumina sequencing. It distributes **all samples across a specified number of lanes** and searches for an assignment that **meets C1/C2 feasibility thresholds and improves overall balance** (≥C3), taking **per‑sample `ratio` weights** into account when provided.

## Key Features
The main capabilities include:
- **Full-lane assignment**: assigns *all* input samples to the requested number of lanes, targeting index balance by optimizing a feasibility‑first objective.
- **Ratio weighting**: optional per‑sample `ratio` scales its contribution to base composition and balance metrics.
- **Constraints support** (JSON/YAML):
  - **Must-include**: specify that certain samples must be placed in particular lanes.
  - **Pre-grouping / lane quotas**: specify that certain lanes must contain a given **count** of samples (useful for pre-batched or pre-grouped inputs).
- **Flexible inputs**: accepts **.xlsx**, **.csv**, or **.tsv**. Required columns: `id,i7,i5`.  
  - In **Group mode**, `ratio` and `type` are **optional** (default `ratio=1.0`, default `type="ALL"`).  
  - In **Check mode**, `ratio` and `lane` are **optional** (default `ratio=1.0`, default `lane="ALL"`).
- **Cycle weighting**: emphasize cycles 3 & 4 via `--w34`, remaining weight spread over ≥C5 (normalized).

---

## Usage
### 1) Group (optimize assignment across N lanes)

**Without constraints:**
```bash
python index_balance_v2.py group \
library.tsv \
--lanes 3 \
--out lanes_assign.csv \
--rc-i5 \
[--cycles 8 --iters 20000 --seed 42]
```

**With constraints (JSON or YAML):**
```bash
python index_balance_v2.py group \
library.tsv \
--constraints constrain.json \
--out lanes_assign.csv \
--rc-i5 \
[--cycles 8 --iters 20000 --seed 42]
```

**Input columns (Group):**

- **Required columns**: `id,i7,i5`
- **Optional:** `ratio` (float; default 1.0), `type` (string; default "ALL")
- `id`: sample ID; `i7`: i7 sequence; `i5`: i5 sequence; `ratio`: number of aimed reads (can be scaled); `type`: labels of pre-groups (if exist)

**Constraints file examples:**

```json
{
  "lanes": {
    "1": { "quotas": { "A": 8, "B": 4, "C": 8 }, "must_include": [] },
    "2": { "quotas": { "A": 7, "B": 10, "C": 9 }, "must_include": ["sampleX"] }
  }
}
```
2 lanes    
assign 8 samples of "type A", 4 samples of "type B", and 8 samples of "type C" to lane 1    
assign 7 samples of "type A", 10 samples of "type B", and 9 samples of "type C" to lane 2. "sampleX" must be included.    

- `must_include`: place the listed sample IDs into the given lane indices before optimization.  
- `quotas`: request that each lane receive the specified number of samples (pre‑grouping). The optimizer fills remaining slots while balancing indices.

### 2) Check (evaluate existing lane assignments)
```bash
python index_balance_v2.py check lanes_assign.csv --rc-i5
```
**Input columns (Check):**

- **Required columns**: `id,i7,i5`
- **Optional:** `ratio` (default 1.0), `lane` (default "ALL")

The report prints per‑cycle **green/blue/dark** fractions and aggregated C1/C2 violations per lane.

### Common options
- `--i7-cycles / --i5-cycles`: override cycle counts; default = minimum observed length.  
- `--min-green / --min-blue / --max-dark`: C1/C2 feasibility thresholds (defaults 0.25, 0.05, 0.40).  
- `--c12-weight`: weight for cycles 1 & 2 (C1/C2) violations (default 1e6; feasibility‑first).  
- `--w34`: total mass assigned to C3/C4 (evenly split), remaining mass spread across ≥C5 (normalized to 1).  
- `--rc-i5`: whether to reverse‑complement i5 before evaluation.

---

## Calculation Method

### Base‑to‑color mapping (Illumina 2‑channel chemistry)
Per cycle, base fractions are converted to channel fractions:
- **Green** = pT + 0.5·pC  
- **Blue**  = pA + 0.5·pC  
- **Dark**  = pG 
so **green + blue + dark = 1**.

### C1/C2 feasibility (cycles 1 and 2 only)
- A violation is added when:
  - `green <= min_green`  
  - `blue <= min_blue` 
  - `dark >= max_dark`  
- The C1/C2 sum is multiplied by `c12_weight` to ensure feasibility dominates optimization.

### ≥C3 balance penalty
- For cycles ≥3, per‑cycle penalties encourage:
  - green ≥ 0.2 and blue ≥ 0.2  
  - dark ≤ 0.4  
  - green and blue not excessively imbalanced (|green − blue| small)
- Cycles are aggregated with the **cycle weights** (via `--w34` policy).
- **Ratio weighting**: when `ratio` is present, per‑sample contributions to base counts are scaled accordingly.

### Lane cost & optimization
- **Lane cost** = `c12_weight * (C1/C2 violations) + (≥C3 penalty)`
- Initialization: round‑robin seeding (after applying constraints such as `must_include` and `lane_quotas` when provided).
- Local search: attempt **swap/move** between lanes; keep changes that reduce total cost, thereby distributing all samples while improving balance.

---

## Notes
Additional remarks:
- If all `ratio=1.0`, results reduce to unweighted balance.  
- Larger `ratio` makes that sample heavier in both feasibility and ≥C3 balance calculations.  
- Constraints let you combine fixed placements (must‑include) and pre‑grouping quotas while still optimizing overall balance.




# Pool Assist Tool

## Overview
This script (`pool_assist.py`) automates pooling calculations for sequencing libraries.  
It calculates sample-specific pipetting volumes based on concentration, fragment size, ratios, and target lane volume.  
The algorithm ensures index balance across lanes while respecting pipetting constraints.

---

## Usage
```bash
python pool_assist.py <input_file.xlsx> <output_file.txt>
```

### Input file format
- accepts **.xlsx**, **.csv**, or **.tsv**
- Expected columns (in order, must include header line):
  1. **Sample ID**
  2. **Qubit Concentration (ng/µL)**
  3. **Average Size (bp)**
  4. **Ratio** (number of aimed reads, can be scaled per lane)
  5. **Lane ID**
  6. **Aimed Total Lane Volume (µL)**

### Output
- A plain text file with recommended pipetting volumes for each sample
- Useful information: 
  - Add water
  - pre-dilute: Yes / No
  - dilute_factor: 1 µl stock solution mix with (dilute_factor – 1) µl water
  - take(µL): take final volume from diluted working solution


---

## Fixed Parameters (in script)
- `MIN_PIPETTE_UL = 1.0` → minimum pipettable volume in µL  
- `SCALE_WHEN_OVERFLOW = False` → whether to rescale volumes if total exceeds aimed volume  
- `ROUND_NDIGITS = 3` → decimal precision for reported volumes
- 10 nM/µl → final concentration of the pooled library

---

## Notes
- Ratios are normalized per lane.  
- All samples within a lane share the same final target volume (`V_aim`).  
- Overflow scaling may lead to pipetting volumes below the minimum threshold.  
