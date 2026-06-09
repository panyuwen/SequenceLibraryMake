- [Index Balance Tool](#index-balance-tool)
  - [Overview](#overview)
  - [Key Features](#key-features)
  - [Usage](#usage)
    - [1) Group (optimize assignment across N lanes)](#1-group-optimize-assignment-across-n-lanes)
    - [2) Check (evaluate existing lane assignments)](#2-check-evaluate-existing-lane-assignments)
    - [Common options](#common-options)
  - [Calculation Method](#calculation-method)
    - [Base‑to‑color mapping (Illumina 2‑channel chemistry)](#basetocolor-mapping-illumina-2channel-chemistry)
    - [C1/C2 feasibility (cycles 1 and 2 only)](#c1c2-feasibility-cycles-1-and-2-only)
    - [≥C3 balance penalty](#c3-balance-penalty)
    - [Lane cost & optimization](#lane-cost--optimization)
  - [Notes](#notes)
- [Index Balance Tool (v3)](#index-balance-tool-v3)
  - [What's New in v3](#whats-new-in-v3)
  - [Usage (v3)](#usage-v3)
  - [Input Columns (v3)](#input-columns-v3)
  - [Feasibility Analysis & Removal Plans](#feasibility-analysis--removal-plans)
  - [Manual Constraints (v3)](#manual-constraints-v3)
  - [Freeze Mode — partial assignment + ballast (v3)](#freeze-mode--partial-assignment--ballast-v3)
  - [Objective Priority](#objective-priority)
  - [New & Changed Options (v3)](#new--changed-options-v3)
  - [Notes (v3)](#notes-v3)
- [Pool Assist Tool](#pool-assist-tool)
  - [Overview](#overview)
  - [Usage](#usage)
    - [Input file format](#input-file-format)
    - [Output](#output)
  - [Fixed Parameters (in script)](#fixed-parameters-in-script)
  - [Notes](#notes)
- [Understanding Index Orientation in NovaSeq X / X Plus](#understanding-index-orientation-in-novaseq-x--x-plus)
  - [Example 1 — *i7* (Index Read 1)](#example-1--i7-index-read-1)
  - [Example 2 — *i5* (Index Read 2)](#example-2--i5-index-read-2)
  - [Summary Table](#summary-table)
- [Index Hopping Counter](#index-hopping-counter)
  - [Overview](#overview)
  - [Usage](#usage)
    - [Parameters](#parameters)
    - [Output Files](#output-files)
  - [Notes](#notes)
    
    
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
[--cycles 8 --iters 20000 --seed 42] ## find the explanation of --rc-i5 below
```

**With constraints (JSON or YAML):**
```bash
python index_balance_v2.py group \
library.tsv \
--constraints constrain.json \
--out lanes_assign.csv \
--rc-i5 \
[--cycles 8 --iters 20000 --seed 42] ## find the explanation of --rc-i5 below
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
python index_balance_v2.py check \
lanes_assign.csv \
--rc-i5 ## find the explanation of --rc-i5 below
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




# Index Balance Tool (v3)

## Overview
`index_balance_v3.py` is a successor to `index_balance_v2.py`. It keeps the same 2‑channel color‑balance machinery (C1/C2 feasibility, ≥C3 weighted balance, ratio weighting, the `check` / `group` subcommands) and adds three things that matter when libraries reuse barcodes or are pre‑pooled:

1. **Automatic duplicate‑index separation.** Two samples sharing the same `i7+i5` barcode cannot be demultiplexed if they land in the same lane. v3 detects every duplicated barcode up front, builds a conflict graph, and produces a **conflict‑free assignment by construction** (graph coloring) before optimizing color balance — no manual specification required.
2. **Atomic bundles.** Pre‑pooled libraries (e.g. an exome capture pooled before sequencing) cannot be split. Add a **`bundle`** column; rows that share a bundle id are one atomic unit and always go to the same lane.
3. **Feasibility analysis.** For a given lane count, v3 says whether all samples fit conflict‑free. If not, it reports the **minimum lanes required** and proposes ranked **sample‑removal plans** to fit the lanes you have.

> **Correctness note:** v3 judges barcode conflicts on the **full** index strings (after optional `--rc-i5`), which is what demultiplexing actually uses. (`--cycles` still only controls the color‑balance evaluation window, exactly as in v2.)

## What's New in v3
| Capability | v2 | v3 |
|---|---|---|
| Duplicate `i7+i5` in one lane | penalized, may persist | detected up front, **separated by construction** |
| Pre‑pooled libraries | not modeled | **`bundle` column** keeps them together (atomic) |
| Too many samples for N lanes | not reported | **min lanes needed + removal plans** |
| Lane loading balance | fixed by round‑robin | **ratio‑sum balancing** (`--balance-weight`) |
| Conflict definition | truncated window | **full barcode** (demux‑correct) |
| `--constraints` (quotas/must‑include) | supported | **supported, with top priority** — honored first, then auto‑optimized within |
| Partial assignment (some lanes pre‑fixed) | not supported | **`--into` freeze mode** — honor existing lanes, fill only the blanks |
| No‑index pre‑pools as read "ballast" | not modeled | **ballast rows** count toward balance, skipped for conflict/color |

## Usage (v3)
### Group — assign all samples into N lanes
```bash
python index_balance_v3.py group \
  library.csv \
  --lanes 3 \
  --out lanes_assign.csv \
  [--rc-i5 --balance-weight 10 --max-plans 3 --auto-expand]
```

### Check — evaluate an existing layout
```bash
python index_balance_v3.py check lanes_assign.csv [--rc-i5]
```
On top of the v2 per‑cycle color report, `check` now lists any **duplicated barcodes** within each lane and the **minimum pairwise Hamming distance** of the combined `i7+i5` (flags collision‑risk pairs at distance ≤ 2).

## Input Columns (v3)
- **Required:** `id, i7, i5`
- **Optional:** `ratio` (float, default 1.0), `type` (default `"ALL"`), `lane` (default empty; read by `check` and by `group --into`), **`bundle`** (default empty)
- **`bundle`** — rows sharing the same non‑empty value form one atomic unit assigned to a single lane. Empty = singleton.

Example with a bundle (the two `POOL_EX1` rows always land in the same lane):
```csv
id,i7,i5,ratio,bundle
WGS_01,ACGTACGT,TGCATGCA,1,
EX_a,GACTGACT,CTGACTGA,1,POOL_EX1
EX_b,TCAGTCAG,AGTCAGTC,1,POOL_EX1
```

## Feasibility Analysis & Removal Plans
Running `group` first prints a feasibility verdict.

**Feasible** → an optimized conflict‑free, color‑balanced, ratio‑balanced assignment is produced (best found by local search, not guaranteed globally optimal; written to `--out`):
```
✅ FEASIBLE: all 11 samples fit into 3 lane(s) with no barcode conflicts.
   1 shared-barcode group(s) auto-separated across lanes.
```

**Infeasible** → (1) the minimum lanes required, and (2) ranked removal plans, each showing per‑lane composition and ratio balance:
```
❌ INFEASIBLE: 6 samples cannot fit into 2 lane(s) without barcode conflicts.
   (1) Minimum lanes required: 3
   (2) Removal options to fit into 2 lane(s):
--- Removal Plan 1: drop 1 unit(s) / 1 sample(s), ratio=1.00 ---
   remove sample T1  (n=1, ratio=1.00)  ids: T1
     Lane 1: 3 samples, ratio=4.00, colorbal=FAIL  ids: T2,U1,S9
     Lane 2: 2 samples, ratio=4.00, colorbal=FAIL  ids: T3,U2
```
- Plans are ranked by least material (`ratio`) dropped, and genuinely distinct alternatives are surfaced (drop `T1` vs `T2` vs `T3`).
- This covers both **over‑subscribed barcodes** (one barcode used in more lanes than exist) **and structural cycles among bundles** — e.g. three bundles that pairwise share three different barcodes need 3 lanes even though every barcode is only used twice.
- With `--out`, each plan's assignment is also written next to your output path with a `_planN` suffix.

**`--auto-expand`** → when the requested lanes are infeasible, automatically solve at the minimum feasible number of lanes instead of printing removal plans:
```
❌ INFEASIBLE: ... cannot fit into 2 lane(s) ...
--auto-expand set: solving at 3 lanes.
✨ All lanes free of barcode conflicts.
(Note: requested 2 lanes was infeasible; used 3 lanes instead.)
```

## Manual Constraints (v3)
For explicit, manual control, pass a JSON/YAML constraints file with `--constraints`. **It takes priority over the automatic engine** and **sets the lane count** (so `--lanes` is not needed):

```bash
python index_balance_v3.py group library.csv --constraints layout.json --out lanes_assign.csv
```

```json
{
  "lanes": {
    "1": { "quotas": { "A": 2, "B": 1 } },
    "2": { "quotas": { "A": 2, "B": 1 }, "must_include": ["A4"] }
  }
}
```
- **`quotas`** — each lane must receive exactly this many samples of each `type` (pre‑grouping). The sum of all quotas must equal the number of input samples, and there must be enough samples of each type.
- **`must_include`** — pin specific samples to a lane *before* optimization. An id may be a **sample id** or a **bundle id** (pinning a bundle pins all its members). Bundles always stay together and consume their members' types from the lane's quota.

**Priority order:** constraints are honored first; the automatic engine then optimizes *within* that fixed structure — it separates duplicate barcodes, satisfies cycle‑1/2 color, and balances ratios using only **quota‑preserving swaps** (it never moves `must_include` units or changes a lane's per‑type counts). Example:
```
Applying --constraints (priority): 2 lanes, 6 quota slots, 1 must_include pin(s).
✅ Constraints satisfied (quotas + must_include).
✅ Lane 1 unique check passed.
✅ Lane 2 unique check passed.
```
If the quotas force a barcode collision that no quota‑preserving swap can resolve, v3 reports the colliding samples and **refuses to write** the file:
```
❌ [CRITICAL] Lane 1 duplicate barcode(s): AAAACCCC+GGGGTTTT -> A1, A2
!! Conflict detected — DO NOT sequence this layout !!
⚠️  Refusing to write out.csv: the assignment contains barcode conflicts. No file written.
```
(Relax the quota or move one of the named samples to a different lane, then re‑run.)

> **YAML** is supported if `PyYAML` is installed; otherwise use JSON.

## Freeze Mode — partial assignment + ballast (v3)
For the common real-world case where **some lanes are already assigned** and you only want to **fill the rest**, pass `--into`:

```bash
python index_balance_v3.py group sample_mix.xlsx --into 1,2,3 --rc-i5 --out lanes.csv
```

How `--into 1,2,3` interprets the input `lane` column:
- **Row already has a lane in `{1,2,3}`** → **frozen** there (not moved).
- **Row pinned to a lane *not* in `--into`** (e.g. `4`, `5`) → **ignored** (already finalized elsewhere).
- **Blank-lane row** → **distributed** into lanes 1/2/3.

Plus two behaviors that make a mixed real sheet "just work":
- **No-index ballast.** A unit with empty `i7`/`i5` (e.g. a pre-pooled complex library) is treated as **ratio-only ballast**: it counts toward its lane's read budget (so balancing is correct) but is **excluded from barcode-conflict and color** scoring. It must carry a lane. *Ballast is decided per unit:* a `bundle` must be either all-indexed or all-no-index — a bundle mixing indexed and no-index members is rejected (so indexed members are never silently dropped). A plain (non-bundle) row is its own unit, so "row with empty index = ballast" still holds.
- **Automatic color scope.** Color (C1/C2 + ≥3) is optimized **only for lanes whose samples are all fully indexed**. A lane containing any pool/ballast is left to the pool's own complexity (reported as `color: N/A (pool-covered)`).

> `--into` and `--constraints` are different assignment mechanisms and **cannot be combined** (the tool errors if both are given).

**Priority in freeze mode** (differs from default group mode — matches pooling intent):

| Priority | Objective |
|---|---|
| 1 (hard) | No barcode conflict among indexed samples in a lane |
| 2 | **Balanced per-lane read totals (ballast included)** — *ensure* |
| 3 | Color of the fully-indexed lanes — *best-effort* |

i.e. read-count balance is **ensured first**; color is optimized within balanced layouts (so the optimizer won't sacrifice balance to "hide" awkward samples in pool-covered lanes). Conflicts are resolved automatically (duplicates separated across the `--into` lanes); if a barcode is over-subscribed for the number of `--into` lanes, it reports the offending sample. The `--out` CSV is the **whole table** with the `lane` column filled (frozen + newly assigned + ignored rows all preserved).

## Objective Priority
The optimizer minimizes a layered objective (high → low):

| Priority | Objective | Knob |
|---|---|---|
| 0 (hard, optional) | User constraints: quotas + must_include | `--constraints` |
| 1 (hard) | No barcode conflict within a lane | — |
| 2 | Cycle 1/2 color feasibility (i7 & i5) | `--c12-weight` (default 1e6) |
| 3 | Balanced per‑lane ratio sums (loading) | `--balance-weight` (default 10) |
| 4 | Cycle ≥3 color balance | `--w34` |

With realistic (color‑balanced) index kits, priority 2 is met in any layout, so lanes come out **ratio‑balanced**; with pathological indexes, color feasibility correctly wins over balance. Set `--balance-weight 0` for v2‑style behavior (fixed round‑robin lane sizes, color‑only optimization).

## New & Changed Options (v3)
| Option | Description | Default |
|---|---|---|
| `--constraints` | JSON/YAML per‑lane quotas + must_include; **takes priority** and sets the lane count | none |
| `--into` | Freeze mode: honor the existing `lane` column (fixed rows + no‑index ballast); distribute blank‑lane rows into these lanes, ignore rows pinned elsewhere | none |
| `--balance-weight` | Weight for balancing per‑lane ratio sums (`0` = v2‑style, color‑only) | 10.0 |
| `--max-plans` | Number of removal plans to show when infeasible | 3 |
| `--auto-expand` | If infeasible, solve at the minimum feasible #lanes instead | off |
| `bundle` (input column) | Atomic pooling unit; rows with the same `bundle` value stay in one lane | empty |

`--lanes` is required **unless** `--constraints` or `--into` is given (each sets the lanes itself). All shared v2 options behave the same: `--cycles`, `--i7-cycles`, `--i5-cycles`, `--min-green`, `--min-blue`, `--max-dark`, `--w34`, `--rc-i5`, `--iters`, `--seed`, `--c12-weight`, `--out`. (In `--into` mode, balance is prioritized over color, so `--balance-weight`/`--c12-weight` are set internally.)

## Notes (v3)
- **Unfixable bundles:** if a pre‑pooled bundle contains two members with the *same* barcode, no lane assignment can demultiplex it — v3 reports `UNFIXABLE`, names the colliding members, and stops (fix the pool, then re‑run).
- **Undetermined verdict:** on very large or tangled bundle graphs the exact colorability search may exceed its step budget; v3 then prints `⚠️ UNDETERMINED` (may still be feasible) rather than a false `INFEASIBLE`.
- **Write safety:** v3 runs a final conflict check and **refuses to write** `--out` if any collision remains.
- For singleton libraries with no shared barcodes, v3 reduces to v2‑style balanced lanes plus color optimization.

---


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
- `10 nM/µl` → final concentration of the pooled library
- `v_raw < MIN_PIPETTE_UL & concentration ≥ 4 nM` → to dilute

---

## Notes
- Ratios are normalized per lane.  
- All samples within a lane share the same final target volume (`V_aim`).  
- Overflow scaling may lead to pipetting volumes below the minimum threshold.  



# Understanding Index Orientation in NovaSeq X / X Plus

Illumina sequencing workflows differ in how they read the two index strands.  
NovaSeq X and X Plus operate exclusively in the **Reverse-Complement (RC) workflow**, which affects the *i5* index orientation but not *i7*.
    
## Example 1 — *i7* (Index Read 1)

**Full adapter sequence**  
```
5′-CAAGCAGAAGACGGCATACGAGAT CTACAGTG GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC*T-3′
```

- The embedded **i7 index** in this adapter is present in the *reverse-complement* orientation (`CTACAGTG`).  
- When sequenced in **Index Read 1**, the NovaSeq X instrument reads it in the **forward direction**, producing the **index read** `CACTGTAG`.  
- Thus, the instrument “sees” the *forward* index sequence.  
- For **index-balance analysis**, use the original index (`CACTGTAG`).  
- For **demultiplexing**, both i7 and i5 need to be reverse-complemented to match the raw BCL reads.

| Context | Direction seen by instrument | Correct form to use |
|:--|:--|:--|
| Adapter embedding | Reverse complement of index | — |
| Sequencing read | Forward (CACTGTAG) | ✅ Use as is for balance |
| Demultiplexing match | — | 🔁 Reverse complement both indices |

---

## Example 2 — *i5* (Index Read 2)

**Full adapter sequence**  
```
5′-AATGATACGGCGACCACCGAGATCTACAC AAGCGACT ACACTCTTTCCCTACACGACGCTCTTCCGATC*T-3′
```

- The **i5 index** is embedded in the *forward* orientation (`AAGCGACT`) within the P5 adapter.  
- Under the **RC workflow**, NovaSeq X reads this region in the *reverse-complement* direction, yielding `AGTCGCTT` during Index Read 2.  
- Therefore, the instrument “sees” the **reverse complement** of the labeled i5 sequence.  
- For **index-balance analysis**, use the reverse-complemented i5 (`AGTCGCTT`).  
- For **demultiplexing**, again reverse-complement both i7 and i5 to match the instrument output.

| Context | Direction seen by instrument | Correct form to use |
|:--|:--|:--|
| Adapter embedding | Forward (AAGCGACT) | — |
| Sequencing read | Reverse complement (AGTCGCTT) | ✅ Use RC for balance |
| Demultiplexing match | — | 🔁 Reverse complement both indices |

---

## Summary Table

| Index | Adapter orientation | Sequencing direction (NovaSeq X) | For balance analysis | For demultiplexing |
|:--|:--|:--|:--|:--|
| i7 | Reverse complement | Forward | Original | Reverse complement |
| i5 | Forward | Reverse complement | Reverse complement | Reverse complement |

---

**Key takeaway:**  
> In NovaSeq X / X Plus (RC workflow), i7 is read in the forward direction, while i5 is read as its reverse complement.  
> Consequently, both indices should be reverse-complemented for demultiplexing, but only i5 should be reverse-complemented for index-balance evaluation.


# Index Hopping Counter


## Overview
This tool calculates **index hopping events** within a lane by comparing observed barcode pairs from Undetermined FASTQ files to a list of valid index pairs.

It supports:
- Dual or single index libraries  
- Automatic or manual trimming of barcode lengths  
- Hamming distance–based mismatch tolerance (`--mismatch`)  
- Classification of reads into:
  - **valid** — matched known index pairs  
  - **invalid** — both indices match known ones, but pair does not exist (index hopping)  
  - **unknown** — cannot uniquely match within mismatch limit or contains `N`


---

## Usage

### Basic command
```bash
python idx_hopping_count.py \
  --r1 Undetermined_S0_L001_R1_001.fastq.gz \
  --valid_pairs SampleSheet.csv \
  --out_prefix lane1_hop
```

### With custom trimming or mismatch
```bash
python hop_from_reads_trimlen.py \
  --r1 Undetermined_S0_L002_R1_001.fastq.gz \
  --valid_pairs my_pairs.csv \
  --len_i7 10 --len_i5 10 \
  --mismatch 0
```

---

### Parameters

| Argument | Description | Default |
|-----------|--------------|----------|
| `--r1` | Undetermined R1 FASTQ file (gzipped) | *Required* |
| `--valid_pairs` | CSV/TSV or Illumina SampleSheet containing valid i7/i5 combinations | *Required* |
| `--len_i7`, `--len_i5` | Expected index lengths (trim observed barcodes before comparison) | 8 |
| `--mismatch` | Allowed mismatches per index (Hamming distance) | 1 |
| `--orientation` | Orientation handling: `auto`, `as_is`, `rc_i7`, `rc_i5`, `rc_both` | `rc_both` |
| `--out_prefix` | Prefix for all output files | `hop_report` |
| `--topn` | Top-N invalid pairs shown in console | 10 |

---

### Output Files

| File | Description |
|------|--------------|
| `<prefix>.summary.tsv` | Summary statistics (valid/invalid/unknown counts, hopping rate) |
| `<prefix>.invalid_pairs.csv` | Index-hopped pairs and their counts |


---

## Notes
- Reads containing `N` in either index are automatically classified as *unknown*.  
- Only dual indices where both sides uniquely match known index sets but form a non-existent pair are counted as *invalid*.  
- Default mismatch logic (`--mismatch 1`) applies **independently** to i7 and i5.  

