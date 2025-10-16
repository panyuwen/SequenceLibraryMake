#!/usr/bin/env python3
import argparse
import csv
import gzip
import sys
import os
from collections import Counter

def hamming_distance(a: str, b: str):
    """Return Hamming distance if same length else None."""
    if len(a) != len(b):
        return None
    return sum(x != y for x, y in zip(a, b))

def match_unique(obs: str, candidates: set, mm: int):
    """
    Return a uniquely matched candidate within <= mm mismatches.
    If none or multiple best matches at same minimal distance, return None.
    """
    hits = []
    if ('N' not in obs) and ('n' not in obs):
        for cand in candidates:
            d = hamming_distance(obs, cand)
            if d is not None and d <= mm:
                hits.append((d, cand))
    if not hits:
        return None
    hits.sort(key=lambda x: (x[0], x[1]))  # sort by distance then lexicographically
    min_d = hits[0][0]
    winners = [c for d, c in hits if d == min_d]
    return winners[0] if len(winners) == 1 else None

def is_valid_pair(i7_obs: str, i5_obs: str, valid_pairs: set, i7_bank: set, i5_bank: set, mm: int):
    """
    Try to uniquely map observed i7/i5 to known banks within <= mm mismatches (each).
    Returns (is_valid_bool_or_None, matched_i7, matched_i5):
      - True  : both uniquely matched and (i7,i5) in valid_pairs
      - False : both uniquely matched but (i7,i5) NOT in valid_pairs
      - None  : at least one side failed unique matching (unknown)
    """
    i7_match = match_unique(i7_obs, i7_bank, mm) if i7_bank else None
    i5_match = match_unique(i5_obs, i5_bank, mm) if i5_bank else None
    if (i7_bank and not i7_match) or (i5_bank and not i5_match):
        return None, i7_match, i5_match
    # If no i5_bank (single-index), treat as matched with empty string
    if not i5_bank:
        i5_match = ''
    return ((i7_match, i5_match) in valid_pairs), i7_match, i5_match


def reverse_complement(seq: str) -> str:
    comp = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
    return seq.translate(comp)[::-1]

def looks_like_samplesheet(path: str) -> bool:
    try:
        with open(path, 'r', newline='') as f:
            for line in f:
                if line.strip().lower() == '[data]':
                    return True
    except Exception:
        return False
    return False

def load_valid_pairs(path: str):
    """
    Load valid (i7,i5) pairs as uppercase strings.
    Returns: set(pairs), set(len_i7), set(len_i5)
    """
    pairs = set()
    len_i7_set = set()
    len_i5_set = set()
    if looks_like_samplesheet(path):
        in_data = False
        headers = None
        with open(path, 'r', newline='') as f:
            r = csv.reader(f)
            for row in r:
                if not row:
                    continue
                if row[0].strip().lower() == '[data]':
                    in_data = True
                    headers = None
                    continue
                if not in_data:
                    continue
                if headers is None:
                    headers = [h.strip() for h in row]
                    def find_col(*names):
                        for n in names:
                            for h in headers:
                                if h.strip().lower() == n.lower():
                                    return headers.index(h)
                        return None
                    idx_i7 = find_col('index','i7_index','i7_index_id','index1')
                    idx_i5 = find_col('index2','i5_index','i5_index_id','index2')
                    if idx_i7 is None:
                        raise ValueError("SampleSheet: cannot find i7 column (tried: index, i7_index, i7_index_id, index1)")
                    continue
                if headers is not None:
                    if len(row) < len(headers):
                        row = row + ['']*(len(headers)-len(row))
                    def get_col(name_opts, default=''):
                        for nm in name_opts:
                            for i,h in enumerate(headers):
                                if h.strip().lower() == nm.lower():
                                    return row[i].strip()
                        return default
                    i7 = get_col(['index','i7_index','i7_index_id','index1']).upper()
                    i5 = get_col(['index2','i5_index','i5_index_id','index2']).upper()
                    if i7:
                        pairs.add((i7, i5))
                        len_i7_set.add(len(i7))
                        len_i5_set.add(len(i5))
    else:
        with open(path, 'r', newline='') as f:
            sample = f.read(4096)
            f.seek(0)
            try:
                dialect = csv.Sniffer().sniff(sample, delimiters=',\t; ')
            except Exception:
                class SimpleDialect(csv.Dialect):
                    delimiter = ','
                    quotechar = '"'
                    escapechar = None
                    doublequote = True
                    lineterminator = '\n'
                    quoting = csv.QUOTE_MINIMAL
                dialect = SimpleDialect()
            r = csv.reader(f, dialect)
            headers = next(r)
            headers_l = [h.strip().lower() for h in headers]
            def idx_of(*cands):
                for c in cands:
                    if c.lower() in headers_l:
                        return headers_l.index(c.lower())
                return None
            i7_idx = idx_of('i7','index','index1')
            i5_idx = idx_of('i5','index2')
            if i7_idx is None:
                raise ValueError("Valid-pairs file: cannot find i7 header (try headers: i7/index/index1)")
            for row in r:
                if not row:
                    continue
                if len(row) < len(headers):
                    row = row + ['']*(len(headers)-len(row))
                i7 = row[i7_idx].strip().upper()
                i5 = row[i5_idx].strip().upper() if i5_idx is not None else ''
                if i7:
                    pairs.add((i7, i5))
                    len_i7_set.add(len(i7))
                    len_i5_set.add(len(i5))
    if not pairs:
        raise ValueError("No valid i7/i5 pairs loaded from: "+path)
    len_i5_set = {l for l in len_i5_set if l > 0}
    return pairs, len_i7_set, len_i5_set

def parse_barcode_from_header(h: str):
    """
    CASAVA 1.8+ header format:
    @<instr>:<run>:<flowcell>:<lane>:<tile>:<x>:<y> <read>:<is_filtered>:<control>:<barcode>
    Returns barcode string like 'ATCGTA+GTCAGT' or 'ATCGTA'.
    """
    h = h.strip()
    if not h or h[0] != '@':
        return None
    parts = h.split(' ')
    if len(parts) < 2:
        return None
    right = parts[1]
    fields = right.split(':')
    if len(fields) < 4:
        return None
    barcode = fields[-1].strip()
    if barcode == '' or barcode == 'N' or all(c == 'N' for c in barcode.replace('+','')):
        return None
    return barcode

def trim_to_lengths(i7: str, i5: str, L7: int = None, L5: int = None):
    if L7 is not None and len(i7) >= L7:
        i7 = i7[:L7]
    if L5 is not None and len(i5) >= L5:
        i5 = i5[:L5]
    return i7, i5

def stream_pairs_from_fastq(r1_path: str, L7: int = None, L5: int = None):
    with gzip.open(r1_path, 'rt') as f1:
        while True:
            h1 = f1.readline()
            if not h1:
                break
            _ = f1.readline(); _ = f1.readline(); _ = f1.readline()
            bc = parse_barcode_from_header(h1)
            if bc is None:
                continue
            if '+' in bc:
                i7, i5 = bc.split('+', 1)
            else:
                i7, i5 = bc, ''
            i7, i5 = trim_to_lengths(i7.upper(), i5.upper(), L7, L5)
            yield (i7, i5)

def choose_orientation(valid_pairs, sample_iter, limit=200000):
    modes = {
        'as_is': lambda i7,i5: (i7, i5),
        'rc_i7': lambda i7,i5: (reverse_complement(i7), i5),
        'rc_i5': lambda i7,i5: (i7, reverse_complement(i5)),
        'rc_both': lambda i7,i5: (reverse_complement(i7), reverse_complement(i5)),
    }
    valid_set = set(valid_pairs)
    samples = []
    for i, pair in enumerate(sample_iter):
        samples.append(pair)
        if i+1 >= limit:
            break
    scores = {}
    for name, fn in modes.items():
        m = 0
        for i7,i5 in samples:
            a,b = fn(i7,i5)
            if (a,b) in valid_set:
                m += 1
        scores[name] = m
    best = max(scores, key=scores.get) if scores else 'as_is'
    return best, modes[best], scores, len(samples)

def main():
    ap = argparse.ArgumentParser(description="Estimate index hopping from Undetermined R1 FASTQ by parsing i7+i5 from headers with length trimming.")
    ap.add_argument("--r1", required=True, help="Undetermined R1 FASTQ.gz")
    ap.add_argument("--valid_pairs", required=True, help="CSV/TSV with i7,i5 OR Illumina SampleSheet.csv")
    ap.add_argument("--out_prefix", default="hop_report", help="Prefix for outputs")
    ap.add_argument("--topn", type=int, default=10, help="Top-N invalid pairs to print")
    ap.add_argument("--orientation", choices=["auto","as_is","rc_i7","rc_i5","rc_both"], default="rc_both",
                    help="Orientation of observed vs valid indices; if you are not sure, try 'auto' on a well demultiplexed sample, which will apply all and picks the best.")
    ap.add_argument("--len_i7", type=int, default=8, help="Expected i7 length to trim observed barcodes (e.g., 8)")
    ap.add_argument("--len_i5", type=int, default=8, help="Expected i5 length to trim observed barcodes (e.g., 8)")
    ap.add_argument("--mismatch", type=int, default=1,
                    help="Allowed mismatches per index (i7 and i5 each). Unique match required within this threshold. Default: 0 (exact).")
    args = ap.parse_args()

    valid_pairs, len_i7_set, len_i5_set = load_valid_pairs(args.valid_pairs)

    L7 = args.len_i7 if args.len_i7 is not None else (next(iter(len_i7_set)) if len(len_i7_set) == 1 else None)
    L5 = args.len_i5 if args.len_i5 is not None else (next(iter(len_i5_set)) if len(len_i5_set) == 1 else None)

    sys.stderr.write(f"[INFO] Valid i7 lengths: {sorted(len_i7_set)}; using L7={L7}\n")
    sys.stderr.write(f"[INFO] Valid i5 lengths: {sorted(len_i5_set)}; using L5={L5}\n")

    orient_name = args.orientation
    transform = lambda i7,i5: (i7,i5)
    if args.orientation == "auto":
        orient_name, transform, orient_scores, n_sampled = choose_orientation(valid_pairs, stream_pairs_from_fastq(args.r1, L7, L5))
    else:
        def make_tf(name):
            if name == "as_is": return lambda i7,i5: (i7,i5)
            if name == "rc_i7": return lambda i7,i5: (reverse_complement(i7), i5)
            if name == "rc_i5": return lambda i7,i5: (i7, reverse_complement(i5))
            if name == "rc_both": return lambda i7,i5: (reverse_complement(i7), reverse_complement(i5))
            raise ValueError("Unknown orientation")
        transform = make_tf(args.orientation)
        orient_scores = {}
        n_sampled = 0

    counts = Counter()
    total_pairs = 0
    for i7,i5 in stream_pairs_from_fastq(args.r1, L7, L5):
        a,b = transform(i7,i5)
        counts[(a,b)] += 1
        total_pairs += 1

    if total_pairs == 0:
        print("No reads found (or no barcodes in headers).", file=sys.stderr)
        sys.exit(2)

    valid_set = set(valid_pairs)
    i7_bank = {i7 for (i7, _i5) in valid_set}
    i5_bank = {i5 for (_i7, i5) in valid_set if i5 != ''}

    valid_count = 0
    invalid_count = 0
    unknown_count = 0
    invalid_set = Counter()

    for (i7_obs, i5_obs), c in counts.items():
        verdict, i7m, i5m = is_valid_pair(i7_obs, i5_obs, valid_set, i7_bank, i5_bank, args.mismatch)
        
        if verdict is True:
            valid_count += c
        elif verdict is False:
            # i7、i5 各自都唯一匹配到了合法库，但组合不在合法对里 → invalid（疑似 hopping）
            invalid_count += c
            invalid_set[(i7_obs, i5_obs)] += c
        else:
            # 其它情况（至少一端无法在容错内唯一匹配）→ unknown
            unknown_count += c

    os.makedirs(os.path.dirname(args.out_prefix), exist_ok=True) if os.path.dirname(args.out_prefix) else None

    with open(args.out_prefix + ".summary.tsv", "w") as f:
        f.write("metric\tvalue\n")
        f.write(f"total_pairs\t{total_pairs}\n")
        f.write(f"valid_pairs\t{valid_count}\n")
        f.write(f"idx_hop_pairs\t{invalid_count}\n")
        f.write(f"unknown_pairs\t{unknown_count}\n")
        f.write(f"orientation\t{orient_name}\n")
        f.write(f"trim_len_i7\t{L7 if L7 is not None else ''}\n")
        f.write(f"trim_len_i5\t{L5 if L5 is not None else ''}\n")
        if n_sampled:
            f.write(f"auto_orientation_samples\t{n_sampled}\n")
        for k,v in orient_scores.items():
            f.write(f"auto_orientation_score_{k}\t{v}\n")

    # with open(args.out_prefix + ".valid_pairs.csv", "w", newline='') as f:
    #     w = csv.writer(f)
    #     w.writerow(["i7","i5","count"])
    #     for (i7,i5), c in counts.most_common():
    #         if (i7,i5) in valid_set:
    #             w.writerow([i7,i5,c])

    with open(args.out_prefix + ".invalid_pairs.csv", "w", newline='') as f:
        w = csv.writer(f)
        w.writerow(["i7","i5","count"])
        for (i7,i5), c in invalid_set.most_common():
            w.writerow([i7,i5,c])

    print("\n=== Index Hopping from R1 FASTQ Headers ===")
    print(f"Total reads parsed: {total_pairs:,}")
    print(f"Orientation: {orient_name}")
    if n_sampled:
        print("Auto-orientation scores:", orient_scores)
    print(f"Primer lengths: i7={L7}, i5={L5}")
    print(f"Valid pairs:   {valid_count:,}")
    print(f"Index hopping pairs: {invalid_count:,}")
    print("\nTop hopping pairs:")
    shown = 0
    for (i7,i5), c in invalid_set.most_common():
        print(f"  ({i7},{i5})\t{c:,}")
        shown += 1
        if shown >= 10:
            break
    print(f"\nWrote: {args.out_prefix}.summary.tsv")
    # print(f"Wrote: {args.out_prefix}.valid_pairs.csv")
    print(f"Wrote: {args.out_prefix}.invalid_pairs.csv")

if __name__ == "__main__":
    main()
