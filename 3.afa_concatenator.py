#!/usr/bin/env python3
"""
Concatenator that uses strain_code_map.txt as canonical strain list.
More tolerant matching: if the numeric code extracted from an alignment header
doesn't exactly equal a code in the map, we attempt substring matches
(both directions) to associate the sequence to a mapped code.

Expected strain_code_map.txt format:
#strain_code\tstrain_name(s)
14693065046\tEcoli_BEN2908;alias

Behavior:
 - Only codes listed in strain_code_map.txt are included in the final output.
 - Sequences are taken from aligned OG files; missing sequences are filled by gap blocks of OG length.
"""
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# --------- Config ----------
input_dir = "output_OGs/Muscle_alignment"
output_dir = "output_concatenated"
strain_map_path = "output_OGs/strain_code_map.txt"
accepted_exts = (".afa", ".aln", ".fasta", ".fa")
MIN_DIGITS = 6
# --------------------------

def read_strain_map(path):
    """Read strain_code_map.txt; return dict code -> safe_filename_name."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"strain map not found: {path}")
    mapping = {}
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            code = parts[0].strip()
            names = parts[1].strip()
            first_name = names.split(";")[0].strip()
            # sanitize filename characters
            safe_name = "".join(ch if (ch.isalnum() or ch in "-._") else "_" for ch in first_name)
            mapping[code] = safe_name
    if not mapping:
        raise RuntimeError(f"No valid entries found in strain map: {path}")
    return mapping

def longest_digit_run(s):
    """Return the longest consecutive run of digits in s, or None."""
    best = ""
    current_digits = []
    for ch in s:
        if '0' <= ch <= '9':
            current_digits.append(ch)
        else:
            if current_digits:
                run = "".join(current_digits)
                if len(run) > len(best):
                    best = run
                current_digits = []
    if current_digits:
        run = "".join(current_digits)
        if len(run) > len(best):
            best = run
    return best if best else None

def extract_code(rec):
    """Extract a numeric code from rec.id or rec.description via simple heuristics."""
    # 1) try second part after split('.')
    parts = rec.id.split(".")
    if len(parts) > 1:
        candidate = parts[1]
        if candidate.isdigit() and len(candidate) >= MIN_DIGITS:
            return candidate
    # 2) longest numeric run in id
    run = longest_digit_run(rec.id)
    if run and len(run) >= MIN_DIGITS:
        return run
    # 3) longest numeric run in description
    desc = getattr(rec, "description", "") or ""
    run = longest_digit_run(desc)
    if run and len(run) >= MIN_DIGITS:
        return run
    return None

def find_best_map_code(code_rec, map_codes):
    """
    Given code extracted from a record (code_rec, string) and an iterable of map_codes (strings),
    try to find the best map code to associate:
      1) exact match
      2) map_code in code_rec (substring)
      3) code_rec in map_code (substring)
    Returns the matched map_code or None.
    """
    if code_rec is None:
        return None
    # exact
    if code_rec in map_codes:
        return code_rec
    # map code is substring of code_rec
    for mc in map_codes:
        if mc in code_rec:
            return mc
    # code_rec is substring of map code
    for mc in map_codes:
        if code_rec in mc:
            return mc
    return None

def main():
    os.makedirs(output_dir, exist_ok=True)

    # 1) load strain map: these are the only codes we will produce
    code_to_name = read_strain_map(strain_map_path)
    target_codes = sorted(code_to_name.keys())
    print(f"Loaded {len(target_codes)} strain codes from '{strain_map_path}'")

    # 2) find alignment files
    files = sorted(f for f in os.listdir(input_dir) if f.lower().endswith(accepted_exts))
    if not files:
        raise FileNotFoundError(f"No alignment files found in '{input_dir}'")
    print(f"Found {len(files)} alignment files in '{input_dir}'")

    # 3) for each OG: build dict map_code -> sequence (only accept sequences that can be associated to map codes)
    ogs = []  # (fname, aln_len, dict map_code->sequence)
    map_code_set = set(target_codes)
    for fname in files:
        fpath = os.path.join(input_dir, fname)
        seqs = list(SeqIO.parse(fpath, "fasta"))
        if not seqs:
            print(f"Skipping empty file: {fname}")
            continue
        aln_len = max(len(rec.seq) for rec in seqs)
        code2seq = {}
        matched_count = 0
        for rec in seqs:
            code_rec = extract_code(rec)
            mapped = find_best_map_code(code_rec, map_code_set)
            if not mapped:
                continue
            # accept the sequence as-is (assume alignment already done)
            code2seq[mapped] = str(rec.seq)
            matched_count += 1
        ogs.append((fname, aln_len, code2seq))
        # small debug line to help you see mapping success per OG
        print(f"{fname}: {matched_count} sequences matched to {len(target_codes)} target codes")

    if not ogs:
        raise RuntimeError("No valid OGs processed.")

    # 4) prepare concatenation containers in the order of the strain map
    concat = {code: [] for code in target_codes}
    missing = {code: [] for code in target_codes}

    # 5) for each OG, append seq or gap block
    for fname, aln_len, code2seq in ogs:
        for code in target_codes:
            if code in code2seq:
                concat[code].append(code2seq[code])
            else:
                concat[code].append("-" * aln_len)
                missing[code].append(fname)

    # 6) write FASTA per strain using the mapped name for filenames
    for code in target_codes:
        seq_full = "".join(concat[code])
        fname_safe = f"{code_to_name[code]}.fasta"
        out_path = os.path.join(output_dir, fname_safe)
        rec = SeqRecord(Seq(seq_full), id=code, description="")
        with open(out_path, "w") as fh:
            SeqIO.write(rec, fh, "fasta")

    print(f"Finished. Wrote {len(target_codes)} concatenated FASTAs to '{output_dir}'.")

if __name__ == "__main__":
    main()
