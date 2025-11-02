#!/usr/bin/env python3
"""
Merge all single-strain FASTA files (one record per file) into a single multifasta
and write a PHYLIP file.

Input directory: output_concatenated/  (files like Ecoli_BEN2908.fasta)
Outputs:
  - all_strains_concatenated.fasta
  - all_strains_concatenated.phy   (PHYLIP relaxed)
"""

import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# ---- CONFIG ----
input_dir = "output_concatenated"
combined_fasta = "all_strains_concatenated.fasta"
combined_phylip = "all_strains_concatenated.phy"
accepted_exts = (".fasta", ".fa", ".fas")
use_phylip_classic = False  # set True to write classic PHYLIP (10-char names), see notes
# ----------------

def find_fasta_files(d):
    files = sorted(
        f for f in os.listdir(d)
        if os.path.isfile(os.path.join(d, f)) and f.lower().endswith(accepted_exts)
    )
    return [os.path.join(d, f) for f in files]

def shorten_name(name, length=10):
    """Make a filename/identifier safe truncated to `length` chars for classic PHYLIP."""
    # remove spaces and problematic chars, then pad/truncate
    safe = "".join(ch if (ch.isalnum() or ch in "-._") else "_" for ch in name)
    if len(safe) >= length:
        return safe[:length]
    else:
        return safe.ljust(length)

def main():
    os.makedirs(".", exist_ok=True)

    fasta_paths = find_fasta_files(input_dir)
    if not fasta_paths:
        raise FileNotFoundError(f"No fasta files found in '{input_dir}'. Expected one FASTA per strain.")

    records = []
    for p in fasta_paths:
        # Expect each file to contain exactly one record (the concatenated sequence)
        try:
            rec = SeqIO.read(p, "fasta")
        except ValueError:
            # If the file contains multiple records, take them all (but warn)
            recs = list(SeqIO.parse(p, "fasta"))
            if not recs:
                print(f"Warning: file {p} is empty; skipping.")
                continue
            if len(recs) > 1:
                print(f"Warning: file {p} contains {len(recs)} records; adding all.")
            records.extend(recs)
            continue
        # If id is numeric or not desired, you can use filename as id:
        # use filename without extension as record id to avoid collisions
        fname = os.path.splitext(os.path.basename(p))[0]
        rec.id = rec.id if rec.id else fname
        rec.id = rec.id.replace(" ", "_")
        rec.description = ""  # keep phylip clean
        # if you prefer the file name as id, uncomment next line:
        # rec.id = fname
        records.append(rec)

    # sanity: ensure all sequences have the same length
    lengths = {len(r.seq) for r in records}
    if len(lengths) != 1:
        raise RuntimeError(f"Not all sequences have the same length: observed lengths = {sorted(lengths)}. "
                           "PHYLIP requires equal lengths. Check your concatenation step.")
    print(f"Collected {len(records)} records, alignment length = {lengths.pop()}")

    # write combined multifasta
    with open(combined_fasta, "w") as fh:
        SeqIO.write(records, fh, "fasta")
    print(f"Wrote combined FASTA: {combined_fasta}")

    # write PHYLIP
    if use_phylip_classic:
        # create copies with truncated IDs to 10 chars (classical PHYLIP)
        classic_records = []
        seen = {}
        for r in records:
            short = shorten_name(r.id, 10)
            # ensure uniqueness
            if short in seen:
                # try to make unique by replacing trailing chars with index
                i = 1
                base = short[:8]
                candidate = f"{base}{i:02d}"
                while candidate in seen:
                    i += 1
                    candidate = f"{base}{i:02d}"
                short = candidate
            seen[short] = True
            newr = SeqRecord(r.seq, id=short, description="")
            classic_records.append(newr)
        with open(combined_phylip, "w") as fh:
            SeqIO.write(classic_records, fh, "phylip")
        print(f"Wrote classic PHYLIP (10-char names): {combined_phylip}")
    else:
        # relaxed PHYLIP supports longer IDs
        with open(combined_phylip, "w") as fh:
            SeqIO.write(records, fh, "phylip-relaxed")
        print(f"Wrote relaxed PHYLIP: {combined_phylip}")

if __name__ == "__main__":
    main()

