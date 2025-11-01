#!/usr/bin/env python3
"""
Move this script to the same directory of:
 - SingleCopyOrthogroups.txt
 - Orthogroups.txt
 - *.faa files

Execution:
    python 1.OG_extractor_2.0.py

Output:
    ./output_OGs/OG000XXXX.fasta for each single copy OG.
"""

import os
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# ----------------------------------------------------------------------

def find_input_files():
    cwd = os.getcwd()
    single_copy = os.path.join(cwd, "SingleCopyOrthogroups.txt")
    orthogroups = os.path.join(cwd, "Orthogroups.txt")

    if not os.path.exists(single_copy):
        raise FileNotFoundError("SingleCopyOrthogroups.txt not found in the current directory.")
    if not os.path.exists(orthogroups):
        raise FileNotFoundError("Orthogroups.txt not found in the current directory.")

    faa_files = [os.path.join(cwd, f) for f in os.listdir(cwd) if f.endswith(".faa")]
    if not faa_files:
        raise FileNotFoundError("No .faa file found on the current directory")

    return single_copy, orthogroups, faa_files

# ----------------------------------------------------------------------

def read_single_copy(single_copy_file):
    with open(single_copy_file) as fh:
        return [line.strip() for line in fh if line.strip()]

def parse_orthogroups(orthogroups_file):
    og_map = {}
    with open(orthogroups_file) as fh:
        for line in fh:
            if ':' not in line:
                continue
            og, rest = line.split(':', 1)
            genes = [g.strip() for g in rest.split() if g.strip()]
            og_map[og.strip()] = genes
    return og_map

# ----------------------------------------------------------------------

def index_fastas(faa_files):
    """
    Create Index {cepa_id: {gene_id: SeqRecord}}
    """
    index = defaultdict(dict)
    print("Indexing Genomes...")
    for fpath in faa_files:
        print(f"  -> Lendo {os.path.basename(fpath)}")
        for rec in SeqIO.parse(fpath, "fasta"):
            try:
                parts = rec.id.split(".")
                cepa_id = int(parts[1])
                index[cepa_id][rec.id] = rec
            except Exception:
                continue
    print(f"Indexing done {len(index)} strains identified")
    return index

# ----------------------------------------------------------------------

def normalize_header(rec_id, strain_name=None):
    """
    Header adjustment: fig|666666.448865.peg.123 -> 448865peg123
    """
    parts = rec_id.split(".")
    if len(parts) >= 4:
        new_id = parts[1] + parts[3]
    else:
        new_id = rec_id
    if strain_name:
        return f"{strain_name}|{new_id}"
    return new_id

def write_multifasta(og, seqs, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"{og}.fasta")
    with open(out_path, "w") as outfh:
        SeqIO.write(seqs, outfh, "fasta")

#----------------------------------------------------------------------
def write_strain_code_map(faa_files, outpath):
    """
    Generates a file relating the RAST handle to the strain name by extracting the record.id of the file:
    """
    code2names = defaultdict(set)

    for fpath in faa_files:
        fname = os.path.basename(fpath)
        strain_name = os.path.splitext(fname)[0]
        try:
            first_rec = next(SeqIO.parse(fpath, "fasta"))
            parts = first_rec.id.split(".")
            if len(parts) > 1 and parts[1].isdigit():
                code2names[int(parts[1])].add(strain_name)
        except Exception:
            continue

    # writes the file
    with open(outpath, "w", encoding="utf-8") as fh:
        fh.write("#strain_code\tstrain_name(s)\n")
        for code in sorted(code2names):
            names = ";".join(sorted(code2names[code]))
            fh.write(f"{code}\t{names}\n")
    print(f"strain_code_map.txt in: {outpath}")

# ----------------------------------------------------------------------

def main():
    single_copy_file, orthogroups_file, faa_files = find_input_files()
    single_ogs = set(read_single_copy(single_copy_file))
    og_map = parse_orthogroups(orthogroups_file)
    fasta_index = index_fastas(faa_files)
    out_dir = os.path.join(os.getcwd(), "output_OGs")
    os.makedirs(out_dir, exist_ok=True)
    write_strain_code_map(faa_files, os.path.join(out_dir, "strain_code_map.txt"))

    missing_log = []
    total_written = 0

    print("Generating multifasta files...")
    for og in single_ogs:
        genes = og_map.get(og)
        if not genes:
            missing_log.append(f"{og}: not found in Orthogroups.txt")
            continue

        seqs_to_write = []
        for gene in genes:
            try:
                parts = gene.split(".")
                cepa_id = int(parts[1])
            except Exception:
                missing_log.append(f"{og}: unexpected format in {gene}")
                continue

            if cepa_id not in fasta_index:
                missing_log.append(f"{og}: cepa_id {cepa_id} not found")
                continue

            rec = fasta_index[cepa_id].get(gene)
            if rec is None:
                missing_log.append(f"{og}: gene {gene} not found")
                continue

            new_id = normalize_header(rec.id)
            seqs_to_write.append(SeqRecord(Seq(str(rec.seq)), id=new_id, description=""))

        if seqs_to_write:
            write_multifasta(og, seqs_to_write, out_dir)
            total_written += 1

    # log final
    print(f"\nProccess done. {total_written} multifastas generated in '{out_dir}'")
    if missing_log:
        with open(os.path.join(out_dir, "missing.log"), "w") as logfh:
            logfh.write("\n".join(missing_log))
        print(f"{len(missing_log)} warns registered in missing.log")

# ----------------------------------------------------------------------

if __name__ == "__main__":
    main()
