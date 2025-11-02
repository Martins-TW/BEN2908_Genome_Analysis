#!/usr/bin/env python3
"""
Concatenador simples (sem padding/truncation e sem checagem de comprimentos).

- Extrai código da cepa do header (heurística).
- Assume que, dentro de cada OG alinhado, todas as sequências têm o mesmo comprimento
  (porque foram alinhadas pelo MUSCLE).
- Para cepas ausentes em um OG, insere gaps do comprimento do OG (aln_len).
- Não realiza truncamento, padding nem checagens adicionais das seqs.
"""
import os
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# ---------- configuração ----------
input_dir = "output_OGs/Muscle_alignment"   # onde estão os .afa
output_dir = "output_concatenated"          # onde serão escritos os FASTA concatenados
accepted_exts = (".afa", ".aln", ".fasta", ".fa")
MIN_DIGITS = 6  # ajuste se os teus códigos tiverem outro tamanho mínimo
# -----------------------------------

def longest_digit_run(s):
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
    # 1) split by '.' and check second part
    parts = rec.id.split(".")
    if len(parts) > 1:
        candidate = parts[1]
        if candidate.isdigit() and len(candidate) >= MIN_DIGITS:
            return candidate
    # 2) longest run of digits in rec.id
    run = longest_digit_run(rec.id)
    if run and len(run) >= MIN_DIGITS:
        return run
    # 3) try description
    desc = getattr(rec, "description", "") or ""
    run = longest_digit_run(desc)
    if run and len(run) >= MIN_DIGITS:
        return run
    return None

def main():
    os.makedirs(output_dir, exist_ok=True)

    files = sorted(f for f in os.listdir(input_dir) if f.lower().endswith(accepted_exts))
    if not files:
        raise FileNotFoundError(f"No alignment files found in '{input_dir}'")

    ogs = []
    strain_codes = set()

    for fname in files:
        fpath = os.path.join(input_dir, fname)
        seqs = list(SeqIO.parse(fpath, "fasta"))
        if not seqs:
            print(f"Skipping empty file: {fname}")
            continue

        # comprimento do alinhamento: assumimos que todas as seqs têm o mesmo tamanho
        aln_len = max(len(rec.seq) for rec in seqs)

        code2seq = {}
        for rec in seqs:
            code = extract_code(rec)
            if not code:
                print(f"Warning: não consegui extrair código de '{rec.id}' em {fname}; registro ignorado.")
                continue
            seq_str = str(rec.seq)   # **sem** padding/truncation
            code2seq[code] = seq_str
            strain_codes.add(code)
        ogs.append((fname, aln_len, code2seq))

    if not ogs:
        raise RuntimeError("Nenhum OG válido processado.")

    concat = {code: [] for code in sorted(strain_codes)}
    missing = {code: [] for code in concat}

    for fname, aln_len, code2seq in ogs:
        for code in concat:
            if code in code2seq:
                concat[code].append(code2seq[code])
            else:
                concat[code].append("-" * aln_len)
                missing[code].append(fname)

    # escreve FASTA por cepa (id = código)
    for code, parts in concat.items():
        seq_full = "".join(parts)
        rec = SeqRecord(Seq(seq_full), id=code, description="")
        out_path = os.path.join(output_dir, f"strain_{code}.fasta")
        with open(out_path, "w") as fh:
            SeqIO.write(rec, fh, "fasta")

    print(f"Finished. Wrote {len(concat)} concatenated FASTAs to '{output_dir}'.")

if __name__ == "__main__":
    main()

