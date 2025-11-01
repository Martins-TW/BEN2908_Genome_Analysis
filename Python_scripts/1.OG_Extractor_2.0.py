#!/usr/bin/env python3
"""
OG_extractor_auto.py — versão 100% automatizada

Coloque este script no mesmo diretório onde estão:
 - SingleCopyOrthogroups.txt
 - Orthogroups.txt
 - Arquivos *.faa (um por genoma anotado)

Execução:
    python OG_extractor_auto.py

Saída:
    ./output_OGs/OG000XXXX.fasta para cada OG de cópia única.
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
        raise FileNotFoundError("SingleCopyOrthogroups.txt não encontrado no diretório atual.")
    if not os.path.exists(orthogroups):
        raise FileNotFoundError("Orthogroups.txt não encontrado no diretório atual.")

    faa_files = [os.path.join(cwd, f) for f in os.listdir(cwd) if f.endswith(".faa")]
    if not faa_files:
        raise FileNotFoundError("Nenhum arquivo .faa encontrado no diretório atual.")

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
    Cria índice {cepa_id: {gene_id: SeqRecord}}
    """
    index = defaultdict(dict)
    print("Indexando genomas...")
    for fpath in faa_files:
        print(f"  -> Lendo {os.path.basename(fpath)}")
        for rec in SeqIO.parse(fpath, "fasta"):
            try:
                parts = rec.id.split(".")
                cepa_id = int(parts[1])
                index[cepa_id][rec.id] = rec
            except Exception:
                # ignora IDs que não seguem o padrão fig|xxx.xxx.peg.xxx
                continue
    print(f"Indexação concluída. {len(index)} cepas identificadas.")
    return index

# ----------------------------------------------------------------------

def normalize_header(rec_id, strain_name=None):
    """
    Ajusta cabeçalho. Exemplo: fig|666666.448865.peg.123 -> 448865peg123
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
# bloco: escreve um .txt relacionando código numérico -> nome do arquivo .faa
from collections import defaultdict

def write_strain_code_map(faa_files, outpath):
    """
    Gera um arquivo com mapeamento: <codigo_int> -> nome do arquivo .faa (sem extensão).
    Extrai o código a partir do primeiro record.id do arquivo .faa.
    """
    code2names = defaultdict(set)

    for fpath in faa_files:
        fname = os.path.basename(fpath)
        strain_name = os.path.splitext(fname)[0]

        # extrai código do primeiro record.id do FASTA
        try:
            first_rec = next(SeqIO.parse(fpath, "fasta"))
            parts = first_rec.id.split(".")
            if len(parts) > 1 and parts[1].isdigit():
                code2names[int(parts[1])].add(strain_name)
        except Exception:
            # ignora arquivos vazios ou com erro de parsing
            continue

    # escreve o arquivo
    with open(outpath, "w", encoding="utf-8") as fh:
        fh.write("#strain_code\tstrain_name(s)\n")
        for code in sorted(code2names):
            names = ";".join(sorted(code2names[code]))
            fh.write(f"{code}\t{names}\n")
    print(f"Index de cepas e IDs armazenados em: {outpath}")

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

    print("Gerando arquivos multifasta...")
    for og in single_ogs:
        genes = og_map.get(og)
        if not genes:
            missing_log.append(f"{og}: não encontrado em Orthogroups.txt")
            continue

        seqs_to_write = []
        for gene in genes:
            try:
                parts = gene.split(".")
                cepa_id = int(parts[1])
            except Exception:
                missing_log.append(f"{og}: formato inesperado em {gene}")
                continue

            if cepa_id not in fasta_index:
                missing_log.append(f"{og}: cepa_id {cepa_id} não encontrado")
                continue

            rec = fasta_index[cepa_id].get(gene)
            if rec is None:
                missing_log.append(f"{og}: gene {gene} não encontrado")
                continue

            new_id = normalize_header(rec.id)
            seqs_to_write.append(SeqRecord(Seq(str(rec.seq)), id=new_id, description=""))

        if seqs_to_write:
            write_multifasta(og, seqs_to_write, out_dir)
            total_written += 1

    # log final
    print(f"\nProcesso concluído. {total_written} multifastas gerados em '{out_dir}'")
    if missing_log:
        with open(os.path.join(out_dir, "missing.log"), "w") as logfh:
            logfh.write("\n".join(missing_log))
        print(f"{len(missing_log)} avisos registrados em missing.log")

# ----------------------------------------------------------------------

if __name__ == "__main__":
    main()
