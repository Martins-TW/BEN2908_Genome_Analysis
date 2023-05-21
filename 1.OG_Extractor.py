#!/usr/bin/env python
# coding: utf-8
""" Input: SingleCopyOrthogroups.txt (output from Orthofinder)
    Output: OGXXXXXXX.fasta (aa multifasta files with each identified single copy orthogroup for Muscle input) """

import os
import csv
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

OGs = [] #list with SingleCopyOthogroups
lst_sog= [] #list of lists with the OGXXXXXXXXX and the group of genes to be aligned in


with open('adress/adress.txt','r') as single_ortho: #change adress
    for og in single_ortho:
        og = og.split('\n')
        OGs.append(og[0])

with open('adress/adress.csv','r') as Ortho: #change adress
    Orthodata = csv.reader(Ortho)
    Orthodata = [[pos-1,row[0]] for pos,row in enumerate(Orthodata)]
    Orthodata.pop(0)
    
for single_og in OGs:
    for item_sog in Orthodata:
        if single_og == item_sog[1][0:9]:
            lst_sog.append(item_sog[1].split('\t'))

for sog in lst_sog:
   num_sog= sog.pop(0)
   for figpeg in sog:
        val_cepa= (int(figpeg.split('.')[1]))
        if val_cepa == 448865:
            with open ('APECO1.faa','r') as apeco1:
                for record in SeqIO.parse(apeco1, "fasta"):
                    if figpeg == record.id:
                        with open ('Muscle_input/' + num_sog + '.fasta', 'a') as file_sog:
                            lst_record_id= record.id.split('.')
                            nrecord= SeqRecord(Seq (str(record.seq)),
                                                   lst_record_id[1] + lst_record_id[3],
                                                   lst_record_id[1] + lst_record_id[3],
                                                   record.description)
                            SeqIO.write(nrecord, file_sog, "fasta")
        elif val_cepa == 541064:
            with open ('BEN2908.faa','r') as ben2908:
                for record in SeqIO.parse(ben2908, "fasta"):
                    if figpeg == record.id:
                        with open ('Muscle_input/' + num_sog + '.fasta', 'a') as file_sog:
                            lst_record_id= record.id.split('.')
                            nrecord= SeqRecord(Seq (str(record.seq)),
                                                   lst_record_id[1] + lst_record_id[3],
                                                   lst_record_id[1] + lst_record_id[3],
                                                   record.description)
                            SeqIO.write(nrecord, file_sog, "fasta")
        elif val_cepa == 290288:
            with open ('CFT073.faa','r') as cft073:
                for record in SeqIO.parse(cft073, "fasta"):
                    if figpeg == record.id:
                        with open ('Muscle_input/' + num_sog + '.fasta', 'a') as file_sog:
                            lst_record_id= record.id.split('.')
                            nrecord= SeqRecord(Seq (str(record.seq)),
                                                   lst_record_id[1] + lst_record_id[3],
                                                   lst_record_id[1] + lst_record_id[3],
                                                   record.description)
                            SeqIO.write(nrecord, file_sog, "fasta")
        elif val_cepa == 555771:
            with open ('chi7122.faa','r') as chi7122:
                for record in SeqIO.parse(chi7122, "fasta"):
                    if figpeg == record.id:
                        with open ('Muscle_input/' + num_sog + '.fasta', 'a') as file_sog:
                            lst_record_id= record.id.split('.')
                            nrecord= SeqRecord(Seq (str(record.seq)),
                                                   lst_record_id[1] + lst_record_id[3],
                                                   lst_record_id[1] + lst_record_id[3],
                                                   record.description)
                            SeqIO.write(nrecord, file_sog, "fasta")
        elif val_cepa == 265702:
            with open ('IMT5155.faa','r') as imt5155:
                for record in SeqIO.parse(imt5155, "fasta"):
                    if figpeg == record.id:
                        with open ('Muscle_input/' + num_sog + '.fasta', 'a') as file_sog:
                            lst_record_id= record.id.split('.')
                            nrecord= SeqRecord(Seq (str(record.seq)),
                                                   lst_record_id[1] + lst_record_id[3],
                                                   lst_record_id[1] + lst_record_id[3],
                                                   record.description)
                            SeqIO.write(nrecord, file_sog, "fasta")
        elif val_cepa == 555769:
            with open ('K12.faa','r') as k12:
                for record in SeqIO.parse(k12, "fasta"):
                    if figpeg == record.id:
                        with open ('Muscle_input/' + num_sog + '.fasta', 'a') as file_sog:
                            lst_record_id= record.id.split('.')
                            nrecord= SeqRecord(Seq (str(record.seq)),
                                                   lst_record_id[1] + lst_record_id[3],
                                                   lst_record_id[1] + lst_record_id[3],
                                                   record.description)
                            SeqIO.write(nrecord, file_sog, "fasta")
        elif val_cepa == 654509:
            with open ('LF82.faa','r') as lf82:
                for record in SeqIO.parse(lf82, "fasta"):
                    if figpeg == record.id:
                        with open ('Muscle_input/' + num_sog + '.fasta', 'a') as file_sog:
                            lst_record_id= record.id.split('.')
                            nrecord= SeqRecord(Seq (str(record.seq)),
                                                   lst_record_id[1] + lst_record_id[3],
                                                   lst_record_id[1] + lst_record_id[3],
                                                   record.description)
                            SeqIO.write(nrecord, file_sog, "fasta")
        elif val_cepa == 713160:
            with open ('IHE3034.faa','r') as ihe3034:
                for record in SeqIO.parse(ihe3034, "fasta"):
                    if figpeg == record.id:
                        with open ('Muscle_input/' + num_sog + '.fasta', 'a') as file_sog:
                            lst_record_id= record.id.split('.')
                            nrecord= SeqRecord(Seq (str(record.seq)),
                                                   lst_record_id[1] + lst_record_id[3],
                                                   lst_record_id[1] + lst_record_id[3],
                                                   record.description)
                            SeqIO.write(nrecord, file_sog, "fasta")
        elif val_cepa == 713159:
            with open ('NRG857c.faa','r') as nrg857c:
                for record in SeqIO.parse(nrg857c, "fasta"):
                    if figpeg == record.id:
                        with open ('Muscle_input/' + num_sog + '.fasta', 'a') as file_sog:
                            lst_record_id= record.id.split('.')
                            nrecord= SeqRecord(Seq (str(record.seq)),
                                                   lst_record_id[1] + lst_record_id[3],
                                                   lst_record_id[1] + lst_record_id[3],
                                                   record.description)
                            SeqIO.write(nrecord, file_sog, "fasta")
        elif val_cepa == 713146:
            with open ('RS218.faa','r') as rs218:
                for record in SeqIO.parse(rs218, "fasta"):
                    if figpeg == record.id:
                        with open ('Muscle_input/' + num_sog + '.fasta', 'a') as file_sog:
                            lst_record_id= record.id.split('.')
                            nrecord= SeqRecord(Seq (str(record.seq)),
                                                   lst_record_id[1] + lst_record_id[3],
                                                   lst_record_id[1] + lst_record_id[3],
                                                   record.description)
                            SeqIO.write(nrecord, file_sog, "fasta")
        elif val_cepa == 713145:
            with open ('SCI07.faa','r') as sci07:
                for record in SeqIO.parse(sci07, "fasta"):
                    if figpeg == record.id:
                        with open ('Muscle_input/' + num_sog + '.fasta', 'a') as file_sog:
                            lst_record_id= record.id.split('.')
                            nrecord= SeqRecord(Seq (str(record.seq)),
                                                   lst_record_id[1] + lst_record_id[3],
                                                   lst_record_id[1] + lst_record_id[3],
                                                   record.description)
                            SeqIO.write(nrecord, file_sog, "fasta")
        elif val_cepa == 871113:
            with open ('SCU-397.faa','r') as scu397:
                for record in SeqIO.parse(scu397, "fasta"):
                    if figpeg == record.id:
                        with open ('Muscle_input/' + num_sog + '.fasta', 'a') as file_sog:
                            lst_record_id= record.id.split('.')
                            nrecord= SeqRecord(Seq (str(record.seq)),
                                                   lst_record_id[1] + lst_record_id[3],
                                                   lst_record_id[1] + lst_record_id[3],
                                                   record.description)
                            SeqIO.write(nrecord, file_sog, "fasta")
        elif val_cepa == 871108:
            with open ('PSUO2.faa','r') as psuo2:
                for record in SeqIO.parse(psuo2, "fasta"):
                    if figpeg == record.id:
                        with open ('Muscle_input/' + num_sog + '.fasta', 'a') as file_sog:
                            lst_record_id= record.id.split('.')
                            nrecord= SeqRecord(Seq (str(record.seq)),
                                                   lst_record_id[1] + lst_record_id[3],
                                                   lst_record_id[1] + lst_record_id[3],
                                                   record.description)
                            SeqIO.write(nrecord, file_sog, "fasta")
        elif val_cepa == 871102:
            with open ('78-Pyelo.faa','r') as _78pyelo: #Na verdade é 78-pyelo a linhagem
                for record in SeqIO.parse(_78pyelo, "fasta"):
                    if figpeg == record.id:
                        with open ('Muscle_input/' + num_sog + '.fasta', 'a') as file_sog:
                            lst_record_id= record.id.split('.')
                            nrecord= SeqRecord(Seq (str(record.seq)),
                                                   lst_record_id[1] + lst_record_id[3],
                                                   lst_record_id[1] + lst_record_id[3],
                                                   record.description)
                            SeqIO.write(nrecord, file_sog, "fasta")
        
