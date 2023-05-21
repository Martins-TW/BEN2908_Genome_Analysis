#!/usr/bin/env python
# coding: utf-8

import os
import csv
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

OGs = [] #lista com os SingleCopyOthogroups
lst_sog= [] #lista de listas com o OGXXXXXXXXX e o grupo de genes para ser alinhado no Muscle

#print ("chewie")
#input()


with open('Results_Oct14/SingleCopyOrthogroups.txt','r') as single_ortho: #mudar endereço
    for og in single_ortho:
        og = og.split('\n')
        OGs.append(og[0])

with open('Results_Oct14/Orthogroups.csv','r') as Ortho: #mudar endereço
    Orthodata = csv.reader(Ortho)
    Orthodata = [[pos-1,row[0]] for pos,row in enumerate(Orthodata)] #lista de listas contendo a posição -1 (por causa que o Orthodata começa com OG000000000) e a os grupos da tabela Orthogroups.csv
    Orthodata.pop(0)
    
for single_og in OGs:
    for item_sog in Orthodata:
        #comparo o OGXXXXXXXXX constante na lista com apenas 1 gene por grupo com os primeiros 9 digitos (OGXXXXXXXXX) da lista com todos orthogrupos
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
        elif val_cepa == 654512:
            with open ('O157H7.faa','r') as sakai:
                for record in SeqIO.parse(sakai, "fasta"):
                    if figpeg == record.id:
                        with open ('Muscle_input/' + num_sog + '.fasta', 'a') as file_sog:
                            lst_record_id= record.id.split('.')
                            nrecord= SeqRecord(Seq (str(record.seq)),
                                                   lst_record_id[1] + lst_record_id[3],
                                                   lst_record_id[1] + lst_record_id[3],
                                                   record.description)
                            SeqIO.write(nrecord, file_sog, "fasta")
        elif val_cepa == 272792:
            with open ('PMV-1.faa','r') as pmv1:
                for record in SeqIO.parse(pmv1, "fasta"):
                    if figpeg == record.id:
                        with open ('Muscle_input/' + num_sog + '.fasta', 'a') as file_sog:
                            lst_record_id= record.id.split('.')
                            nrecord= SeqRecord(Seq (str(record.seq)),
                                                   lst_record_id[1] + lst_record_id[3],
                                                   lst_record_id[1] + lst_record_id[3],
                                                   record.description)
                            SeqIO.write(nrecord, file_sog, "fasta")
        elif val_cepa == 67883:
            with open ('EC958.faa','r') as ec958:
                for record in SeqIO.parse(ec958, "fasta"):
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
        elif val_cepa == 379313:
            with open ('IAL34.faa','r') as ial34:
                for record in SeqIO.parse(ial34, "fasta"):
                    if figpeg == record.id:
                        with open ('Muscle_input/' + num_sog + '.fasta', 'a') as file_sog:
                            lst_record_id= record.id.split('.')
                            nrecord= SeqRecord(Seq (str(record.seq)),
                                                   lst_record_id[1] + lst_record_id[3],
                                                   lst_record_id[1] + lst_record_id[3],
                                                   record.description)
                            SeqIO.write(nrecord, file_sog, "fasta")
        elif val_cepa == 380256:
            with open ('IAL35.faa','r') as ial35:
                for record in SeqIO.parse(ial35, "fasta"):
                    if figpeg == record.id:
                        with open ('Muscle_input/' + num_sog + '.fasta', 'a') as file_sog:
                            lst_record_id= record.id.split('.')
                            nrecord= SeqRecord(Seq (str(record.seq)),
                                                   lst_record_id[1] + lst_record_id[3],
                                                   lst_record_id[1] + lst_record_id[3],
                                                   record.description)
                            SeqIO.write(nrecord, file_sog, "fasta")
        elif val_cepa == 379321:
            with open ('IAL36.faa','r') as ial36:
                for record in SeqIO.parse(ial36, "fasta"):
                    if figpeg == record.id:
                        with open ('Muscle_input/' + num_sog + '.fasta', 'a') as file_sog:
                            lst_record_id= record.id.split('.')
                            nrecord= SeqRecord(Seq (str(record.seq)),
                                                   lst_record_id[1] + lst_record_id[3],
                                                   lst_record_id[1] + lst_record_id[3],
                                                   record.description)
                            SeqIO.write(nrecord, file_sog, "fasta")
        elif val_cepa == 379322:
            with open ('IAL37.faa','r') as ial37:
                for record in SeqIO.parse(ial37, "fasta"):
                    if figpeg == record.id:
                        with open ('Muscle_input/' + num_sog + '.fasta', 'a') as file_sog:
                            lst_record_id= record.id.split('.')
                            nrecord= SeqRecord(Seq (str(record.seq)),
                                                   lst_record_id[1] + lst_record_id[3],
                                                   lst_record_id[1] + lst_record_id[3],
                                                   record.description)
                            SeqIO.write(nrecord, file_sog, "fasta")
        elif val_cepa == 379323:
            with open ('IAL38.faa','r') as ial38:
                for record in SeqIO.parse(ial38, "fasta"):
                    if figpeg == record.id:
                        with open ('Muscle_input/' + num_sog + '.fasta', 'a') as file_sog:
                            lst_record_id= record.id.split('.')
                            nrecord= SeqRecord(Seq (str(record.seq)),
                                                   lst_record_id[1] + lst_record_id[3],
                                                   lst_record_id[1] + lst_record_id[3],
                                                   record.description)
                            SeqIO.write(nrecord, file_sog, "fasta")
        elif val_cepa == 379324:
            with open ('IAL39.faa','r') as ial39:
                for record in SeqIO.parse(ial39, "fasta"):
                    if figpeg == record.id:
                        with open ('Muscle_input/' + num_sog + '.fasta', 'a') as file_sog:
                            lst_record_id= record.id.split('.')
                            nrecord= SeqRecord(Seq (str(record.seq)),
                                                   lst_record_id[1] + lst_record_id[3],
                                                   lst_record_id[1] + lst_record_id[3],
                                                   record.description)
                            SeqIO.write(nrecord, file_sog, "fasta")
        elif val_cepa == 379328:
            with open ('IAL41.faa','r') as ial41:
                for record in SeqIO.parse(ial41, "fasta"):
                    if figpeg == record.id:
                        with open ('Muscle_input/' + num_sog + '.fasta', 'a') as file_sog:
                            lst_record_id= record.id.split('.')
                            nrecord= SeqRecord(Seq (str(record.seq)),
                                                   lst_record_id[1] + lst_record_id[3],
                                                   lst_record_id[1] + lst_record_id[3],
                                                   record.description)
                            SeqIO.write(nrecord, file_sog, "fasta")
        elif val_cepa == 379331:
            with open ('IAL42.faa','r') as ial42:
                for record in SeqIO.parse(ial42, "fasta"):
                    if figpeg == record.id:
                        with open ('Muscle_input/' + num_sog + '.fasta', 'a') as file_sog:
                            lst_record_id= record.id.split('.')
                            nrecord= SeqRecord(Seq (str(record.seq)),
                                                   lst_record_id[1] + lst_record_id[3],
                                                   lst_record_id[1] + lst_record_id[3],
                                                   record.description)
                            SeqIO.write(nrecord, file_sog, "fasta")
        elif val_cepa == 379338:
            with open ('IAL56.faa','r') as ial56:
                for record in SeqIO.parse(ial56, "fasta"):
                    if figpeg == record.id:
                        with open ('Muscle_input/' + num_sog + '.fasta', 'a') as file_sog:
                            lst_record_id= record.id.split('.')
                            nrecord= SeqRecord(Seq (str(record.seq)),
                                                   lst_record_id[1] + lst_record_id[3],
                                                   lst_record_id[1] + lst_record_id[3],
                                                   record.description)
                            SeqIO.write(nrecord, file_sog, "fasta")
        elif val_cepa == 475482:
            with open ('UTI89.faa','r') as uti89:
                for record in SeqIO.parse(uti89, "fasta"):
                    if figpeg == record.id:
                        with open ('Muscle_input/' + num_sog + '.fasta', 'a') as file_sog:
                            lst_record_id= record.id.split('.')
                            nrecord= SeqRecord(Seq (str(record.seq)),
                                                   lst_record_id[1] + lst_record_id[3],
                                                   lst_record_id[1] + lst_record_id[3],
                                                   record.description)
                            SeqIO.write(nrecord, file_sog, "fasta")
        elif val_cepa == 871217:
            with open ('536.faa','r') as _536: #Na verdade é 536 a linhagem
                for record in SeqIO.parse(_536, "fasta"):
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
        elif val_cepa == 871109:
            with open ('SB0258h1.faa','r') as sb0258h1:
                for record in SeqIO.parse(sb0258h1, "fasta"):
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
        elif val_cepa == 871107:
            with open ('NMBU_W05E18.faa','r') as nmbuw05e18:
                for record in SeqIO.parse(nmbuw05e18, "fasta"):
                    if figpeg == record.id:
                        with open ('Muscle_input/' + num_sog + '.fasta', 'a') as file_sog:
                            lst_record_id= record.id.split('.')
                            nrecord= SeqRecord(Seq (str(record.seq)),
                                                   lst_record_id[1] + lst_record_id[3],
                                                   lst_record_id[1] + lst_record_id[3],
                                                   record.description)
                            SeqIO.write(nrecord, file_sog, "fasta")
        elif val_cepa == 871106:
            with open ('KC-Dl-1.faa','r') as kcdl1:
                for record in SeqIO.parse(kcdl1, "fasta"):
                    if figpeg == record.id:
                        with open ('Muscle_input/' + num_sog + '.fasta', 'a') as file_sog:
                            lst_record_id= record.id.split('.')
                            nrecord= SeqRecord(Seq (str(record.seq)),
                                                   lst_record_id[1] + lst_record_id[3],
                                                   lst_record_id[1] + lst_record_id[3],
                                                   record.description)
                            SeqIO.write(nrecord, file_sog, "fasta")
        elif val_cepa == 871105:
            with open ('A1_136.faa','r') as a1136:
                for record in SeqIO.parse(a1136, "fasta"):
                    if figpeg == record.id:
                        with open ('Muscle_input/' + num_sog + '.fasta', 'a') as file_sog:
                            lst_record_id= record.id.split('.')
                            nrecord= SeqRecord(Seq (str(record.seq)),
                                                   lst_record_id[1] + lst_record_id[3],
                                                   lst_record_id[1] + lst_record_id[3],
                                                   record.description)
                            SeqIO.write(nrecord, file_sog, "fasta")
        elif val_cepa == 871104:
            with open ('631.faa','r') as _631: #Na verdade é 631 a linhagem
                for record in SeqIO.parse(_631, "fasta"):
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
        elif val_cepa == 871101:
            with open ('28Eco12.faa','r') as _28eco12: #Na verdade é 28Eco12 a linhagem
                for record in SeqIO.parse(_28eco12, "fasta"):
                    if figpeg == record.id:
                        with open ('Muscle_input/' + num_sog + '.fasta', 'a') as file_sog:
                            lst_record_id= record.id.split('.')
                            nrecord= SeqRecord(Seq (str(record.seq)),
                                                   lst_record_id[1] + lst_record_id[3],
                                                   lst_record_id[1] + lst_record_id[3],
                                                   record.description)
                            SeqIO.write(nrecord, file_sog, "fasta")
