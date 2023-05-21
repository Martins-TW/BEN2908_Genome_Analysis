import os
from Bio import SeqIO

for filen in os.listdir("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output"):
    afa_records =SeqIO.parse("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output/" + filen, "fasta")
    for record in afa_records:
        if int(record.id[0:6]) == 555771:
            #print ("ok")
            if "chi7122.fasta" in os.listdir("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat"):
                seq = record.seq
                #input()
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/chi7122.fasta", "a") as strain_afa:
                    strain_afa.write(str(seq))
            else:
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/chi7122.fasta", "a") as strain_afa:
                    SeqIO.write(record, strain_afa, "fasta")
        elif int(record.id[0:6]) == 448865:
            #print ("ok")
            if "APECO1.fasta" in os.listdir("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat"):
                seq = record.seq
                #input()
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/APECO1.fasta", "a") as strain_afa:
                    strain_afa.write(str(seq))
            else:
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/APECO1.fasta", "a") as strain_afa:
                    SeqIO.write(record, strain_afa, "fasta")
        elif int(record.id[0:6]) == 541064:
            #print ("ok")
            if "BEN2908.fasta" in os.listdir("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat"):
                seq = record.seq
                #input()
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/BEN2908.fasta", "a") as strain_afa:
                    strain_afa.write(str(seq))
            else:
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/BEN2908.fasta", "a") as strain_afa:
                    SeqIO.write(record, strain_afa, "fasta")
        elif int(record.id[0:6]) == 290288:
            #print ("ok")
            if "CFT073.fasta" in os.listdir("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat"):
                seq = record.seq
                #input()
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/CFT073.fasta", "a") as strain_afa:
                    strain_afa.write(str(seq))
            else:
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/CFT073.fasta", "a") as strain_afa:
                    SeqIO.write(record, strain_afa, "fasta")
        elif int(record.id[0:6]) == 265702:
            #print ("ok")
            if "IMT5155.fasta" in os.listdir("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat"):
                seq = record.seq
                #input()
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/IMT5155.fasta", "a") as strain_afa:
                    strain_afa.write(str(seq))
            else:
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/IMT5155.fasta", "a") as strain_afa:
                    SeqIO.write(record, strain_afa, "fasta")
        elif int(record.id[0:6]) == 555769:
            #print ("ok")
            if "K12.fasta" in os.listdir("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat"):
                seq = record.seq
                #input()
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/K12.fasta", "a") as strain_afa:
                    strain_afa.write(str(seq))
            else:
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/K12.fasta", "a") as strain_afa:
                    SeqIO.write(record, strain_afa, "fasta")
        elif int(record.id[0:6]) == 654509:
            #print ("ok")
            if "LF82.fasta" in os.listdir("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat"):
                seq = record.seq
                #input()
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/LF82.fasta", "a") as strain_afa:
                    strain_afa.write(str(seq))
            else:
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/LF82.fasta", "a") as strain_afa:
                    SeqIO.write(record, strain_afa, "fasta")
        elif int(record.id[0:6]) == 654512:
            #print ("ok")
            if "Sakai.fasta" in os.listdir("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat"):
                seq = record.seq
                #input()
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/Sakai.fasta", "a") as strain_afa:
                    strain_afa.write(str(seq))
            else:
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/Sakai.fasta", "a") as strain_afa:
                    SeqIO.write(record, strain_afa, "fasta")
        elif int(record.id[0:6]) == 272792:
            #print ("ok")
            if "PMV1.fasta" in os.listdir("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat"):
                seq = record.seq
                #input()
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/PMV1.fasta", "a") as strain_afa:
                    strain_afa.write(str(seq))
            else:
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/PMV1.fasta", "a") as strain_afa:
                    SeqIO.write(record, strain_afa, "fasta")
        elif int(record.id[0:6]) == 713145:
            #print ("ok")
            if "SCI07.fasta" in os.listdir("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat"):
                seq = record.seq
                #input()
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/SCI07.fasta", "a") as strain_afa:
                    strain_afa.write(str(seq))
            else:
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/SCI07.fasta", "a") as strain_afa:
                    SeqIO.write(record, strain_afa, "fasta")
        elif int(record.id[0:6]) == 713160:
            #print ("ok")
            if "IHE3034.fasta" in os.listdir("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat"):
                seq = record.seq
                #input()
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/IHE3034.fasta", "a") as strain_afa:
                    strain_afa.write(str(seq))
            else:
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/IHE3034.fasta", "a") as strain_afa:
                    SeqIO.write(record, strain_afa, "fasta")
        elif int(record.id[0:6]) == 713159:
            #print ("ok")
            if "NRG857c.fasta" in os.listdir("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat"):
                seq = record.seq
                #input()
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/NRG857c.fasta", "a") as strain_afa:
                    strain_afa.write(str(seq))
            else:
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/NRG857c.fasta", "a") as strain_afa:
                    SeqIO.write(record, strain_afa, "fasta")
        elif int(record.id[0:6]) == 713146:
            #print ("ok")
            if "RS218.fasta" in os.listdir("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat"):
                seq = record.seq
                #input()
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/RS218.fasta", "a") as strain_afa:
                    strain_afa.write(str(seq))
            else:
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/RS218.fasta", "a") as strain_afa:
                    SeqIO.write(record, strain_afa, "fasta")
        elif int(record.id[0:5]) == 67883:
            #print ("ok")
            if "EC958.fasta" in os.listdir("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat"):
                seq = record.seq
                #input()
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/EC958.fasta", "a") as strain_afa:
                    strain_afa.write(str(seq))
            else:
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/EC958.fasta", "a") as strain_afa:
                    SeqIO.write(record, strain_afa, "fasta")
        elif int(record.id[0:6]) == 379313:
            #print ("ok")
            if "IAL34.fasta" in os.listdir("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat"):
                seq = record.seq
                #input()
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/IAL34.fasta", "a") as strain_afa:
                    strain_afa.write(str(seq))
            else:
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/IAL34.fasta", "a") as strain_afa:
                    SeqIO.write(record, strain_afa, "fasta")
        elif int(record.id[0:6]) == 380256:
            #print ("ok")
            if "IAL35.fasta" in os.listdir("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat"):
                seq = record.seq
                #input()
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/IAL35.fasta", "a") as strain_afa:
                    strain_afa.write(str(seq))
            else:
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/IAL35.fasta", "a") as strain_afa:
                    SeqIO.write(record, strain_afa, "fasta")
        elif int(record.id[0:6]) == 379321:
            #print ("ok")
            if "IAL36.fasta" in os.listdir("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat"):
                seq = record.seq
                #input()
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/IAL36.fasta", "a") as strain_afa:
                    strain_afa.write(str(seq))
            else:
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/IAL36.fasta", "a") as strain_afa:
                    SeqIO.write(record, strain_afa, "fasta")
        elif int(record.id[0:6]) == 379322:
            #print ("ok")
            if "IAL37.fasta" in os.listdir("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat"):
                seq = record.seq
                #input()
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/IAL37.fasta", "a") as strain_afa:
                    strain_afa.write(str(seq))
            else:
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/IAL37.fasta", "a") as strain_afa:
                    SeqIO.write(record, strain_afa, "fasta")
        elif int(record.id[0:6]) == 379323:
            #print ("ok")
            if "IAL38.fasta" in os.listdir("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat"):
                seq = record.seq
                #input()
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/IAL38.fasta", "a") as strain_afa:
                    strain_afa.write(str(seq))
            else:
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/IAL38.fasta", "a") as strain_afa:
                    SeqIO.write(record, strain_afa, "fasta")
        elif int(record.id[0:6]) == 379324:
            #print ("ok")
            if "IAL39.fasta" in os.listdir("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat"):
                seq = record.seq
                #input()
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/IAL39.fasta", "a") as strain_afa:
                    strain_afa.write(str(seq))
            else:
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/IAL39.fasta", "a") as strain_afa:
                    SeqIO.write(record, strain_afa, "fasta")
        elif int(record.id[0:6]) == 379328:
            #print ("ok")
            if "IAL41.fasta" in os.listdir("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat"):
                seq = record.seq
                #input()
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/IAL41.fasta", "a") as strain_afa:
                    strain_afa.write(str(seq))
            else:
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/IAL41.fasta", "a") as strain_afa:
                    SeqIO.write(record, strain_afa, "fasta")
        elif int(record.id[0:6]) == 379331:
            #print ("ok")
            if "IAL42.fasta" in os.listdir("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat"):
                seq = record.seq
                #input()
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/IAL42.fasta", "a") as strain_afa:
                    strain_afa.write(str(seq))
            else:
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/IAL42.fasta", "a") as strain_afa:
                    SeqIO.write(record, strain_afa, "fasta")
        elif int(record.id[0:6]) == 379338:
            #print ("ok")
            if "IAL56.fasta" in os.listdir("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat"):
                seq = record.seq
                #input()
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/IAL56.fasta", "a") as strain_afa:
                    strain_afa.write(str(seq))
            else:
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/IAL56.fasta", "a") as strain_afa:
                    SeqIO.write(record, strain_afa, "fasta")
        elif int(record.id[0:6]) == 475482:
            #print ("ok")
            if "UTI89.fasta" in os.listdir("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat"):
                seq = record.seq
                #input()
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/UTI89.fasta", "a") as strain_afa:
                    strain_afa.write(str(seq))
            else:
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/UTI89.fasta", "a") as strain_afa:
                    SeqIO.write(record, strain_afa, "fasta")
        elif int(record.id[0:6]) == 871113:
            #print ("ok")
            if "SCU-397.fasta" in os.listdir("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat"):
                seq = record.seq
                #input()
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/SCU-397.fasta", "a") as strain_afa:
                    strain_afa.write(str(seq))
            else:
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/SCU-397.fasta", "a") as strain_afa:
                    SeqIO.write(record, strain_afa, "fasta")
        elif int(record.id[0:6]) == 871102:
            #print ("ok")
            if "78-Pyelo.fasta" in os.listdir("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat"):
                seq = record.seq
                #input()
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/78-Pyelo.fasta", "a") as strain_afa:
                    strain_afa.write(str(seq))
            else:
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/78-Pyelo.fasta", "a") as strain_afa:
                    SeqIO.write(record, strain_afa, "fasta")
        elif int(record.id[0:6]) == 871108:
            #print ("ok")
            if "PSUO2.fasta" in os.listdir("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat"):
                seq = record.seq
                #input()
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/PSUO2.fasta", "a") as strain_afa:
                    strain_afa.write(str(seq))
            else:
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/PSUO2.fasta", "a") as strain_afa:
                    SeqIO.write(record, strain_afa, "fasta")
        elif int(record.id[0:6]) == 871217:
            #print ("ok")
            if "536.fasta" in os.listdir("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat"):
                seq = record.seq
                #input()
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/536.fasta", "a") as strain_afa:
                    strain_afa.write(str(seq))
            else:
                with open("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output_concat/536.fasta", "a") as strain_afa:
                    SeqIO.write(record, strain_afa, "fasta")
        
print ("done")
