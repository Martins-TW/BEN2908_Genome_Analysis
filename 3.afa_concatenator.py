import os
from Bio import SeqIO
""" Input: OGXXXXXXX.afa (output from Muscle_exec)
    Output: NameOfTheStrain.fasta (concatenated .afa files of each strain)"""

for filen in os.listdir("address/address"): #change address
    afa_records =SeqIO.parse("address/address" + filen, "fasta") #change address
    for record in afa_records:
        if int(record.id[0:6]) == 555771:
            if "chi7122.fasta" in os.listdir("address/address"): #change address
                seq = record.seq
                with open("address/address/Strain.fasta", "a") as strain_afa: #change address
                    strain_afa.write(str(seq))
            else:
                with open("address/address/Strain.fasta", "a") as strain_afa: #change address
                    SeqIO.write(record, strain_afa, "fasta")
        elif int(record.id[0:6]) == 448865:
            if "APECO1.fasta" in os.listdir("address/address"): #change address
                seq = record.seq
                with open("address/address/Strain.fasta", "a") as strain_afa: #change address
                    strain_afa.write(str(seq))
            else:
                with open("address/address/Strain.fasta", "a") as strain_afa: #change address
                    SeqIO.write(record, strain_afa, "fasta")
        elif int(record.id[0:6]) == 541064:
            if "BEN2908.fasta" in os.listdir("address/address"): #change address
                seq = record.seq
                with open("address/address/Strain.fasta", "a") as strain_afa: #change address
                    strain_afa.write(str(seq))
            else:
                with open("address/address/Strain.fasta", "a") as strain_afa: #change address
                    SeqIO.write(record, strain_afa, "fasta")
        elif int(record.id[0:6]) == 290288:
            if "CFT073.fasta" in os.listdir("address/address"): #change address
                seq = record.seq
                with open("address/address/Strain.fasta", "a") as strain_afa: #change address
                    strain_afa.write(str(seq))
            else:
                with open("address/address/Strain.fasta", "a") as strain_afa: #change address
                    SeqIO.write(record, strain_afa, "fasta")
        elif int(record.id[0:6]) == 265702:
            if "IMT5155.fasta" in os.listdir("address/address"): #change address
                seq = record.seq
                with open("address/address/Strain.fasta", "a") as strain_afa: #change address
                    strain_afa.write(str(seq))
            else:
                with open("address/address/Strain.fasta", "a") as strain_afa: #change address
                    SeqIO.write(record, strain_afa, "fasta")
        elif int(record.id[0:6]) == 555769:
            if "K12.fasta" in os.listdir("address/address"): #change address
                seq = record.seq
                with open("address/address/Strain.fasta", "a") as strain_afa: #change address
                    strain_afa.write(str(seq))
            else:
                with open("address/address/Strain.fasta", "a") as strain_afa: #change address
                    SeqIO.write(record, strain_afa, "fasta")
        elif int(record.id[0:6]) == 654509:
            if "LF82.fasta" in os.listdir("address/address"): #change address
                seq = record.seq
                with open("address/address/Strain.fasta", "a") as strain_afa: #change address
                    strain_afa.write(str(seq))
            else:
                with open("address/address/Strain.fasta", "a") as strain_afa: #change address
                    SeqIO.write(record, strain_afa, "fasta")
        elif int(record.id[0:6]) == 713145:
            if "SCI07.fasta" in os.listdir("address/address"): #change address
                seq = record.seq
                with open("address/address/Strain.fasta", "a") as strain_afa: #change address
                    strain_afa.write(str(seq))
            else:
                with open("address/address/Strain.fasta", "a") as strain_afa: #change address
                    SeqIO.write(record, strain_afa, "fasta")
        elif int(record.id[0:6]) == 713160:
            if "IHE3034.fasta" in os.listdir("address/address"): #change address
                seq = record.seq
                with open("address/address/Strain.fasta", "a") as strain_afa: #change address
                    strain_afa.write(str(seq))
            else:
                with open("address/address/Strain.fasta", "a") as strain_afa: #change address
                    SeqIO.write(record, strain_afa, "fasta")
        elif int(record.id[0:6]) == 713159:
            if "NRG857c.fasta" in os.listdir("address/address"): #change address
                seq = record.seq
                with open("address/address/Strain.fasta", "a") as strain_afa: #change address
                    strain_afa.write(str(seq))
            else:
                with open("address/address/Strain.fasta", "a") as strain_afa: #change address
                    SeqIO.write(record, strain_afa, "fasta")
        elif int(record.id[0:6]) == 713146:
            if "RS218.fasta" in os.listdir("address/address"): #change address
                seq = record.seq
                with open("address/address/Strain.fasta", "a") as strain_afa: #change address
                    strain_afa.write(str(seq))
            else:
                with open("address/address/Strain.fasta", "a") as strain_afa: #change address
                    SeqIO.write(record, strain_afa, "fasta")
        elif int(record.id[0:6]) == 871113:
            if "SCU-397.fasta" in os.listdir("address/address"): #change address
                seq = record.seq
                with open("address/address/Strain.fasta", "a") as strain_afa: #change address
                    strain_afa.write(str(seq))
            else:
                with open("address/address/Strain.fasta", "a") as strain_afa: #change address
                    SeqIO.write(record, strain_afa, "fasta")
        elif int(record.id[0:6]) == 871102:
            if "78-Pyelo.fasta" in os.listdir("address/address"): #change address
                seq = record.seq
                with open("address/address/Strain.fasta", "a") as strain_afa: #change address
                    strain_afa.write(str(seq))
            else:
                with open("address/address/Strain.fasta", "a") as strain_afa: #change address
                    SeqIO.write(record, strain_afa, "fasta")
        elif int(record.id[0:6]) == 871108:
            if "PSUO2.fasta" in os.listdir("address/address"): #change address
                seq = record.seq
                with open("address/address/Strain.fasta", "a") as strain_afa: #change address
                    strain_afa.write(str(seq))
            else:
                with open("address/address/Strain.fasta", "a") as strain_afa: #change address
                    SeqIO.write(record, strain_afa, "fasta")
print ("done")
