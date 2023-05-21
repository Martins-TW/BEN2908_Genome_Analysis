import os
from Bio import SeqIO 
""" Input: NameOfTheStrain.fasta (merged multifasta file with all the outputs from afa_conatenator)
    Output: NameOfTheStrain.phy (a .phylip file for RAxML input)"""

afa_records = SeqIO.parse("address/address/File.txt","fasta") #change address
phy_records = SeqIO.write(afa_records, "address/address/File.phy", "phylip") #change address

print ("\ndone")
