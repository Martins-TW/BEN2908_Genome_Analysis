import os
from Bio import SeqIO 

#muscle_exe= "muscle"

afa_records = SeqIO.parse("Single_gene/fimH_tree_strains_protein_sequences.txt","fasta")
phy_records = SeqIO.write(afa_records, "Raxml_input_SG/fimH_tree_strains_protein_sequences.phy", "phylip")


#C:\Users\Tobias\Desktop\LAB\Genome_Annoucement_MT78\Phylogenetic_Analysis\Genomes_Ortho\Python_Procedures\Single_gene

print ("\ndone")
