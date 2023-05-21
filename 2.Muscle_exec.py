from Bio.Align.Applications import MuscleCommandline
import os

print ("chewie")

muscle_exe= "muscle"

for filen in os.listdir("/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_input"):
    #print (filen)
    cline = MuscleCommandline(muscle_exe, input='/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_input/' + filen,
                                              out= '/mnt/c/Users/Tobias/Desktop/LAB/Genome_Annoucement_MT78/Phylogenetic_Analysis/Genomes_Ortho/Python_Procedures/Muscle_output/' + filen[0:9] + ".afa")
    os.system(str(cline))

#print (cline)
print ("\ndone")
