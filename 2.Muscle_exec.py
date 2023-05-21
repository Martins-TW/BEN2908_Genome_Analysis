from Bio.Align.Applications import MuscleCommandline
import os
""" Input: OGXXXXXXX.fasta (output from OG_Extractor)
    Output: OGXXXXXXX.afa (aligned multifasta files)"""

muscle_exe= "muscle"

for filen in os.listdir("adress/adress") #change adress
    cline = MuscleCommandline(muscle_exe, input='adress/address' + filen,
                                              out= 'address/address' + filen[0:9] + ".afa") #change addresses
    os.system(str(cline))
    
print ("\ndone")
