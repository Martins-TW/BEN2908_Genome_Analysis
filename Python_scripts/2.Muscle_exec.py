#!/usr/bin/env python3
"""
2.Muscle_exec.py
Aligns all OG*.fasta files in ./output_OGs using MUSCLE.

Input:  OGXXXXXXX.fasta  (output from OG_extractor_auto.py)
Output: OGXXXXXXX.afa    (aligned multifasta files)
"""

import os
import subprocess

# Path to the MUSCLE executable (adjust if needed)
MUSCLE_EXE = "/usr/bin/muscle"

# Define directories
input_dir = "output_OGs"
output_dir = os.path.join(input_dir, "Muscle_alignment")
os.makedirs(output_dir, exist_ok=True)

# Loop through each FASTA file
for filename in os.listdir(input_dir):
    if not filename.endswith(".fasta"):
        continue

    input_path = os.path.join(input_dir, filename)
    output_name = os.path.splitext(filename)[0] + ".afa"
    output_path = os.path.join(output_dir, output_name)

    print(f"Aligning {filename} -> {output_name}")

    # Run MUSCLE via subprocess 
    subprocess.run(
        [MUSCLE_EXE, "-align", input_path, "-output", output_path],
        check=True
    )

print("\nAll alignments completed successfully!")
