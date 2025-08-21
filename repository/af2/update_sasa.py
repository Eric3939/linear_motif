import os
import subprocess
from Bio import SeqIO
from Bio.PDB import PDBParser
import pandas as pd

# Function to calculate average SASA across the proteome
def calculate_average_sasa(sasa_values):
    return sum(sasa_values) / len(sasa_values) if sasa_values else 0

# Function to run DSSP using subprocess and parse the output
def run_dssp(pdb_file):
    # Run DSSP and capture the output
    result = subprocess.run(["dssp", "-i", pdb_file], capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"DSSP failed to run on {pdb_file}: {result.stderr}")
    
    # Parse the DSSP output to extract SASA values
    sasa_values = []
    dssp_output = result.stdout.splitlines()
    in_residue_section = False
    for line in dssp_output:
        if line.startswith("  #  RESIDUE AA STRUCTURE"):
            in_residue_section = True
            continue
        if in_residue_section and line.strip():
            # Extract the SASA value (column 35-38 in DSSP output)
            sasa = float(line[34:38].strip())
            sasa_values.append(sasa)
    return sasa_values

# Initialize variables
p = 0
l = []
sasa_values = []

# Paths
fasta_path = "uniprotkb_proteome_UP000005640_AND_revi_2024_05_15.fasta"
alphafold_dir = "../AlphaFold/UP000005640_9606_HUMAN_v4"
output_tsv = "residues_characteristics_table_sasa.tsv"

# Parse the fasta file
for s in SeqIO.parse(fasta_path, 'fasta'):
    id = s.id.split("|")[1]
    pdb_file = f"AF-{id}-F1-model_v4.pdb"
    pdb_path = os.path.join(alphafold_dir, pdb_file)
    
    if os.path.exists(pdb_path):
        parser = PDBParser()
        structure = parser.get_structure(id, pdb_path)
        
        # Check if the lengths match
        if len(list(structure.get_residues())) == len(s):
            # Run DSSP using subprocess
            try:
                sasa_residues = run_dssp(pdb_path)
                for i, residue in enumerate(s.seq):
                    res_id = f"{id}_{i+1}"
                    plddt = list(structure.get_residues())[i]['C'].get_bfactor()
                    sasa = sasa_residues[i]
                    sasa_values.append(sasa)
                    l.append([res_id, residue, plddt, sasa])
            except Exception as e:
                print(f"Error processing {pdb_file}: {e}")
                for i, residue in enumerate(s.seq):
                    l.append([f"{id}_{i+1}", residue, "NaN", "NaN"])
        else:
            for i, residue in enumerate(s.seq):
                l.append([f"{id}_{i+1}", residue, "NaN", "NaN"])
    else:
        for i, residue in enumerate(s.seq):
            l.append([f"{id}_{i+1}", residue, "NaN", "NaN"])
    
    p += 1
    print(f"Processed {p} proteins")

# Calculate average SASA
average_sasa = calculate_average_sasa(sasa_values)

# Replace "NaN" SASA values with the average SASA
for entry in l:
    if entry[3] == "NaN":
        entry[3] = average_sasa

# Load the existing TSV file
df_existing = pd.read_csv("residues_characteristics_table_corrected.tsv", sep='\t')

# Create a new DataFrame from the list
df_new = pd.DataFrame(l, columns=['residueID', 'RLC', 'pLDDT', 'SASA'])

# Merge the new DataFrame with the existing one
df_merged = pd.merge(df_existing, df_new[['residueID', 'SASA']], on='residueID', how='left')

# Reorder columns to place SASA between embedding and elm
cols = df_merged.columns.tolist()
cols = cols[:3] + [cols[-1]] + cols[3:-1]
df_merged = df_merged[cols]

# Save the new TSV file
df_merged.to_csv(output_tsv, sep='\t', index=False)

print(f"New TSV file saved as {output_tsv}")