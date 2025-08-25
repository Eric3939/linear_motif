# import pandas as pd

# df = pd.read_csv('BIOGRID-PTM-4.4.247.ptmtab.txt', sep='\t')
# print(df)

from Bio import SeqIO
import pandas as pd

# Load PTM data (human-only subset)
ptm_df = pd.read_csv("BIOGRID-PTM-4.4.247.ptmtab.txt", sep='\t')  # Replace with your PTM file
human_ptm = ptm_df[ptm_df["Organism Name"] == "Homo sapiens"]

# Parse UniProt FASTA to create a gene-to-UniProt mapping
uniprot_mapping = {}
for record in SeqIO.parse("/mnt/d/McGill/Xia's Lab/Data/UniProt_Human_Proteome/uniprotkb_proteome_UP000005640_AND_revi_2024_05_15.fasta/uniprotkb_proteome_UP000005640_AND_revi_2024_05_15.fasta", "fasta"):
    # Extract UniProt ID (e.g., "A0A087X1C5" from "sp|A0A087X1C5|CP2D7_HUMAN")
    uniprot_id = record.id.split("|")[1]
    
    # Extract gene name from header (e.g., "GN=CYP2D7" -> "CYP2D7")
    gene_name = None
    for part in record.description.split():
        if part.startswith("GN="):
            gene_name = part.split("=")[1]
            break
    
    if gene_name:
        uniprot_mapping[gene_name] = uniprot_id

# Map PTM data to UniProt IDs
human_ptm["UniProt_ID"] = human_ptm["Official Symbol"].map(uniprot_mapping)

# Display results
print(human_ptm)
# print(human_ptm[["PTM ID", "Official Symbol", "Entrez Gene ID", "UniProt_ID"]])



#          #PTM ID  Entrez Gene ID  BioGRID ID Systematic Name Official Symbol  ... Organism Name Has Relationships                                              Notes Source Database UniProt_ID
# 43190     676304            5371      111384               -             PML  ...  Homo sapiens              True                                                  -         BIOGRID     P29590
# 43192     676306            6613      112497               -           SUMO2  ...  Homo sapiens              True  Figure 3A; SUMO2 chains formed in the absence ...         BIOGRID     P61956
# 43196     676310            7153      113006               -           TOP2A  ...  Homo sapiens              True                                                  -         BIOGRID     P11388
# 43197     676311            7155      113008               -           TOP2B  ...  Homo sapiens              True                                                  -         BIOGRID     Q02880
# 43198     676312            8065      113743               -            CUL5  ...  Homo sapiens              True                                                  -         BIOGRID     Q93034
# ...          ...             ...         ...             ...             ...  ...           ...               ...                                                ...             ...        ...
# 1128334  1765221            6774      112651               -           STAT3  ...  Homo sapiens             False                                                  -         BIOGRID     P40763
# 1128335  1765222            6774      112651               -           STAT3  ...  Homo sapiens             False                                                  -         BIOGRID     P40763
# 1128336  1765223            7169      113022    RP11-112J3.4            TPM2  ...  Homo sapiens             False                                                  -         BIOGRID     P07951
# 1128337  1765224            7169      113022    RP11-112J3.4            TPM2  ...  Homo sapiens             False                                                  -         BIOGRID     P07951
# 1128338  1765225            7169      113022    RP11-112J3.4            TPM2  ...  Homo sapiens             False                                                  -         BIOGRID     P07951

# [999211 rows x 19 columns]
