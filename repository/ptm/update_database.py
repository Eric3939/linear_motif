import pandas as pd
import numpy as np
from Bio import SeqIO
import re
import pickle
# import sys
# import os
# project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../../'))    # to the project root (linear_motif_cc2/)
# sys.path.append(project_root)
# from linear_motif_cc2.data.build_protein_data import Protein    

class Protein:
    def __init__(self, 
                 id: str, 
                 description: str, 
                 length: int,
                 sequence: np.ndarray,
                 pfam_domains: np.ndarray = None,
                 elm_motifs: np.ndarray = None,
                 plm: np.ndarray = None,
                 disorder: np.ndarray = None,
                 rlc: np.ndarray = None,
                 rsa: np.ndarray = None):
        """
        Initialize a Protein object with the given attributes.
        
        Parameters:
        - id: Protein identifier (string)
        - description: Protein description (string)
        - length: Length of the protein sequence (int)
        - sequence: 1D array of characters representing amino acids
        - pfam_domains: 2D array of [domain_name(str), start(int), end(int)]
        - elm_motifs: 2D array of [motif_name(str), start(int), end(int)]
        - plm: 1D array of float16 values (per-residue)
        - disorder: 1D array of float16 values (per-residue)
        - rlc: 1D array of float16 values (per-residue)
        - rsa: 1D array of float16 values (per-residue)
        """
        self.id = id
        self.description = description
        self.length = length
        self.sequence = sequence
        
        # Initialize domain and motif labels as False arrays
        self.domain_label = np.zeros(length, dtype=bool)
        self.motif_label = np.zeros(length, dtype=bool)
        
        # Set domains and update domain labels
        self.pfam_domains = pfam_domains if pfam_domains is not None else np.empty((0, 3), dtype=object)
        self._update_domain_labels()
        
        # Set motifs and update motif labels
        self.elm_motifs = elm_motifs if elm_motifs is not None else np.empty((0, 3), dtype=object)
        self._update_motif_labels()
        
        # Initialize per-residue properties
        self.plm = plm if plm is not None else np.zeros(length, dtype=np.float16)
        self.disorder = disorder if disorder is not None else np.zeros(length, dtype=np.float16)
        self.rlc = rlc if rlc is not None else np.zeros(length, dtype=np.float16)
        self.rsa = rsa if rsa is not None else np.zeros(length, dtype=np.float16)
        
        # Initialize the NEW post-translational modification sites (from BIOGRID-PTM)
        self.ptm = []   # 2d list of [PTM(str), position(int, 1-based)]

        # Validate lengths
        self._validate_lengths()
    
    def _update_domain_labels(self):
        """Update the domain_label array based on pfam_domains."""
        self.domain_label[:] = False
        for domain in self.pfam_domains:
            _, start, end = domain
            self.domain_label[start-1:end] = True
    
    def _update_motif_labels(self):
        """Update the motif_label array based on elm_motifs."""
        self.motif_label[:] = False
        for motif in self.elm_motifs:
            _, start, end = motif
            self.motif_label[start-1:end] = True
    
    def _validate_lengths(self):
        """Validate that all per-residue arrays match the protein length."""
        per_residue_arrays = [
            self.sequence,
            self.domain_label,
            self.motif_label,
            self.plm,
            self.disorder,
            self.rlc,
            self.rsa
        ]
        
        for arr in per_residue_arrays:
            if len(arr) != self.length:
                raise ValueError(f"{self.id}: Array length {len(arr)} does not match protein length {self.length}")
    
    def add_pfam_domain(self, domain_name: str, start: int, end: int):
        """Add a new PFAM domain and update labels."""
        new_domain = np.array([[domain_name, start, end]], dtype=object)
        self.pfam_domains = np.vstack([self.pfam_domains, new_domain])
        self._update_domain_labels()
    
    def add_elm_motif(self, motif_name: str, start: int, end: int):
        """Add a new ELM motif and update labels."""
        new_motif = np.array([[motif_name, start, end]], dtype=object)
        self.elm_motifs = np.vstack([self.elm_motifs, new_motif])
        self._update_motif_labels()
    
    def __repr__(self):
        return f"Protein(ID='{self.id}', Description='{self.description}', Length={self.length})"




# n=0
# for record in SeqIO.parse("/mnt/d/McGill/Xia's Lab/Data/UniProt_Human_Proteome/uniprotkb_proteome_UP000005640_AND_revi_2024_05_15.fasta/uniprotkb_proteome_UP000005640_AND_revi_2024_05_15.fasta", "fasta"):
#     try:
#         uniprot_id = record.id.split("|")[1]            # UniProt id

#         # gene_name = (record.description.split()[8]).split('=')[1]   
#         gene_name = re.search(r"GN=(\w+)", record.description).group(1)     # HGNC gene name

#     except Exception as e:
#         # print(record)
#         n += 1
#     # break



# Dictionary of Gene names mapped to UniProt ID
d_gene_prot = dict()

for record in SeqIO.parse("/mnt/d/McGill/Xia's Lab/Data/UniProt_Human_Proteome/uniprotkb_proteome_UP000005640_AND_revi_2024_05_15.fasta/uniprotkb_proteome_UP000005640_AND_revi_2024_05_15.fasta", "fasta"):
    try:
        uniprot_id = record.id.split("|")[1]            # UniProt id
        gene_name = re.search(r"GN=(\w+)", record.description).group(1)     # HGNC gene name
        
        d_gene_prot[gene_name] = uniprot_id
    except Exception as e:
        continue

with open('../../protein_database.pickle', 'rb') as f:
    protein_data = pickle.load(f)

df_ptm = pd.read_csv("BIOGRID-PTM-4.4.247.ptmtab.txt", sep='\t')

# initialize the NEW ptm attribute
for protein in protein_data.values():
    protein.ptm = []


# go through PTM to update protein database
for i, row in df_ptm.iterrows():
    
    gene_name = row['Official Symbol']
    position = row['Position']  # 1-based
    residue = row['Residue']
    sequence = row['Sequence']
    ptm = row['Post Translational Modification']
    organism = row['Organism Name']

    if gene_name == '-' or position == '-' or residue == '-' or sequence == '-' or ptm == '-':
        continue

    position = int(position)

    if organism != 'Homo sapiens':
        continue

    # Map gene name to uniprot id
    if gene_name not in d_gene_prot:
        continue
    uniprot_id = d_gene_prot[gene_name]
    protein = protein_data[uniprot_id]

    # check sequence and site: if there is a mismatch in sequence length or ptm residue, then we filter away this ptm
    if len(sequence) != len(protein.sequence):
        continue
    if residue != protein.sequence[position - 1]: # since position 1-based
        continue
    
    # update protein database
    if len(protein.ptm) == 0:
        protein.ptm = [[ptm, position]]
    else:
        if [ptm, position] in protein.ptm:
            continue               
        protein.ptm.append([ptm, position])

    protein_data[uniprot_id] = protein
        

with open('../../protein_database_1.pickle', 'wb') as f:
    pickle.dump(protein_data, f)




