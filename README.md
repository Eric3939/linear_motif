# linear_motif

Github repository for 
Automating Linear Motif Predictions to Map Human Signaling Networks.
Yitao (Eric) Sun, Yu (Brandon) Xia, Jasmin Coulombe-Huntington
Department of Bioengineering, McGill University

# Abstract
Short linear motifs (SLiMs) are short sequence patterns that mediate transient protein-protein interactions, often within disordered regions of proteins. SLiMs play central roles in signaling, trafficking, and post-translational regulation, but their short length and low complexity make them difficult to identify both experimentally and computationally. Since the release of motif discovery tools like MEME Suite, the availability of protein-protein interaction data (e.g., BioGRID) has increased by more than five-fold, and recent advances in machine learning offer new opportunities for large-scale, high-resolution motif discovery. Here, we present a new Gibbs sampling-based SLiM discovery method that introduces two key innovations: First, we replace the traditional position-specific scoring matrix (PSSM) with a Hidden Markov Model (HMM) to better accommodate insertions and deletions common in disordered regions; Second, we introduce biased sampling guided by pre-annotated residue-level features derived from Protein Language Models (PLMs), AlphaFold2-derived predictions (disorder, solvent accessibility), and evolutionary conservation. We evaluate our approach using the ELM database and show improved recovery of known motif instances compared to existing tools. We presented three case studies at the end to showcase potential applications of our algorithm. 

<img width="975" height="665" alt="image" src="https://github.com/user-attachments/assets/dbfec2a9-375c-4780-bd26-f3030768b41e" />


# Running the algorithm
To run the algorithm, first locate search.py in the script folder. Run the following command:

python search.py [LMBD protein]

The LMBD protein is the interacting protein that has many partners in BioGRID. The algorithm will access the BioGRID data to find all interacting partners and run the motif dicovery algorithm on them. Filters are applied as described in our paper (either low throughput or at least 2 sources). 
The parameters can be changed in the main function inside search.py


# Predicted Motif Database
The motif database will be available shortly under the results folder. This contains all the motifs found in the human proteome. They are organized into LMBD proteins. Users can either enter their desired LMBD protein as keys to look for all motifs interacting with this LMBD protein, or can use the script we provide to look for all motifs predicted with a certain protein, as contrast to the LMBD protein.  
