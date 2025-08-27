# linear_motif

**Automating Linear Motif Predictions to Map Human Signaling Networks**
Yitao (Eric) Sun, Yu (Brandon) Xia, Jasmin Coulombe-Huntington
Department of Bioengineering, McGill University


## Abstract
Short linear motifs (SLiMs) are short sequence patterns that mediate transient protein-protein interactions, often within disordered regions of proteins. SLiMs play central roles in signaling, trafficking, and post-translational regulation, but their short length and low complexity make them difficult to identify both experimentally and computationally. Since the release of motif discovery tools like MEME Suite, the availability of protein-protein interaction data (e.g., BioGRID) has increased by more than five-fold, and recent advances in machine learning offer new opportunities for large-scale, high-resolution motif discovery. Here, we present a new Gibbs sampling-based SLiM discovery method that introduces two key innovations: First, we replace the traditional position-specific scoring matrix (PSSM) with a Hidden Markov Model (HMM) to better accommodate insertions and deletions common in disordered regions; Second, we introduce biased sampling guided by pre-annotated residue-level features derived from Protein Language Models (PLMs), AlphaFold2-derived predictions (disorder, solvent accessibility), and evolutionary conservation. We evaluate our approach using the ELM database and show improved recovery of known motif instances compared to existing tools. We presented three case studies at the end to showcase potential applications of our algorithm. 

<img width="1200" height="800" alt="image" src="https://github.com/user-attachments/assets/debcdcef-579a-43c2-a7ad-ec1cf1b0d472" />


## Running the algorithm
Before running the algorithm:
We call a hub protein (a protein that interacts with several other proteins) a linear motif binding domain protein (LMBD protein). The algorithm (method 1 and 3) finds all the interacting proteins of a LMBD protein and search for a motif in these interacting motifs. We call these interacting proteins a LMBD protein network. The protein-protein interactions are obtained from BioGRID, wiht a filter applied asdescribed in our paper (either low throughput or at least 2 sources). 

### Step 1 - Clone repo
Download this Github repository in user's local device. 
<code>
  git clone https://github.com/Eric3939/linear_motif.git
  cd linear_motif
</code>

### Step 2 - Download data
Download these files from Google Drive and place under `data/`:
https://drive.google.com/drive/folders/1472iWG8U6g5XaJBz2bdI_UI-kOFbpU-n?usp=sharing
- protein_database_1.pickle
- biogrid_net.gpickle
The first pickle file contains all the human proteins we included in our study, with their feature scores (PLM, disorder, solvent accessibility, conservation) pre-annotated.
The second gpickel file contains all the human PPIs in BioGRID, their number of citations and throughput are annotated. 
Download the two files below in https://drive.google.com/drive/folders/1472iWG8U6g5XaJBz2bdI_UI-kOFbpU-n?usp=sharing
protein_database_1.pickle
This contains all the human proteins we included in our study, with their feature scores (PLM, disorder, solvent accessibility, conservation) pre-annotated.

biogrid_net.gpickle
This contains all the human PPIs in BioGRID, their number of citations and throughput are annotated. 

Save the two files directly under the data/ folder. 

Step 3
There are a few ways to run our algorithm. You can either run it on one single LMBD protein network, or on all the LMBD protein networks in BioGRID.
1. Single LMBD protein network with proteins we identified in BioGRID
Locate to script/ folder.
Run code:
python search.py [LMBD protein] [output directory]

In this method, user is asked to only give a LMBD protein. The algorithm will find all the proteins that interact with the LMBD protein from BioGRID and serch motifs on these interacting proteins. All proteins using this method are pre-annotated. This is the easiest way to use our algorithm, however, it does not offer fluxibility to define what proteins to input. 

2. User defined proteins
Locate to script/ folder.
Create a text file that contains all the protein user wants to search for motifs. One protein's UniProt ID per line.
Run code:
python run_users_proteins.py [proteins_file_path]

In this method, user is asked to input the proteins that they believe might contain a common motif. This is a common way that most linear motif discovery algorithms work. The algorithm retrieves the pre-annotated proteins from our database. Please note that our algorithm currently does not support proteins that are not yet recorded in our database.

3. All LMBD protein networks in human proteome.
Locate to script/ folder.
Run code:
python run_proteome.py [biogrid_path] [output_folder]

Running this code will recreate our study in the paper. User can also explore running the algorithm in other parameters by modifying the search.py script.

If user is using a High Performance Computing (HPC) system with SLURM scheduler, we also made a tool that submit motif searching jobs automatically (user might need to edit the script for it to be compatible with their own pipeline):
Locate to script/ folder.
Run code:
python submit_slurm.py [biogrid_path] [results_folder] [dataframe_path]


Step 4
The results will be saved inside the output folder. Each search on a LMBD protein network outputs one pickle file, which stores the run result of that network. This means that a full proteome run (e.g., method 3) will output multiple pickle files.
User can use the script result_table.py under script/ folder to convert all pickle files into one single table listing all the discovered motifs.



# Predicted Motif Database
We ran our algorithm on all human proteins as described in our paper. The results (221840 predicted motifs in human proteome) are listed in the table results_table.csv.

# Paper repository
All other codes used in our paper are in the repostory/ folder, especially the codes we used to calculate the four feature scores (PLM, disorder, solvent accessibility, conservation).
