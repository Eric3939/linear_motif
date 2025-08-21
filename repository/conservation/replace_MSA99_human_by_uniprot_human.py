from Bio import SeqIO
import sys


n=0
m=0

d = {}  # UniProtID and EnsemblID pair
with open ("uniprot_blasted_MSA99_fmt6.txt", 'r') as f:
    for row in f:
        l_r = row.split("\t")
        uniprotID = l_r[0].split('|')[1]
        ensemblID = l_r[1].split(" ")[1]
        d[uniprotID] = ensemblID
del uniprotID, ensemblID, l_r, row


for s in SeqIO.parse("../UniProt_Human_Proteome/uniprotkb_proteome_UP000005640_AND_revi_2024_05_15.fasta/uniprotkb_proteome_UP000005640_AND_revi_2024_05_15.fasta", 'fasta'):
    uniprotid = s.id.split('|')[1]

    # UniProt Human sequence
    seqs = [s]
    
    # Find uniprot's corresponding ensembl id (blasted result)
    try:
        ensemblid = d[uniprotid]
    except KeyError:
        print(f"protein {uniprotid} has no hit in UCBC MSA99")
        m+=1
    
    # Add MSA99 animal alignments (excluding human)
    i = 0
    for msa_s in SeqIO.parse(f"concatenated_MSA99/{ensemblid}.fasta", 'fasta'):
        i+=1
        if i == 1:
            continue    # skip the first sequence (human sequence)
        seqs.append(msa_s)

    # Write into fasta file
    SeqIO.write(seqs, f"MSA99_replaced_uniprot_human/{uniprotid}.fasta", 'fasta')
    
    # n+=1
    # if n>5:
    #     break

print(f"In total, {m} uniprot proteins have no hit in UCBC MSA99")


