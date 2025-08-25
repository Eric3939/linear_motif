import pandas as pd
from Bio import SeqIO

df = pd.read_csv("BIOGRID-PTM-4.4.247.ptmtab.txt", sep='\t')


# print(df['Post Translational Modification'].unique())
print(df['Post Translational Modification'].value_counts())
# print(df['Systematic Name'].value_counts())



# for col in df.columns:
#     print(f'{col}\t{df[col].iloc[0]}')      # 43190 is the first human


# Official Symbol


# n=0
# for i, row in df.iterrows():
#     if row['Organism Name'] == 'Homo sapiens':
#         print(n)
#         break
#     n+=1


    




        # sequence = row['Sequence']
        # position = row['Position']
        # residue = row['Residue']

        # print(f'site: {residue}')
        # print(f'residue found at position {position}: {sequence[position-1]}')
        # break


# Residue
# K    1035798
# S      73854
# T      16504
# Y       1911
# P         46
# R         38
# D         38
# G         30
# V         28
# N         26
# A         23
# E         16
# L          8
# F          6
# C          4
# I          4
# Q          3
# W          2
