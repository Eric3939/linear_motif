# predict missing plddt data (1061612, 9.3%) with iupred3
# update residues_characteristics_table_blosum.tsv

import pandas as pd
import subprocess
from Bio import SeqIO
from time import time


# Read original table
c1 = time()
df = pd.read_csv("final_result/residues_characteristics_table_blosum.tsv", sep='\t')
d1 = {}     # {residueID: [plDDT, con_10, delta_con_10, con_5, delta_con_5, elm, pfam_1e-4]}
for t in df.itertuples():
    residueID, plDDT, con_10, delta_con_10, con_5, delta_con_5, elm, pfam = t[1], t[2], t[3], t[4], t[5], t[6], t[7], t[8]
    d1[residueID] = [plDDT, con_10, delta_con_10, con_5, delta_con_5, elm, pfam]
c2 = time()
print(f"Original Table read. Time Used: {c2-c1}")


d2 = {}     # iupred {residueID: disorderness}


d3 = {}     # UniProt Proteome
c1 = time()
for s in SeqIO.parse("../Data/UniProt_Human_Proteome/uniprotkb_proteome_UP000005640_AND_revi_2024_05_15.fasta/uniprotkb_proteome_UP000005640_AND_revi_2024_05_15.fasta", 'fasta'):
    p_id = s.id.split("|")[1]
    d3[p_id] = s
c2 = time()
print(f"UniProt Proteome read. Time Used: {c2-c1}")


# Update plDDT with iupred
c1 = time()
n=0
for k, v in d1.items():
    if not v[0] == 72.63121929111676:
        continue
    else:
        if k in d2:
            v[0] = d2[k]
            n+=1
        else:
            # run iupred
            # create protein fasta file
            p_id, r_i = k.split("_")
            SeqIO.write(d3[p_id], 'update_plddt_iupred_temp_files/temp_protein.fasta', 'fasta')

            # run
            try:
                c = ['python', '../Disorderness/iupred3/iupred3.py', 'update_plddt_iupred_temp_files/temp_protein.fasta', 'long']
                result = subprocess.run(c, text=True, capture_output=True)
                result = result.stdout
                # Handle error (use 'no smoothing' option)
                if result.strip() == "":
                    c = ['python', '../Disorderness/iupred3/iupred3.py', 'update_plddt_iupred_temp_files/temp_protein.fasta', 'long', '-s', 'no']
                    result = subprocess.run(c, text=True, capture_output=True)
                    result = result.stdout
                    print(f"{p_id} is predicted with 'no smoothing' option")

                # store result in d2
                with open("update_plddt_iupred_temp_files/temp_iupred_result.txt", 'w') as f:
                    f.write(result)
                df = pd.read_csv("update_plddt_iupred_temp_files/temp_iupred_result.txt", sep='\t')
                for t in df.itertuples():
                    r_i, disorder = t[1], t[3]
                    plDDT_pred = 100 - disorder * 100
                    plDDT_pred = round(plDDT_pred, 2)
                    residueID = f"{p_id}_{r_i}"
                    d2[residueID] = plDDT_pred

                # update d1
                v[0] = d2[k]
                n+=1
            except Exception as e:
                print(f"Unexpected error occured for {p_id}\n{e}")

        c2 = time()
        print(f"\rUpdate Progress: {n}/1061612. Time Used: {round(c2-c1, 0)}s", end='', flush=True)    
c2 = time()
print(f"\nplDDT updated using iupred. Time Used: {c2-c1}")

# Write
c1 = time()
with open("residues_characteristics_table_iupred.tsv", 'w') as f:
    f.write("residueID\tplDDT\tcon_10\tdelta_con_10\tcon_5\tdelta_con_5\telm\tpfam_1e-4\n")
with open("residues_characteristics_table_iupred.tsv", 'a') as f:
    for k, v in d1.items():
        f.write(f"{k}\t{v[0]}\t{v[1]}\t{v[2]}\t{v[3]}\t{v[4]}\t{v[5]}\t{v[6]}\n")
c2 = time()
print(f"File Written. Time Used: {c2-c1}")

# Verify
c1 = time()
df = pd.read_csv("residues_characteristics_table_iupred.tsv", sep='\t')
print(df)
c2 = time()
print(f"Verified. Time Used: {c2-c1}")




# Original Table read. Time Used: 86.19511079788208
# UniProt Proteome read. Time Used: 0.6546914577484131
# Update Progress: 983466/1061612. Time Used: 434.0sC0HMA1 is predicted with 'no smoothing' option
# Update Progress: 1029465/1061612. Time Used: 479.0sA0A0A0MT78 is predicted with 'no smoothing' option
# Update Progress: 1029480/1061612. Time Used: 480.0sA0A0A0MT89 is predicted with 'no smoothing' option
# Update Progress: 1029492/1061612. Time Used: 481.0sA0A0A0MTA4 is predicted with 'no smoothing' option
# Update Progress: 1029507/1061612. Time Used: 483.0sA0A0J9YX06 is predicted with 'no smoothing' option
# Update Progress: 1029522/1061612. Time Used: 485.0sA0A0J9YXA8 is predicted with 'no smoothing' option
# Update Progress: 1035864/1061612. Time Used: 501.0sP0DOY5 is predicted with 'no smoothing' option
# Update Progress: 1041076/1061612. Time Used: 504.0sP0DPI4 is predicted with 'no smoothing' option
# Update Progress: 1041080/1061612. Time Used: 505.0sP0DPR3 is predicted with 'no smoothing' option
# Update Progress: 1061612/1061612. Time Used: 549.0s
# plDDT updated using iupred. Time Used: 549.3501536846161
# File Written. Time Used: 52.77397036552429
#              residueID  plDDT   con_10  delta_con_10   con_5  delta_con_5  elm  pfam_1e-4
# 0         A0A087X1C5_1  53.74 -10.0000           0.0 -5.0000          0.0    0          0
# 1         A0A087X1C5_2  58.20 -10.0000           0.0 -5.0000          0.0    0          0
# 2         A0A087X1C5_3  68.13 -10.0000           0.0 -5.0000          0.0    0          0
# 3         A0A087X1C5_4  72.93 -10.0000           0.0 -5.0000          0.0    0          0
# 4         A0A087X1C5_5  75.49 -10.0000           0.0 -5.0000          0.0    0          0
# ...                ...    ...      ...           ...     ...          ...  ...        ...
# 11405806     Q9Y6Z2_53  55.13   1.5981           0.0  2.5556          0.0    0          0
# 11405807     Q9Y6Z2_54  54.72   1.5981           0.0  2.5556          0.0    0          0
# 11405808     Q9Y6Z2_55  54.47   1.5981           0.0  2.5556          0.0    0          0
# 11405809     Q9Y6Z2_56  52.41   1.5981           0.0  2.5556          0.0    0          0
# 11405810     Q9Y6Z2_57  50.23   1.5981           0.0  2.5556          0.0    0          0

# [11405811 rows x 8 columns]
# Verified. Time Used: 10.904101848602295










# c = ['python', 'iupred3/iupred3.py', 'iupred3/test_protein.fasta', 'long']
# result = subprocess.run(c, text=True, capture_output=True)
# # print(result.stdout)
# result = result.stdout

# with open("test_result.txt", 'w') as f:
#     f.write(result)

# df = pd.read_csv("test_result.txt", sep='\t')
# print(df)