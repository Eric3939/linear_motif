import os
from datetime import datetime
  

os.chdir("D:/McGill/Xia's Lab/ClustalO/clustal-omega-1.2.2-win64/clustal-omega-1.2.2-win64")

l = os.listdir("D:\McGill\Xia's Lab\Data\MSA_99_vertebrates\ClustalO_results")
l.append("Q8NBB2.fasta")
l.append("Q8NBR9.fasta")
l.append("Q8WZ42.fasta")

for f in os.listdir("D:\McGill\Xia's Lab\Data\MSA_99_vertebrates\MSA99_replaced_uniprot_human_exclude_empty"):
    if f not in l:
        try:
            c = datetime.now().strftime('%H:%M:%S')
            print(f"\nProcessing {f}      Time: {c}")
            os.system(f"clustalo -i ../../../Data/MSA_99_vertebrates/MSA99_replaced_uniprot_human_exclude_empty/{f} -o ../../../Data/MSA_99_vertebrates/ClustalO_results/{f} --outfmt=fasta --threads=16")
            c = datetime.now().strftime('%H:%M:%S')
            print(f"Clustal finished on {f}     Time: {c}")
        except:
            with open("ErrorLog.txt", 'a') as f2:
                f2.write(f"Error occured for {f}\n")
            c = datetime.now().strftime('%H:%M:%S')
            print(f"Error occured for {f}       Time: {c}")
