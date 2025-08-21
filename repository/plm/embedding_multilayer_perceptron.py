import h5py
import numpy as np
import pandas as pd
# import csv
import random
# from sklearn.linear_model import LogisticRegression
# from sklearn.model_selection import train_test_split
# from sklearn.metrics import confusion_matrix#, classification_report, ConfusionMatrixDisplay
# # import matplotlib.pyplot as plt



# Data Sampling

# n_negative_labels = 17154   # number of negatively labelled embeddings I want to learn from

# Store all motifs' coordinates in a set
s = set()

df1 = pd.read_csv("../../Data/ELM/elm_corrected.tsv", sep='\t')
df1 = df1[['UniProtID', 'Start', 'End']]
for t in df1.itertuples():
    p_id = t[1]
    start = t[2]
    end = t[3]
    start = int(start)
    end = int(end)
    for i in range(start, end + 1):
        s.add(tuple([p_id, i]))
print("Stored all motifs coordinates in a set")

n_negative_labels = len(s)*3
print(n_negative_labels)

# Write motifs (positive label)'s embeddings
l = []  # output
with h5py.File("../../Data/Embeddings/per-residue.h5", 'r') as f1:
    for t in s:
        p_id, r_id = t
        try:
            e = list(f1[p_id][r_id - 1])
            l.append(e)
        except:
            print(f"protein residue pair {p_id} {r_id} not in UniProt embeddings")

print(len(l))
print(len(l[0]))
print(len(l[1]))
print(len(l[100]))
print(type(l[1]))
print("Wrote motifs embeddings")
    
# Ramdomly select non-motifs (negative label)
with h5py.File("../../Data/Embeddings/per-residue.h5", 'r') as f2:
    # stores id, h5 object pairs
    d = {}
    for id, e in f2.items():
        d[id] = e
    l_h5 = list(d.items())
    del d
    print("Stored id, h5 object pairss")
    
    n = 0
    while n < n_negative_labels:
        pair = random.choice(l_h5)
        id = pair[0]
        p_e = np.array(pair[1])                 # protein embeddings
        r_i = random.randint(0, len(p_e) - 1)   # residue 
        
        if (id, r_i) not in s:
            r_e = p_e[r_i]                          # residue embeddings
            r_e = r_e.reshape(1, -1)                # reshape to a horizontal row
            r_e = r_e.tolist()
            l.append(r_e[0])
            
            if n ==1:
                print(r_e)
                print(r_e[0])
                print(r_e[-1])
                print(type(r_e))
                print(len(r_e))
                print(len(r_e[0]))

        n += 1
        # print("negative label:", n)


print(len(l))
print(len(l[-1]))
print(type(l[-1]))
del s
# df2 = pd.DataFrame(l, columns=list(range(1, 1025)))
# del l

# df2.to_csv("learning_embeddings.tsv", sep="\t", index=False)

# with open('learning_embeddings.tsv', 'w', newline='') as file:
#     writer = csv.writer(file, delimiter='\t')
#     writer.writerows(l)

X = np.array(l, dtype=np.float32)
from sys import getsizeof
print(f"Array size: {getsizeof(l)}")
memory_usage = getsizeof(l) + sum(getsizeof(sublist) for sublist in l) + sum(getsizeof(item) for sublist in l for item in sublist)
print(f"List size (sublists): {memory_usage}")

print(f"numpy memory: {X.nbytes}")
del l
y1 = np.ones(n_negative_labels)
y2 = np.zeros(len(X) - n_negative_labels)
y = np.concatenate((y1, y2))
del y1, y2

from sklearn.neural_network import MLPClassifier
# from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix#, classification_report, ConfusionMatrixDisplay
# import matplotlib.pyplot as plt


X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

mlp = MLPClassifier(hidden_layer_sizes=(512, 256, 128), max_iter=10000, random_state=42, verbose=1, solver='sgd', learning_rate='adaptive')

mlp.fit(X_train, y_train)
print(mlp.score(X_test, y_test))

y_pred = mlp.predict(X_test)

# # Try Logistic Regression
# lr = LogisticRegression(max_iter=10000, verbose=1)
# lr.fit(X_train, y_train)
# print(lr.score(X_test, y_test))
# y_pred = lr.predict(X_test)

cm = confusion_matrix(y_test, y_pred)
print(cm)
tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel()
print(tn, fp, fn, tp)

print("precision=", tp/(tp+fp))
print("sensitivity(recall)=", tp/(tp+fn))
print("specificity=", tn/(tn+fp))