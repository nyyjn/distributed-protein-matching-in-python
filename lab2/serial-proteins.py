import time
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re

pattern = input("Enter pattern to search in proteins:").upper()

start_time = time.time()

df = pd.read_csv("proteins.csv")

# Reads the protein patterns from the file in pairs (id, sequence)
X = df.iloc[:, [0, 3]].to_numpy()

# Find occurences of the pattern string inside each protein sequence
protein_matches = []
for prot_id, seq in X:
    matches = re.findall(pattern, seq)
    if matches:
        protein_matches.append((prot_id, len(matches)))

total_occurrences = sum(count for _, count in protein_matches)

print(f"Proteins containing pattern: {len(protein_matches)}")
print(f"Total occurences: {total_occurrences}")
print(f"Execution time: {time.time() - start_time:.4f} seconds")

sorted_matches = protein_matches.sort(key=lambda protein_matches: protein_matches[1], reverse=True)

# Print a barchart of occurrences using protein id as X and number of occurrences as Y
# Represent the 10 proteins with more matches. 
plt.bar(protein_matches)
plt.xlabel(protein_matches[0])
plt.ylabel(protein_matches[1])

# If several exist with same number of matches, use them
# ordered by the maximum Hydrofob.
# 9. Print in the screen the id of the protein with max occurrences. If several
# exist with same number, print the one with max. hydrofob.
