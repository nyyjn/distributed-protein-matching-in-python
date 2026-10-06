import time
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re

# Get pattern input and capitalize it
pattern = input("Enter pattern to search in proteins:").upper()

# Start time measurement
start_time = time.time()

# Read proteins.csv file
df = pd.read_csv("proteins.csv")

# Read the protein patterns from the file in pairs (id, sequence, hydrofob)
X = df.iloc[:, [0, 3, 2]].to_numpy()

# Find occurences of the pattern string inside each protein sequence
protein_matches = []
for prot_id, seq, hydrofob in X:
    matches = re.findall(pattern, seq)
    if matches:
        protein_matches.append((prot_id, len(matches), float(hydrofob)))

print(f"Execution time: {time.time() - start_time:.4f} seconds")

# Sort in a descending order by matches and hydrofobs
# If matches are equal, python moves to sorting by hydrofobs
protein_matches.sort(key=lambda x: (x[1], x[2]), reverse=True)

# Select top 10 proteins
top_10 = protein_matches[:10]

# Separate by IDs (X) and occurrences (Y) for plotting
protein_ids = [str(item[0]) for item in top_10]
occurences = [item[1] for item in top_10]

# Plot the bar chart
# If several proteins exist with same number of matches, order them by the maximum Hydrofob
plt.figure(figsize=(10, 5))
plt.bar(protein_ids, occurences, color='skyblue')
plt.xlabel('Protein ID')
plt.ylabel('Number of Occurrences')
plt.title('Top 10 Proteins with Most Pattern Matches')
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()

# Print in the screen the id of the protein with max occurrences
# If several exist with same number, print the one with max. hydrofob
if protein_matches:
    top_protein = protein_matches[0]
    print(f"\nTop protein: ID = {top_protein[0]} | Matches = {top_protein[1]} | Hydrophobicity = {top_protein[2]}")

