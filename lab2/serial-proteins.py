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

# Finish measuring sequential part of the code
# time_seq_1 = time.time() - start_time
# Find occurences of the pattern string inside each protein sequence

# Start measuring parallel part of the code
# par_start = time.time()

protein_matches = []
for prot_id, seq, hydrofob in X:
    matches = re.findall(pattern, seq)
    if matches:
        protein_matches.append((prot_id, len(matches), float(hydrofob)))

# Finish measuring parallel part of the code
#time_par = time.time() - par_start

# Sort in a descending order by matches and hydrofobs
# If matches are equal, python moves to sorting by hydrofobs
seq_start_2 = time.time()
protein_matches.sort(key=lambda x: (x[1], x[2]), reverse=True)

# Select top 10 proteins
top_10 = protein_matches[:10]

# Separate by IDs (X) and occurrences (Y) for plotting
protein_ids = [str(item[0]) for item in top_10]
occurences = [item[1] for item in top_10]

# Start measuring sequential part of the code
# time_seq_2 = time.time() - seq_start_2
print(f"Execution time: {time.time() - start_time:.4f} seconds")

# Finish measuring total execution time
# total_time = time.time() - start_time

# Finish measuring total sequential part of the code
# time_seq_total = time_seq_1 + time_seq_2

# Calculate P and 1 - P
# P = time_par / total_time
# seq_fraction = time_seq_total / total_time
# upper_bound = 1 / seq_fraction

# Print all of Amdahl's Law Metrics
# print(f"Amdahl's Law Metrics:")
# print(f"Total Execution Time: {total_time:.4f} seconds")
# print(f"Parallel Time (T_par): {time_par:.4f} seconds")
# print(f"Sequential Time (T_seq): {time_seq_total:.4f} seconds")
# print(f"Calculated P (Parallel Fraction): {P:.4f} ({P*100:.1f}%)")
# print(f"Calculated 1-P (Sequential Fraction): {seq_fraction:.4f} ({seq_fraction*100:.1f}%)")
# print(f"Theoretical Upper Bound (Max Speedup): {upper_bound:.2f}x")

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

# Plot Amdahl's Curves
# N = np.arange(1, 65)

# P_values = [1.0, P] 

# plt.figure(figsize=(10, 5))
# for p_val in P_values:
#     S = 1 / ((1 - p_val) + (p_val / N))
    
#     if abs(p_val - P) < 1e-5:
#         plt.plot(N, S, label=f'Measured P = {P:.2f}', linewidth=3, linestyle='solid', color='red')
#     else:
#         plt.plot(N, S, label='Ideal Linear Speedup (P = 1.0)', linewidth=1.5, linestyle=':', color='gray')

# plt.title("Amdahl's Law Speedup Curve vs. Ideal Linear Speedup", fontsize=13)
# plt.xlabel("Number of Processors (N)", fontsize=11)
# plt.ylabel("Theoretical Speedup (S)", fontsize=11)
# plt.axhline(y=upper_bound, color='green', linestyle='solid', label=f'Max Upper Bound ({upper_bound:.2f}x)')
# plt.legend(fontsize=10)
# plt.grid(True, linestyle=':', alpha=0.7)
# plt.tight_layout()
# plt.show()