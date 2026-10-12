import time
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
from mpi4py import MPI

# Initialize MPI communicator, rank (process ID), and size (total process count)
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

pattern = None
X = None

# Start measuring sequential time
# t_seq_start = time.time()

# Read pattern input and CSV file on rank 0
if rank == 0:
    pattern = input("Enter pattern to search in proteins: ").upper()
    df = pd.read_csv("proteins.csv")
    X = df.iloc[:, [0, 3, 2]].to_numpy()
    # Split data into chunks
    chunks = np.array_split(X, size)
else:
    chunks = None

# Broadcast pattern to all ranks
pattern = comm.bcast(pattern, root=0)

# Compile pattern once so it is not compiled on every single iteration
pattern_re = re.compile(pattern)

# End measuring sequential time
# time_seq_total = time.time() - t_seq_start

# Synchronize processes and start timer for time execution measurement
comm.Barrier()

# Start measuring parallel time
start_time = time.time()

# Distribute chunks
local_data = comm.scatter(chunks, root=0)

# Process chunks in parallel
# Find occurences of the pattern string inside each protein sequence
local_matches = []
for prot_id, seq, hydrofob in local_data:
    matches = pattern_re.findall(seq)
    if matches:
        local_matches.append((prot_id, len(matches), float(hydrofob)))

# Locally sort matches; reduces amt of data sent to rank 0
local_matches.sort(key=lambda x: (x[1], x[2]), reverse=True)
local_top_10 = local_matches[:10]

# Gather results back to rank 0
all_matches = comm.gather(local_top_10, root=0)

# End measuring parallel time
# time_par = time.time() - start_time
# total_time = time_seq_total + time_par

# P = time_par / total_time if total_time > 0 else 0.0
# seq_fraction = 1.0 - P
# upper_bound = 1.0 / seq_fraction if seq_fraction > 0 else float('inf')

# Finish processing on rank 0
if rank == 0:
    # Flatten the list with protein matches
    protein_matches = [item for sublist in all_matches for item in sublist]

    # Sort in a descending order by matches and hydrofobs
    # If matches are equal, python moves to sorting by hydrofobs
    protein_matches.sort(key=lambda x: (x[1], x[2]), reverse=True)
    
    # Select top 10 proteins
    top_10 = protein_matches[:10]
    
    # Separate by IDs (X) and occurrences (Y) for plotting
    if top_10:
        protein_ids = [str(item[0]) for item in top_10]
        occurrences = [item[1] for item in top_10]

        print(f"Execution time: {time.time() - start_time:.4f} seconds")
        
        # Plot the bar chart
        # If several proteins exist with same number of matches, order them by the maximum Hydrofob
        plt.figure(figsize=(10, 5))
        plt.bar(protein_ids, occurrences, color='skyblue')
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
            print(f"Top protein: ID = {top_protein[0]} | Matches = {top_protein[1]} | Hydrophobicity = {top_protein[2]}")
    else:
        print("No matches found for the given pattern.")

    # Print all of Amdahl's Law Metrics
    # print(f"Amdahl's Law Metrics:")
    # print(f"Total Running Time: {total_time:.4f} seconds")
    # print(f"Parallel Time (T_par): {time_par:.4f} seconds")
    # print(f"Sequential Time (T_seq): {time_seq_total:.4f} seconds")
    # print(f"Calculated P (Parallel Fraction): {P:.4f} ({P*100:.1f}%)")
    # print(f"Calculated 1-P (Sequential Fraction): {seq_fraction:.4f} ({seq_fraction*100:.1f}%)")
    # print(f"Theoretical Upper Bound (Max Speedup): {upper_bound:.2f}x")

    # # Plot Amdahl's Curves
    # N = np.arange(1, 16)

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

# Below is the version before optimization

# # Initialize MPI communicator, rank (process ID), and size (total process count)
# comm = MPI.COMM_WORLD
# rank = comm.Get_rank()
# size = comm.Get_size()

# pattern = None
# X = None

# # Read pattern input and CSV file on rank 0
# if rank == 0:
#     pattern = input("Enter pattern to search in proteins: ").upper()
#     df = pd.read_csv("proteins.csv")
#     X = df.iloc[:, [0, 3, 2]].to_numpy()
#     # Split data into chunks
#     chunks = np.array_split(X, size)
# else:
#     chunks = None

# # Broadcast pattern to all ranks
# pattern = comm.bcast(pattern, root=0)

# # Synchronize processes and start timer for time execution measurement
# comm.Barrier()
# start_time = time.time()

# # Distribute chunks
# local_data = comm.scatter(chunks, root=0)

# # Process chunks in parallel
# # Find occurences of the pattern string inside each protein sequence
# local_matches = []
# for prot_id, seq, hydrofob in local_data:
#     matches = re.findall(pattern, seq)
#     if matches:
#         local_matches.append((prot_id, len(matches), float(hydrofob)))

# # Gather results back to rank 0
# all_matches = comm.gather(local_matches, root=0)

# # Finish processing on rank 0
# if rank == 0:
#     # Flatten the list with protein matches
#     protein_matches = [item for sublist in all_matches for item in sublist]

#     # Sort in a descending order by matches and hydrofobs
#     # If matches are equal, python moves to sorting by hydrofobs
#     protein_matches.sort(key=lambda x: (x[1], x[2]), reverse=True)
    
#     # Select top 10 proteins
#     top_10 = protein_matches[:10]
    
#     # Separate by IDs (X) and occurrences (Y) for plotting
#     protein_ids = [str(item[0]) for item in top_10]
#     occurrences = [item[1] for item in top_10]

#     print(f"Execution time: {time.time() - start_time:.4f} seconds")
    
#     # Plot the bar chart
#     # If several proteins exist with same number of matches, order them by the maximum Hydrofob
#     plt.figure(figsize=(10, 5))
#     plt.bar(protein_ids, occurrences, color='skyblue')
#     plt.xlabel('Protein ID')
#     plt.ylabel('Number of Occurrences')
#     plt.title('Top 10 Proteins with Most Pattern Matches')
#     plt.xticks(rotation=45)
#     plt.tight_layout()
#     plt.show()
    
#     # Print in the screen the id of the protein with max occurrences
#     # If several exist with same number, print the one with max. hydrofob
#     if protein_matches:
#         top_protein = protein_matches[0]
#         print(f"Top protein: ID = {top_protein[0]} | Matches = {top_protein[1]} | Hydrophobicity = {top_protein[2]}")