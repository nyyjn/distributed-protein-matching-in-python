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

# Synchronize processes and start timer for time execution measurement
comm.Barrier()
start_time = time.time()

# Distribute chunks
local_data = comm.scatter(chunks, root=0)

# Process chunks in parallel
# Find occurences of the pattern string inside each protein sequence
local_matches = []
for prot_id, seq, hydrofob in local_data:
    matches = re.findall(pattern, seq)
    if matches:
        local_matches.append((prot_id, len(matches), float(hydrofob)))

# Gather results back to rank 0
all_matches = comm.gather(local_matches, root=0)

# Finish processing on rank 0
if rank == 0:
    # Flatten the list with protein matches
    protein_matches = [item for sublist in all_matches for item in sublist]

    print(f"Execution time: {time.time() - start_time:.4f} seconds")

    # Sort in a descending order by matches and hydrofobs
    # If matches are equal, python moves to sorting by hydrofobs
    protein_matches.sort(key=lambda x: (x[1], x[2]), reverse=True)
    
    # Select top 10 proteins
    top_10 = protein_matches[:10]
    
    # Separate by IDs (X) and occurrences (Y) for plotting
    protein_ids = [str(item[0]) for item in top_10]
    occurrences = [item[1] for item in top_10]
    
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