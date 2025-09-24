import pandas as pd

# Example sequences
seq1 = "ACGGCTC"   # horizontal
seq2 = "ATGGCCTC"  # vertical

# Scoring scheme
match = 1
mismatch = -3
gap = -4

# Initialize matrix with dimensions (len(seq2)+1) x (len(seq1)+1)
rows = len(seq2) + 1
cols = len(seq1) + 1
matrix = [[0 for _ in range(cols)] for _ in range(rows)]

# Initialize first row and first column with gap penalties
for i in range(1, rows):
    matrix[i][0] = matrix[i-1][0] + gap
for j in range(1, cols):
    matrix[0][j] = matrix[0][j-1] + gap

# Fill matrix using Needleman-Wunsch DP algorithm
for i in range(1, rows):
    for j in range(1, cols):
        diag = matrix[i-1][j-1] + (match if seq2[i-1] == seq1[j-1] else mismatch)
        up = matrix[i-1][j] + gap
        left = matrix[i][j-1] + gap
        matrix[i][j] = max(diag, up, left)

# Convert to DataFrame for better visualization
df = pd.DataFrame(matrix)

# Add sequence labels
df.columns = ["-"] + list(seq1)
df.index = ["-"] + list(seq2)

import caas_jupyter_tools
caas_jupyter_tools.display_dataframe_to_user("Needleman-Wunsch Global Alignment Matrix", df)
