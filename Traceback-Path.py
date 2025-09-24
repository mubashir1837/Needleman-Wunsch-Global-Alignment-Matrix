import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Sequences
seq1 = "GATTACA"
seq2 = "GCATGCU"

# Scoring scheme
match = 1
mismatch = -1
gap = -1

# Initialize scoring matrix
n = len(seq1) + 1
m = len(seq2) + 1
score_matrix = np.zeros((m, n), dtype=int)

# Initialize first row and column
for i in range(m):
    score_matrix[i][0] = i * gap
for j in range(n):
    score_matrix[0][j] = j * gap

# Fill the matrix
for i in range(1, m):
    for j in range(1, n):
        diag = score_matrix[i-1][j-1] + (match if seq1[j-1] == seq2[i-1] else mismatch)
        delete = score_matrix[i-1][j] + gap
        insert = score_matrix[i][j-1] + gap
        score_matrix[i][j] = max(diag, delete, insert)

# Traceback path
i, j = m-1, n-1
path = [(i, j)]
while i > 0 or j > 0:
    current = score_matrix[i][j]
    if i > 0 and j > 0 and current == score_matrix[i-1][j-1] + (match if seq1[j-1] == seq2[i-1] else mismatch):
        i, j = i-1, j-1
    elif i > 0 and current == score_matrix[i-1][j] + gap:
        i -= 1
    else:
        j -= 1
    path.append((i, j))

# Convert to DataFrame
df = pd.DataFrame(score_matrix,
                  index=["-"] + list(seq2),
                  columns=["-"] + list(seq1))

# Plot heatmap
plt.figure(figsize=(8,6))
sns.heatmap(df, annot=True, fmt="d", cmap="YlGnBu", cbar=False)

# Draw traceback path (arrows)
for (y1, x1), (y2, x2) in zip(path[::-1], path[-2::-1]):
    plt.arrow(x1+0.5, y1+0.5, x2-x1, y2-y1,
              head_width=0.2, head_length=0.2,
              fc="red", ec="red")

plt.title("Needleman-Wunsch Alignment Matrix with Traceback Path")
plt.savefig("needleman_wunsch_traceback.png", dpi=300, bbox_inches="tight")
plt.show()
