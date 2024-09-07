import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import ward
from scipy.spatial.distance import euclidean
import os
vafMatrix = pd.DataFrame()
totMatrix = pd.DataFrame()
cloneTitle = "C1"

#Add the VAF data into the table
def VAFadder(matrix):
    global vafMatrix
    matrix = matrix.set_index('position')
    vafMatrix[cloneTitle] = matrix['VAF']
    totMatrix[cloneTitle] = matrix['totalCount']
    #Zero out NANs, may need to come back to this later


folder_path = r'C:\Users\Eugene Cho\Desktop\BioinformaticsFiles\ase_out\\'

# Loop through every file in the folder
for filename in os.listdir(folder_path):
    file_path = os.path.join(folder_path, filename)
    if os.path.isfile(file_path):
        data_matrix = pd.read_csv(file_path, delimiter='\t')
        x = filename.split("_")
        cloneTitle = x[1]
        data_matrix['VAF'] = data_matrix['altCount'] / data_matrix['totalCount']
        VAFadder(data_matrix)
initialSize = vafMatrix.shape[0]
mask = vafMatrix.isnull().any(axis=1)
mask1 = (totMatrix < 5).any(axis=1)
mask = mask | mask1
# List containing all positions with null VAF values
# print(nullPositions)
vafMatrix = vafMatrix[~mask]
print(len(vafMatrix))
postSize = vafMatrix.shape[0]
print("Initial Size: ", str(initialSize))
print("Post Size: ", str(postSize))
print("Rows Lost: ", str(initialSize - postSize))
# Matrix print statements for mean, vafMatrix, and the original matrix
# print(vafMatrix.mean())
# print(vafMatrix)
# print(data_matrix)
g = sns.clustermap(vafMatrix, method='ward', metric='euclidean',
                    cmap='viridis', col_cluster=True, figsize= (15, 10),
                    row_cluster=True, linewidths=.001, cbar_kws={'label': 'VAF'})
ax = g.ax_heatmap
ax.set_yticks([])
ax.tick_params(axis='y', which='major', pad=8)
plt.title("")
plt.savefig('heatmap.png', dpi=600)
# plt.show()

# sns.scatterplot(data=vafMatrix.mean())
# plt.ylim(0, 1)
# plt.title("Scatter Plot Example")
# plt.xlabel("X-axis Label")
# plt.ylabel("Y-axis Label")
# plt.show()

