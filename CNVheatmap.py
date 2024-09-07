import os
import re
import dash
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import ward, linkage, dendrogram, leaves_list
from scipy.spatial.distance import euclidean
import plotly.graph_objs as go
import plotly.io as pio
import plotly.figure_factory as ff
from plotly.subplots import make_subplots
import dash_bio as dashbio
from sklearn.cluster import KMeans
from dash import dcc
from dash import html


cnvMatrix = pd.DataFrame()
cloneTitle = ''
numberSample = 0

def create_pos(row):
    return f"{row['chromosome']}_{row['start']}"

def CNVadder(matrix):
    global cnvMatrix
    matrix = matrix.set_index('Position')
    cnvMatrix[numberSample] = matrix['log2']


folder_path = r'C:\Users\Eugene Cho\Desktop\BioinformaticsFiles\soma_out\\'
file_extension = '.cns'
all_files = os.listdir(folder_path)
specific_files = [f for f in all_files if f.endswith(file_extension)]
for filename in specific_files:
    file_path = os.path.join(folder_path, filename)
    data_matrix = pd.read_csv(file_path, delimiter='\t')
    data_matrix = data_matrix[data_matrix['chromosome'] != 'chrX']
    data_matrix = data_matrix[data_matrix['chromosome'] != 'chrY']
    x = filename.split("_")
    cloneTitle = x[1]
    numberSample = cloneTitle[1:3]
    if numberSample.isalpha():
        numberSample = 0
    data_matrix['SampleNum'] = numberSample
    data_matrix['Position'] = data_matrix.apply(create_pos, axis=1)
    CNVadder(data_matrix)
print(cnvMatrix)
mask = cnvMatrix.isnull().any(axis=1)
cnvMatrix = cnvMatrix[~mask]
df = pd.DataFrame(cnvMatrix)
print(df.index)
x = df.index.str.split('_').to_list()

# Extract the first column (chromosome part)
first_column = [row[0].replace('chr', '') for row in x]
second_column = [row[1] for row in x]
df['Chromosome'] = first_column
df['Pos'] = second_column
df['Pos'] = pd.to_numeric(df['Pos'])
df['Chromosome'] = pd.to_numeric(df['Chromosome'])
print(df)

cnvMatrix = df.sort_values(by=['Chromosome', 'Pos'], kind='mergesort').drop(columns=['Chromosome','Pos'])
cnvMatrix.columns = pd.to_numeric(cnvMatrix.columns)
cnvMatrix = cnvMatrix.reindex(sorted(cnvMatrix.columns), axis=1)
cnvMatrix.columns = cnvMatrix.columns.astype(str)
cnvMatrix.rename(columns={'0': 'Mixed'}, inplace=True)
print(cnvMatrix)
cnvMatrix.to_csv('outputCNV.txt', sep=' ', index=True)
custom_colorscale = [
    [0, 'blue'], [0.1, 'cyan'],  # 0 to 1
    [0.1, 'yellow'], [0.5, 'orange'],  # 1 to 5
    [0.5, 'red'], [1, 'darkred']  # 5 and higher
]
custom_colorscale1 = [
    [0, 'blue'],      # 0 value
    [0.1, 'cyan'],    # transition between 0 and 1
    [0.2, 'yellow'],  # transition between 1 and 2
    [0.5, 'orange'],  # transition between 2 and 5
    [0.75, 'red'],    # transition between 5 and higher
    [1, 'darkred']    # highest value
]


clustergram = dashbio.Clustergram(
        data=cnvMatrix.values,
        row_labels=cnvMatrix.index.tolist(),
        column_labels=cnvMatrix.columns.tolist(),
        center_values=False,
        height=900,
        width=1100,
        hidden_labels='row',
        color_map= 'spectral',
    )
app = dash.Dash(__name__)
app.layout = html.Div([
    dcc.Graph(figure=clustergram)  # Use dcc.Graph to display the clustergram
])

if __name__ == '__main__':
    app.run_server(debug=True)


# Create a Heatmap with custom color scale
heatmap = go.Heatmap(
    z=cnvMatrix.values,
    x=cnvMatrix.columns,
    y=cnvMatrix.index,
    colorscale='Spectral',
    zmin=-3,  # Minimum value
    zmax=3  # Maximum value (adjust as needed)
)

# Plot the heatmap
fig = go.Figure(data=[heatmap])
fig.update_layout(
    title='CNV Heatmap with Custom Color Scales',
    xaxis_title='Sample',
    yaxis_title='Mutations',
    yaxis_nticks=36,
    yaxis_showticklabels=False

)
fig.write_html('cnvTest.html') 
fig.show()