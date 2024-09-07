import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import ward
from scipy.spatial.distance import euclidean
import plotly.graph_objs as go
import plotly.io as pio
import os
import re

oncoFrame = pd.read_csv('homo.txt', delimiter='\t')
oncoFrame = oncoFrame.set_index('Name')
print(oncoFrame)
# Oncoprint Designer
heatmap = go.Heatmap(
    z=oncoFrame['AF'],
    x=oncoFrame['Sample'],
    y=oncoFrame.index,
    colorscale='viridis'
)
fig = go.Figure(data=[heatmap])
fig.update_layout(
    title='Oncoprint for Mutations',
    xaxis_title='Sample',
    yaxis_title='Mutations',
    yaxis_nticks=36,
)
fig.update_yaxes(tickfont=dict(size=9))  # Set the font size here
fig.write_html('oncoTest1.html') 
fig.show()