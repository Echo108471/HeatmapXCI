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

def create_name(row):
    return f"{row['Chromosome']}_{row['Pos']}_{row['Ref']}_{row['Alt']}"

file_path = r'C:\Users\Eugene Cho\Desktop\BioinformaticsFiles\CGI_analysis\alterations.tsv'
oncoFrame = pd.read_csv(file_path, delimiter='\t')
cgiMatrix = pd.DataFrame({
    'Chromosome': oncoFrame['chr'],
    'Pos': oncoFrame['pos'],
    'Ref': oncoFrame['ref'],
    'Alt': oncoFrame['alt'],
    'Oncogenic Prediction': oncoFrame['CGI-Oncogenic Prediction']
})
cgiMatrix['Name'] = cgiMatrix.apply(create_name, axis=1)
cgiMatrix = cgiMatrix.set_index('Name')
conditions = (
    (cgiMatrix['Oncogenic Prediction'] == 'driver (boostDM: tissue specific model)') |
    (cgiMatrix['Oncogenic Prediction'] == 'driver (oncodriveMUT)')
)
cgiMatrix = cgiMatrix[conditions]
