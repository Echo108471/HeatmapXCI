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


file_path = r'C:\Users\Eugene Cho\Desktop\BioinformaticsFiles\MeTel-master\GENIE\data_mutations_extended.txt'
data_mut = pd.read_csv(file_path, delimiter='\t')
file_path1 = r'C:\Users\Eugene Cho\Desktop\BioinformaticsFiles\MeTel-master\GENIE\data_STAD.txt'
STAD = pd.read_csv(file_path1, delimiter='\t')
# print(STAD)
combine = data_mut[data_mut['Tumor_Sample_Barcode'].isin(STAD['Sample-Identifier'])]
data_cleaned = pd.DataFrame({
    'Symbol': combine['Hugo_Symbol'],
    'HGVSc': combine['HGVSc'],
    'Barcode': combine['Tumor_Sample_Barcode'],
})

# data_cleaned['HGVSc'] = data_cleaned['HGVSc'].astype(str).apply(lambda x: x.split(':', 1)[-1] if ':' in x else x)
data_cleaned['Combined'] = data_cleaned['Symbol'] + '\t' + data_cleaned['HGVSc']
distribution = data_cleaned['Combined'].value_counts()
distribution.to_csv('mutDist.txt', sep='\t', index=True)
print(distribution)
data_cleaned.to_csv('STADMut.txt', sep='\t', index=True)
# print(combine)