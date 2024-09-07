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



somMatrix = pd.DataFrame()
cloneTitle = ''
numberSample = 0

# Function to create the names in the format of Chr_Pos_Ref_Alt
def create_name(row):
    return f"{row['Chromosome']}_{row['StartPos']}_{row['Ref']}_{row['Alt']}"

def create_name1(row):
    return f"{row['Chromosome']} {row['StartPos']} {row['Ref']} {row['Alt']}"


# # Add the somatic mutation data into the table
# def SOMadder(matrix):
#     global somMatrix
#     new_data = pd.DataFrame({
#         'Sample': [cloneTitle] * len(matrix),
#         'Chromosome': matrix['Chr'],
#         'StartPos': matrix['Start'],
#         'EndPos': matrix['End'],
#         'Ref': matrix['Ref'],
#         'Alt': matrix['Alt'],
#         'OtherInfo9': matrix['Otherinfo9'],
#         'OtherInfo10': matrix['Otherinfo10'],
#         'SampleNum': numberSample
#     })
#     # print(new_data)
#     somMatrix = somMatrix._append(new_data, ignore_index=True)


# # Iteratively finds all of the annotated sample data files in the folder
# folder_path = r'C:\Users\Eugene Cho\Desktop\BioinformaticsFiles\HC\\'
# file_extension = '.txt'
# all_files = os.listdir(folder_path)
# specific_files = [f for f in all_files if f.endswith(file_extension)]
# # print(specific_files)
# for filename in specific_files:
#     file_path = os.path.join(folder_path, filename)
#     data_matrix = pd.read_csv(file_path, delimiter='\t')
#     data_matrix = data_matrix[data_matrix['Chr'].str.match(r'chr[1-9]$|chr1[0-9]$|chr2[0-2]$')]
#     data_matrix = data_matrix[data_matrix['Func_refGene'].str.match('exonic')]
#     data_matrix = data_matrix[data_matrix['Otherinfo7'].str.match('PASS')]
#     x = filename.split("_")
#     cloneTitle = x[1]
#     print(cloneTitle)
#     print(data_matrix.shape)
#     numberSample = cloneTitle[1:3]
#     if numberSample.isalpha():
#         numberSample = 0
#     data_matrix['SampleNum'] = numberSample
#     SOMadder(data_matrix)

# # Sorting samples
# somMatrix['SampleNum'] = pd.to_numeric(somMatrix['SampleNum'])
# somMatrix = somMatrix.sort_values(by=['SampleNum'], kind='mergesort').drop(columns='SampleNum')
# somMatrix['Ref'] = somMatrix['Ref'].apply(lambda x: x[:4])
# somMatrix['Alt'] = somMatrix['Alt'].apply(lambda x: x[:4])
# print(somMatrix.shape)
# somMatrix.to_csv('output.txt', sep='|', index=True)



# make a second dataframe to organize the data, splitting the values
data_matrix = pd.read_csv('output.txt', delimiter='|')
print(data_matrix)
oncoFrame = pd.DataFrame(columns=['GT', 'AD', 'DP', 'GQ', 'PL'])
split_data = data_matrix['OtherInfo10'].str.split(':', expand=True)
split_data = split_data.filter([0, 1, 2, 3, 4])
split_data.columns = ['GT', 'AD', 'DP', 'GQ', 'PL']
split_data['AD'] = split_data['AD'].str.split(',').str[1]
split_data = split_data[split_data['GT'].str.match(r'^(0/1|0\|1|1/1|1\|1|1/0|1\|0)$')]
split_data['DP'] = pd.to_numeric(split_data['DP'])
split_data = split_data[split_data['DP'] >= 10]
split_data['AF'] = split_data['AD'].astype(float) / split_data['DP'].astype(float)
split_data['Name'] = data_matrix.apply(create_name, axis=1)
split_data['Sample'] = data_matrix['Sample']
split_data['Name1'] = data_matrix.apply(create_name1, axis=1)
print(split_data)
oncoFrame = pd.concat([oncoFrame, split_data], ignore_index=True)
oncoFrame.to_csv('output1.txt', sep='\t', index=False)
oncoFrame['Name1'].to_csv('output2.txt', sep='|', index=False)

# filtered_df = oncoFrame[oncoFrame['AF'] >= 0.8]
filtered_df = oncoFrame[(oncoFrame['Sample'] == 'Mixed') & (oncoFrame['AF'] >= 0.8)]
filtered_df1 = oncoFrame[(oncoFrame['Sample'] == 'Mixed') & (oncoFrame['AF'] <= 0.65) & (oncoFrame['AF'] >= 0.35)]
# filtered_df1 = oncoFrame[oncoFrame['AF'] <= 0.65]
# filtered_df1 = filtered_df1[filtered_df1['AF'] >= 0.35]
filtered_df = oncoFrame[oncoFrame['Name'].isin(filtered_df['Name'])]
filtered_df1 = oncoFrame[oncoFrame['Name'].isin(filtered_df1['Name'])]

file_path = r'C:\Users\Eugene Cho\Desktop\BioinformaticsFiles\CGI_analysis\alterations.tsv'
altFrame = pd.read_csv(file_path, delimiter='\t')
cgiMatrix = pd.DataFrame({
    'Chromosome': altFrame['chr'],
    'StartPos': altFrame['pos'],
    'Ref': altFrame['ref'],
    'Alt': altFrame['alt'],
    'Oncogenic Prediction': altFrame['CGI-Oncogenic Prediction'],
    'ID': altFrame['Input ID'],
    'Gene': altFrame['CGI-Gene'],
    'HGVSc': altFrame['CGI-HGVSc'],
    'HGVSp': altFrame['CGI-HGVSp']
})

cgiMatrix['Name'] = cgiMatrix.apply(create_name, axis=1)
conditions = (
    (cgiMatrix['Oncogenic Prediction'] == 'driver (boostDM: tissue specific model)') |
    (cgiMatrix['Oncogenic Prediction'] == 'driver (oncodriveMUT)')
)
cgiMatrix = cgiMatrix[conditions]
# print(cgiMatrix)
filtered_df2 = oncoFrame[oncoFrame['Name'].isin(cgiMatrix['Name'])]
df2_subset = pd.DataFrame({
    'Name': cgiMatrix['Name'],
    'Gene': cgiMatrix['Gene']
})
filtered_df2 = filtered_df2.set_index('Name')
df2_subset = df2_subset.set_index('Name')
filtered_df2 = filtered_df2.join(df2_subset)
filtered_df2.to_csv('output3.txt', sep='\t', index=True)

print(filtered_df2)
metelMatrix = pd.DataFrame({
    'Patient ID': cgiMatrix['Name'],
    'Gene': cgiMatrix['Gene'],
    'HGVSc': cgiMatrix['HGVSc'],
    'HGVSp': cgiMatrix['HGVSp'],
})
metelMatrix = metelMatrix.set_index('Patient ID')
metelMatrix1 = filtered_df2
metelMatrix1['Gene'] = metelMatrix['Gene']
metelMatrix1['HGVSc'] = metelMatrix['HGVSc']
metelMatrix1['HGVSp'] = metelMatrix['HGVSp']

filtered_samp1 = metelMatrix1[metelMatrix1['Sample'] == 'C7'].drop(columns='Name1')
filtered_samp2 = metelMatrix1[metelMatrix1['Sample'] == 'C9'].drop(columns='Name1')
# print(filtered_samp1)
# print(filtered_samp2)
# filtered_samp1 = filtered_samp1.merge(filtered_samp2, how='left', copy=False)
filtered_samp1.to_csv('outputCGI1.txt', sep='\t', index=True)

mergedMatrix = pd.DataFrame({
    'Gene': filtered_samp1['Gene'],
    'HGVSc': filtered_samp1['HGVSc'],
    'HGVSp': filtered_samp1['HGVSp'],
})
mergedMatrix['VAF1'] = filtered_samp1['AF']
mergedMatrix['VAF2'] = filtered_samp2['AF']
mask = mergedMatrix.isnull().any(axis=1)
mergedMatrix = mergedMatrix[~mask]
# print(mergedMatrix)

mergedMatrix.to_csv('outputCGI.txt', sep='\t', index=True)


# print(filtered_df['Name'])
# print(filtered_df1['Name'])
# print(filtered_df2['Name'])

oncoFrame = filtered_df
print(oncoFrame)
oncoFrame = oncoFrame.reset_index()
duplicates = oncoFrame[oncoFrame.duplicated(subset=['Name', 'Sample'], keep=False)]
print(duplicates)
oncoFrame = oncoFrame.drop_duplicates(subset=['Name', 'Sample'])
secondFrame = oncoFrame.pivot(index='Name', columns='Sample', values='AF')
nan_counts = secondFrame.isna().sum(axis=1)
rows_with_four_or_more_nans = (nan_counts >= 4).sum()
rows_with_three_or_more_nans = (nan_counts >= 3).sum()
rows_with_two_or_more_nans = (nan_counts >= 2).sum()

print(f'Number of rows with four or more NaNs: {rows_with_four_or_more_nans}')
print(f'Number of rows with three or more NaNs: {rows_with_three_or_more_nans}')
print(f'Number of rows with two or more NaNs: {rows_with_two_or_more_nans}')

secondFrame_cleaned = secondFrame.dropna()
mask = secondFrame.isnull()

print(secondFrame)
print(secondFrame_cleaned)


# clustergram = dashbio.Clustergram(
#         data=secondFrame_cleaned.values,
#         row_labels=secondFrame_cleaned.index.tolist(),
#         column_labels=secondFrame_cleaned.columns.tolist(),
#         center_values=False,
#         height=900,
#         width=1100,
#         hidden_labels='row',
#         color_map= 'viridis'
#     )
# app = dash.Dash(__name__)
# app.layout = html.Div([
#     dcc.Graph(figure=clustergram)  # Use dcc.Graph to display the clustergram
# ])

# if __name__ == '__main__':
#     app.run_server(debug=True)
# print(oncoFrame)
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
    yaxis_showticklabels=False
)
fig.update_yaxes(tickfont=dict(size=9))  # Set the font size here
fig.write_html('oncoTest.html') 
fig.show()