import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind


# Read the data from the CSV file
data = pd.read_csv('/Users/drake/Desktop/tf_trial/analysis/analysis2/untitled_folder_2/error_snp_mismatch_results.csv')

# Box Plots for Error Rates
plt.figure(figsize=(8, 6))
sns.boxplot(data=data[['MinION1_ErrorRate', 'MinION2_ErrorRate']])
plt.ylabel('Error Rate')
plt.title('Error Rates for MinION1 and MinION2')
plt.xticks(ticks=[0, 1], labels=['MinION1', 'MinION2'])
plt.show()

# Bar Charts for Total Indels and SNPs
plt.figure(figsize=(10, 6))
sns.barplot(x='Sample_ID', y='MinION1_Indels', data=data, color='red', label='MinION1 Indels')
sns.barplot(x='Sample_ID', y='MinION2_Indels', data=data, color='blue', label='MinION2 Indels', alpha=0.7)
sns.barplot(x='Sample_ID', y='MinION1_SNPs', data=data, color='orange', label='MinION1 SNPs')
sns.barplot(x='Sample_ID', y='MinION2_SNPs', data=data, color='green', label='MinION2 SNPs', alpha=0.7)
plt.xlabel('Sample ID')
plt.ylabel('Counts')
plt.title('Total Indels and SNPs for MinION1 and MinION2')
plt.legend()
plt.xticks(rotation=90)
plt.show()

# Stacked Bar Charts for Error Types
data_melted = data.melt(id_vars='Sample_ID', value_vars=['MinION1_Indels', 'MinION2_Indels', 'MinION1_SNPs', 'MinION2_SNPs'])
plt.figure(figsize=(10, 6))
sns.barplot(x='Sample_ID', y='value', hue='variable', data=data_melted)
plt.xlabel('Sample ID')
plt.ylabel('Counts')
plt.title('Error Types for MinION1 and MinION2')
plt.legend(title='Error Types')
plt.xticks(rotation=90)
plt.show()

# Line Plots for Error Rates
plt.figure(figsize=(10, 6))
sns.lineplot(data=data[['MinION1_ErrorRate', 'MinION2_ErrorRate']], marker='o')
plt.xlabel('Sample ID')
plt.ylabel('Error Rate')
plt.title('Error Rate Trends for MinION1 and MinION2')
plt.legend(labels=['MinION1', 'MinION2'])
plt.xticks(range(len(data['Sample_ID'])), data['Sample_ID'], rotation=90)
plt.show()

# Heatmap for Pairwise Identities
plt.figure(figsize=(10, 6))
sns.heatmap(data[['MinION1_PairwiseIdentity', 'MinION2_PairwiseIdentity']], annot=True, cmap='coolwarm')
plt.xlabel('Pipeline')
plt.ylabel('Sample ID')
plt.title('Pairwise Identities for MinION1 and MinION2')
plt.xticks(ticks=[0.5, 1.5], labels=['MinION1', 'MinION2'])
plt.show()

# Box Plots for Error Rates
plt.figure(figsize=(10, 6))
sns.boxplot(data=data[['MinION1_ErrorRate', 'MinION2_ErrorRate']])
plt.ylabel('Error Rate')
plt.title('Error Rates for MinION1 and MinION2')
plt.xticks(ticks=[0, 1], labels=['MinION1', 'MinION2'])
plt.show()

# Violin Plots for Error Rates
plt.figure(figsize=(10, 6))
sns.violinplot(data=data[['MinION1_ErrorRate', 'MinION2_ErrorRate']])
plt.ylabel('Error Rate')
plt.title('Error Rate Distribution for MinION1 and MinION2')
plt.xticks(ticks=[0, 1], labels=['MinION1', 'MinION2'])
plt.show()

# Violin Plots for Indels and SNPs Error Profiles for MinION1 and MinION2
plt.figure(figsize=(12, 6))

plt.subplot(1, 2, 1)
sns.violinplot(data=data[['MinION1_Indels', 'MinION1_SNPs']], inner='quart', palette='Set1')
plt.ylabel('Counts')
plt.title('MinION1 Indels and SNPs Error Profiles')
plt.xticks(ticks=[0, 1], labels=['Indels', 'SNPs'])

plt.subplot(1, 2, 2)
sns.violinplot(data=data[['MinION2_Indels', 'MinION2_SNPs']], inner='quart', palette='Set2')
plt.ylabel('Counts')
plt.title('MinION2 Indels and SNPs Error Profiles')
plt.xticks(ticks=[0, 1], labels=['Indels', 'SNPs'])

plt.tight_layout()
plt.show()
# Melt the data for easier visualization
data_melted = data.melt(id_vars='Sample_ID', value_vars=['MinION1_Indels', 'MinION1_SNPs', 'MinION2_Indels', 'MinION2_SNPs'],
                        var_name='Error Type', value_name='Counts')

# Extract the pipeline label (MinION1 or MinION2) from the 'Error Type' column
data_melted['Pipeline'] = data_melted['Error Type'].str.extract(r'(MinION\d)')

# Create a single violin plot for MinION1 and MinION2
plt.figure(figsize=(10, 6))
sns.violinplot(x='Error Type', y='Counts', hue='Pipeline', data=data_melted, palette='pastel', split=True, inner='quart')
plt.xlabel('Error Type')
plt.ylabel('Counts')
plt.title('Error Profiles (Indels and SNPs) for MinION1 and MinION2')
plt.legend(title='Pipeline', labels=['MinION1', 'MinION2'])
plt.show()


# Melt the data for easier visualization
data_melted = data.melt(id_vars='Sample_ID', value_vars=['MinION1_Indels', 'MinION1_SNPs', 'MinION2_Indels', 'MinION2_SNPs'],
                        var_name='Error Type', value_name='Counts')

# Extract the pipeline label (MinION1 or MinION2) and the error type (Indels or SNPs) from the 'Error Type' column
data_melted['Pipeline'] = data_melted['Error Type'].str.extract(r'(MinION\d)')
data_melted['Error Type'] = data_melted['Error Type'].str.extract(r'_(\w+)$')

# Create a single violin plot for MinION1 and MinION2, showing SNPs on one side and indels on the other side
plt.figure(figsize=(12, 6))
sns.violinplot(x='Error Type', y='Counts', hue='Pipeline', data=data_melted, palette='pastel', split=True, inner='quart')
plt.xlabel('Error Type')
plt.ylabel('Counts')
plt.title('Error Profiles (Indels and SNPs) for MinION1 and MinION2')
plt.legend(title='Pipeline', labels=['MinION1', 'MinION2'])
plt.show()


# Melt the data for easier visualization
data_melted = data.melt(id_vars='Sample_ID', value_vars=['MinION1_Indels', 'MinION1_SNPs', 'MinION2_Indels', 'MinION2_SNPs'],
                        var_name='Error Type', value_name='Counts')

# Extract the pipeline label (MinION1 or MinION2) and the error type (Indels or SNPs) from the 'Error Type' column
data_melted['Pipeline'] = data_melted['Error Type'].str.extract(r'(MinION\d)')
data_melted['Error Type'] = data_melted['Error Type'].str.extract(r'_(\w+)$')

# Create a single violin plot for MinION1 and MinION2, showing SNPs on one side and indels on the other side
plt.figure(figsize=(12, 6))
sns.violinplot(x='Error Type', y='Counts', hue='Pipeline', data=data_melted, palette='pastel', split=True, inner=None)
plt.xlabel('Error Type')
plt.ylabel('Counts')
plt.title('Error Profiles (Indels and SNPs) for MinION1 and MinION2')
plt.legend(title='Pipeline', labels=['MinION1', 'MinION2'])
plt.show()

# Melt the data for easier visualization
data_melted = data.melt(id_vars='Sample_ID', value_vars=['MinION1_Indels', 'MinION1_SNPs', 'MinION2_Indels', 'MinION2_SNPs'],
                        var_name='Error Type', value_name='Counts')

# Extract the pipeline label (MinION1 or MinION2) and the error type (Indels or SNPs) from the 'Error Type' column
data_melted['Pipeline'] = data_melted['Error Type'].str.extract(r'(MinION\d)')
data_melted['Error Type'] = data_melted['Error Type'].str.extract(r'_(\w+)$')

# Create a figure with two subplots (1 row, 2 columns)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 6))


# Melt the data for easier visualization
data_melted = data.melt(id_vars='Sample_ID', value_vars=['MinION1_Indels', 'MinION1_SNPs', 'MinION2_Indels', 'MinION2_SNPs'],
                        var_name='Error Type', value_name='Counts')

# Extract the pipeline label (MinION1 or MinION2) and the error type (Indels or SNPs) from the 'Error Type' column
data_melted['Pipeline'] = data_melted['Error Type'].str.extract(r'(MinION\d)')
data_melted['Error Type'] = data_melted['Error Type'].str.extract(r'_(\w+)$')

# Create a figure with two subplots (1 row, 2 columns)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 6))

# Melt the data for easier visualization
data_melted = data.melt(id_vars='Sample_ID', value_vars=['MinION1_Indels', 'MinION1_SNPs', 'MinION2_Indels', 'MinION2_SNPs'],
                        var_name='Error Type', value_name='Counts')

# Extract the pipeline label (MinION1 or MinION2) and the error type (Indels or SNPs) from the 'Error Type' column
data_melted['Pipeline'] = data_melted['Error Type'].str.extract(r'(MinION\d)')
data_melted['Error Type'] = data_melted['Error Type'].str.extract(r'_(\w+)$')

# Create a single violin plot for MinION1 and MinION2, showing SNPs on one side and indels on the other side
plt.figure(figsize=(12, 6))
sns.violinplot(x='Error Type', y='Counts', hue='Pipeline', data=data_melted, palette='pastel', split=True, inner='stick')
plt.xlabel('Error Type')
plt.ylabel('Counts')
plt.title('Error Profiles (Indels and SNPs) for MinION1 and MinION2')
plt.legend(title='Pipeline', labels=['MinION1', 'MinION2'])
plt.show()


# Melt the data for easier visualization
data_melted = data.melt(id_vars='Sample_ID', value_vars=['MinION1_Indels', 'MinION1_SNPs', 'MinION2_Indels', 'MinION2_SNPs'],
                        var_name='Error Type', value_name='Counts')

# Extract the pipeline label (MinION1 or MinION2) and the error type (Indels or SNPs) from the 'Error Type' column
data_melted['Pipeline'] = data_melted['Error Type'].str.extract(r'(MinION\d)')
data_melted['Error Type'] = data_melted['Error Type'].str.extract(r'_(\w+)$')

# Create a single violin plot for MinION1 and MinION2, showing SNPs on one side and indels on the other side
plt.figure(figsize=(12, 6))
sns.violinplot(x='Error Type', y='Counts', hue='Pipeline', data=data_melted, palette='pastel', split=True, inner='stick')
plt.xlabel('Error Type')
plt.ylabel('Counts')
plt.title('Error Profiles (Indels and SNPs) for MinION1 and MinION2')

# Customize the legend
handles, labels = plt.gca().get_legend_handles_labels()
unique_labels = sorted(set(labels), key=labels.index)
plt.legend(handles=handles[:len(unique_labels)], labels=unique_labels, title='Pipeline')

plt.show()

# Melt the data for easier visualization
data_melted = data.melt(id_vars='Sample_ID', value_vars=['MinION1_Indels', 'MinION1_SNPs', 'MinION2_Indels', 'MinION2_SNPs'],
                        var_name='Error Type', value_name='Counts')

# Extract the pipeline label (MinION1 or MinION2) and the error type (Indels or SNPs) from the 'Error Type' column
data_melted['Pipeline'] = data_melted['Error Type'].str.extract(r'(MinION\d)')
data_melted['Error Type'] = data_melted['Error Type'].str.extract(r'_(\w+)$')

# Create a single violin plot for MinION1 and MinION2, showing SNPs on one side and indels on the other side
plt.figure(figsize=(12, 6))
sns.violinplot(x='Error Type', y='Counts', hue='Pipeline', data=data_melted, palette='pastel', split=True, inner='stick')

# Customize the labels for the KDE plots inside the violin
for i, violin in enumerate(ax.collections):
    if i % 2 == 0:
        violin.set_label('KDE MinION1')
    else:
        violin.set_label('KDE MinION2')

plt.xlabel('Error Type')
plt.ylabel('Counts')
plt.title('Error Profiles (Indels and SNPs) for MinION1 and MinION2')

# Customize the legend
handles, labels = plt.gca().get_legend_handles_labels()
unique_labels = sorted(set(labels), key=labels.index)
plt.legend(handles=handles[:len(unique_labels)], labels=unique_labels, title='Pipeline')

plt.show()
