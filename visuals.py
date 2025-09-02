import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Read the CSV file into a pandas DataFrame
data = pd.read_csv('/Users/drake/Desktop/tf_trial/analysis/analysis2/results.csv')

# Remove any rows with missing values
data = data.dropna()

# Box Plot
sns.boxplot(data=data[['MinION 1 SNPs/Mismatches', 'MinION 1 Indels', 'MinION 2 SNPs/Mismatches', 'MinION 2 Indels']])
plt.title('Box Plot of SNPs/Mismatches and Indels')
plt.xlabel('Technology')
plt.ylabel('Count')
plt.show()

# Grouped Bar Chart
sns.barplot(data=data[['MinION 1 SNPs/Mismatches', 'MinION 1 Indels', 'MinION 2 SNPs/Mismatches', 'MinION 2 Indels']], ci=None)
plt.title('Grouped Bar Chart of SNPs/Mismatches and Indels')
plt.xlabel('Sample ID')
plt.ylabel('Count')
plt.legend(['MinION 1 SNPs/Mismatches', 'MinION 1 Indels', 'MinION 2 SNPs/Mismatches', 'MinION 2 Indels'])
plt.show()

# Line Plot with Error Bars
sns.lineplot(data=data[['MinION 1 Error Rate', 'MinION 2 Error Rate']])
plt.title('Line Plot of Error Rates')
plt.xlabel('Sample ID')
plt.ylabel('Error Rate')
plt.legend(['MinION 1', 'MinION 2'])
plt.show()

# Stacked Area Chart
stacked_data = data[['MinION 1 SNPs/Mismatches', 'MinION 1 Indels', 'MinION 2 SNPs/Mismatches', 'MinION 2 Indels']]
stacked_data['Sample ID'] = data['Sample ID']
stacked_data = stacked_data.set_index('Sample ID')
stacked_data.plot(kind='area', stacked=True)
plt.title('Stacked Area Chart of SNPs/Mismatches and Indels')
plt.xlabel('Sample ID')
plt.ylabel('Count')
plt.show()

# Scatter Plot Matrix
sns.pairplot(data=data[['MinION 1 SNPs/Mismatches', 'MinION 1 Indels', 'MinION 2 SNPs/Mismatches', 'MinION 2 Indels']])
plt.show()

import pandas as pd
import matplotlib.pyplot as plt

# Calculate the total SNPs and Indels for MinION 1 and MinION 2
data['MinION 1 Total'] = data['MinION 1 SNPs/Mismatches'] + data['MinION 1 Indels']
data['MinION 2 Total'] = data['MinION 2 SNPs/Mismatches'] + data['MinION 2 Indels']

# Set the sample IDs as the x-axis labels
x_labels = data['Sample ID']

# Set the width of the bars
bar_width = 0.35

# Set the positions of the bars on the x-axis
bar_positions_1 = range(len(data))
bar_positions_2 = [x + bar_width for x in bar_positions_1]

# Create the grouped bar chart
plt.bar(bar_positions_1, data['MinION 1 Total'], width=bar_width, label='MinION 1')
plt.bar(bar_positions_2, data['MinION 2 Total'], width=bar_width, label='MinION 2')

# Set the labels and title
plt.xlabel('Sample ID')
plt.ylabel('Count')
plt.title('Total SNPs and Indels for MinION 1 and MinION 2')
plt.xticks([r + bar_width/2 for r in range(len(data))], x_labels, rotation='vertical')

# Add a legend
plt.legend()

# Show the chart
plt.tight_layout()
plt.show()


# Calculate the overall totals of SNPs and Indels for MinION 1 and MinION 2
minion1_total = data['MinION 1 SNPs/Mismatches'].sum() + data['MinION 1 Indels'].sum()
minion2_total = data['MinION 2 SNPs/Mismatches'].sum() + data['MinION 2 Indels'].sum()

# Create the grouped bar chart
plt.bar(['MinION 1', 'MinION 2'], [minion1_total, minion2_total])

# Set the labels and title
plt.xlabel('Technology')
plt.ylabel('Total Count')
plt.title('Overall Total SNPs and Indels for MinION 1 and MinION 2')

# Show the chart
plt.show()