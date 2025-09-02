import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

def read_clusters(file_path):
    clusters = []
    with open(file_path, 'rb') as file:  # Open the file in binary mode
        for line in file:
            try:
                decoded_line = line.decode('utf-8').strip()
                cluster_id, read_id = decoded_line.split()
                clusters.append([cluster_id, read_id])
            except UnicodeDecodeError as e:
                print(f"Skipping malformed line: {line.strip()} ({e})")
            except ValueError as e:
                print(f"Skipping malformed line: {decoded_line} ({e})")
    return pd.DataFrame(clusters, columns=['ClusterID', 'ReadID'])

def plot_clusters(df):
    plt.figure(figsize=(10, 6))
    plt.scatter(df['ClusterID'], df['ReadID'], alpha=0.5)
    plt.xlabel('Cluster ID')
    plt.ylabel('Read ID')
    plt.title('Scatter Plot of Clusters')
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python visualize_clusters.py <clusters.out>")
        sys.exit(1)

    file_path = sys.argv[1]
    df = read_clusters(file_path)
    plot_clusters(df)
