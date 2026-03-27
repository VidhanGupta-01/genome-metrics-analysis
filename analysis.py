import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

# Load data
df = pd.read_csv("data/sample_genome_metrics.csv")

# Select numeric features
X = df[["genes", "transcripts", "exons", "completeness"]]

# PCA
pca = PCA(n_components=2)
components = pca.fit_transform(X)

# Plot
plt.figure()
plt.scatter(components[:, 0], components[:, 1])

for i, genome in enumerate(df["genome"]):
    plt.text(components[i, 0], components[i, 1], genome)

plt.title("Genome PCA Analysis")
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.savefig("pca_plot.png")
