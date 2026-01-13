import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


CSV_FILES = [
    "/gpfs/scratch/bsc18/bsc665560/Docs/hsapiens/GRCh38_noN.fna_map-ont_minimizers_key_value_sorted_index_table.csv",
    "/gpfs/scratch/bsc18/bsc665560/Docs/drosophila/drosophila.fasta_map-ont_minimizers_key_value_sorted_keys_ordered.csv",
    "/gpfs/scratch/bsc18/bsc665560/Docs/saccharomyces/saccharomyces.fasta_map-ont_minimizers_key_value_sorted_keys_ordered.csv",
]

LABELS = [
    "GRCh38",
    "Drosophila",
    "Saccharomyces",
]

OUTPUT_FILE = "CDF_comparison.png"
MAX_POINTS = 200_000

def load_cdf_from_csv(file_path, max_points):
    df = pd.read_csv(file_path)

    if "key" not in df.columns:
        raise ValueError(f"{file_path} no tiene columna 'key'")

    keys = np.sort(df["key"].to_numpy())
    n = len(keys)

    # Subsampling
    if n > max_points:
        step = n // max_points
        keys = keys[::step]
        n = len(keys)

    # CDF empírica normalizada
    cdf_y = np.arange(1, n + 1) / n
    return keys, cdf_y


def plot_multiple_cdfs():
    plt.figure(figsize=(16, 7), dpi=200)

    colors = ["#5b8def", "#e76f51", "#2a9d8f", "#9b5de5"]

    for i, (file, label) in enumerate(zip(CSV_FILES, LABELS)):
        x, y = load_cdf_from_csv(file, MAX_POINTS)
        plt.plot(
            x, y,
            linewidth=1.5,
            color=colors[i % len(colors)],
            label=label
        )


    plt.title("CDF comparison (normalized)", fontsize=16)
    plt.xlabel("Key", fontsize=16)
    plt.ylabel("Cumulative Probability", fontsize=16)

    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.legend(fontsize=16)

    plt.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig(OUTPUT_FILE, dpi=200)
    plt.close()

    print(f"CDFs guardadas en {OUTPUT_FILE}")


if __name__ == "__main__":
    plot_multiple_cdfs()
