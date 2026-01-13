import matplotlib.pyplot as plt
from collections import Counter

def load_error_cdf_from_csv(path):
    errors = []

    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("model"):
                continue

            parts = line.split(",")
            if len(parts) < 2:
                continue

            try:
                err = int(parts[1].strip())
            except ValueError:
                continue

            errors.append(err)

    if not errors:
        print(f"[WARNING] File empty or invalid: {path}")
        return [], [], []

    counter = Counter(errors)
    xs = sorted(counter.keys())

    cumulative = []
    running = 0
    total = len(errors)

    for e in xs:
        running += counter[e]
        cumulative.append(running / total)

    return xs, cumulative, errors


baseline_file = '/gpfs/scratch/bsc18/bsc665560/Docs/serials/LISA/LISA_param_error_GRCh38.csv'
file_seg2     = '/gpfs/scratch/bsc18/bsc665560/Docs/serials/LISA/LISA_SEQ_param_error_GRCh38.csv'
file_seg4     = '/gpfs/scratch/bsc18/bsc665560/Docs/serials/LISA/LISA_4SEQ_param_error_GRCh38.csv'
file_seg8     = '/gpfs/scratch/bsc18/bsc665560/Docs/serials/LISA/LISA_8SEQ_param_error_GRCh38.csv'
file_seg16    = '/gpfs/scratch/bsc18/bsc665560/Docs/serials/LISA/LISA_16SEQ_param_error_GRCh38.csv'

files_config = [
    (baseline_file, "Baseline-LISA",                "#1f77b4"),
    (file_seg2,     "Segmented-LISA (2 segments)",  "#ffbf00"),
    (file_seg4,     "Segmented-LISA (4 segments)",  "#2ca02c"),
    (file_seg8,     "Segmented-LISA (8 segments)",  "#d62728"),
    (file_seg16,    "Segmented-LISA (16 segments)", "#9467bd"),
]

plt.figure(figsize=(16, 8), dpi=220)

for path, label, color in files_config:
    xs, ys, _ = load_error_cdf_from_csv(path)
    if xs:
        plt.plot(xs, ys, linewidth=2.5, color=color, label=label)

plt.title(
    "Baseline-LISA vs Segmented-LISA: CDF of Leaf Prediction Error (GRCh38)",
    fontsize=16
)
plt.xlabel("Leaf Prediction Error (absolute)", fontsize=16)
plt.ylabel("Cumulative Distribution Function (normalized)", fontsize=16)

plt.xscale("log")
plt.ylim(0, 1.03)

plt.grid(True, linestyle='--', alpha=0.5)
plt.legend(fontsize=16)

plt.tight_layout()
plt.savefig("Baseline_vs_Segmented_error_logscale.png", dpi=300)
plt.close()

print("Generated figure: Baseline_vs_Segmented_error_logscale.png")
