import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import os
import numpy as np

BASELINE_CSV = '/gpfs/scratch/bsc18/bsc665560/Docs/serials/LISA/LISA_param_error_GRCh38.csv'
SEG16_CSV    = '/gpfs/scratch/bsc18/bsc665560/Docs/serials/LISA/LISA_16SEQ_param_error_GRCh38.csv'

OUTPUT_DIR = 'graphs'
OUTPUT_FILE = os.path.join(OUTPUT_DIR, '16SEQ_error_count.png')

MAX_ERR = 170
REMOVE_OUTLIER = True
OUTLIER_THRESHOLD = 1_000_000

YMAX = 5_200_000  

def load_error_counts(csv_path):
    df = pd.read_csv(csv_path)

    if not {'model', 'err'}.issubset(df.columns):
        raise ValueError("El CSV debe tener columnas 'model' y 'err'")

    if REMOVE_OUTLIER:
        df = df[df['err'] < OUTLIER_THRESHOLD]

    counts = (
        df['err']
        .value_counts()
        .sort_index()
    )

    counts = counts[counts.index <= MAX_ERR]
    return counts


def plot_error_bars(ax, counts, title, bar_color, ymax=None):
    x_labels = counts.index.astype(str)
    x_pos = np.arange(len(x_labels))

    bars = ax.bar(
        x_pos,
        counts.values,
        color=bar_color,
        edgecolor='black',
        linewidth=0.3,
        width=0.8
    )

    ax.set_title(title, fontsize=16)
    ax.set_xlabel('Error value', fontsize=16)
    ax.set_ylabel('Number of models', fontsize=16)

    ax.set_xticks(x_pos)
    ax.set_xticklabels(x_labels, fontsize=14, rotation=45, ha='right')
    ax.tick_params(axis='y', labelsize=14)

    ax.grid(axis='y', linestyle='--', alpha=0.5)

    ax.yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    ax.ticklabel_format(style='plain', axis='y', useOffset=False)

    if ymax is not None:
        ax.set_ylim(0, ymax)

    y_offset = (ax.get_ylim()[1] if ymax is not None else max(counts.values)) * 0.01
    for bar, err_value in zip(bars, counts.index):
        height = bar.get_height()
        if err_value >= 10 and height > 0:
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                height + y_offset,
                f'{int(height)}',
                ha='center',
                va='bottom',
                fontsize=11,
                rotation=90
            )

baseline_counts = load_error_counts(BASELINE_CSV)
seg16_counts    = load_error_counts(SEG16_CSV)

os.makedirs(OUTPUT_DIR, exist_ok=True)

fig, (ax1, ax2) = plt.subplots(
    2, 1,
    figsize=(16, 12),
    dpi=400
)

plot_error_bars(
    ax1,
    baseline_counts,
    'Baseline-LISA: Number of models per error value (GRCh38)',
    bar_color='#1f77b4',
    ymax=YMAX
)

plot_error_bars(
    ax2,
    seg16_counts,
    '16 Segments-LISA: Number of models per error value (GRCh38)',
    bar_color="#fd633d",
    ymax=YMAX
)

plt.tight_layout()
plt.savefig(OUTPUT_FILE, dpi=400)
plt.close()

print(f"Gráfico guardado en: {OUTPUT_FILE}")
