import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, ScalarFormatter
import os
import numpy as np

baseline_df = pd.read_csv(
    '/gpfs/scratch/bsc18/bsc665560/Docs/serials/LISA/LISA_param_error_GRCh38.csv'
)

seg16_df = pd.read_csv(
    '/gpfs/scratch/bsc18/bsc665560/Docs/serials/LISA/LISA_16SEQ_param_error_GRCh38.csv'
)

baseline_df = baseline_df[baseline_df['err'] < 1_000_000]
seg16_df    = seg16_df[seg16_df['err'] < 1_000_000]

os.makedirs('graphs', exist_ok=True)

plt.figure(figsize=(14, 14), dpi=400)


plt.scatter(
    baseline_df['model'],
    baseline_df['err'],
    color='#1f77b4',
    s=10,
    alpha=0.5,
    edgecolors='none',
    label='Baseline-LISA'
)


plt.scatter(
    seg16_df['model'],
    seg16_df['err'],
    color='#fd633d',
    s=10,
    alpha=0.5,
    edgecolors='none',
    label='16 Segments-LISA'
)

plt.title(
    'Baseline-LISA vs 16 Segments-LISA: Distribution of error value per leaf (GRCh38)',
    fontsize=16
)
plt.xlabel('Leaf index (model)', fontsize=16)
plt.ylabel('Error (err)', fontsize=16)

ax = plt.gca()

ax.yaxis.set_major_locator(MultipleLocator(5))
ax.yaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
ax.ticklabel_format(style='plain', axis='y', useOffset=False)

ax.xaxis.set_major_formatter(ScalarFormatter(useOffset=False))
ax.ticklabel_format(style='plain', axis='x', useOffset=False)

ax.tick_params(axis='both', which='major', labelsize=16)
ax.tick_params(axis='both', which='minor', labelsize=16)

plt.xticks(rotation=45, ha='right')

plt.grid(True, linestyle='--', alpha=0.5)

plt.ylim(0, 135)

plt.legend(fontsize=14)
plt.tight_layout()

output = '/gpfs/scratch/bsc18/bsc665560/Docs/graphs/16SEQ_error_scatter.png'
plt.savefig(output, dpi=400)
plt.close()

print(f"Gráfico guardado en: {output}")
