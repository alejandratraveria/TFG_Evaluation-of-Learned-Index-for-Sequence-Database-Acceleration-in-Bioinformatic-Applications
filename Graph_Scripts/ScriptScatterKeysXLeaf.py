import pandas as pd
import matplotlib.pyplot as plt
import os

baseline_file = '/gpfs/scratch/bsc18/bsc665560/Docs/RMI_Baseline_LISA_keys_per_leaf.csv'

OUTPUT_DIR = 'graphs'
OUTPUT_FILE = os.path.join(OUTPUT_DIR, 'Baseline_LISA_keys_per_leaf.png')

MAX_KEYS_TO_PLOT = 10_000

print("Cargando datos...")

try:
    df_baseline = pd.read_csv(baseline_file)


except FileNotFoundError as e:
    print(f"Error: No se encontró un archivo.\n{e}")
    exit(1)


df_baseline = df_baseline[df_baseline['keys'] <= MAX_KEYS_TO_PLOT]


os.makedirs(OUTPUT_DIR, exist_ok=True)

plt.figure(figsize=(18, 9), dpi=600)

s_size = 10
alpha_val = 0.6


plt.scatter(
    df_baseline['leaf_id'],
    df_baseline['keys'],
    s=s_size,
    alpha=alpha_val,
    color='#1f77b4',
    label='Baseline LISA',
    edgecolor='none'
)

plt.title(
    'Baseline-LISA: Key Distribution per Leaf',
    fontsize=14
)
plt.xlabel('Leaf id', fontsize=14)
plt.ylabel('Number of assigned keys', fontsize=14)

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

plt.legend(fontsize=14, framealpha=0.9, loc='upper right')
plt.grid(True, linestyle='--', alpha=0.35)

plt.ylim(0, 250)
plt.tight_layout()

plt.savefig(OUTPUT_FILE, dpi=600)
plt.close()

print(f"Scatter plot (solo Baseline LISA) guardado en: {OUTPUT_FILE}")
