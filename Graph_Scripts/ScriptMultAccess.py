import re
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import os


# Baseline
file_baseline = '/gpfs/scratch/bsc18/bsc665560/Docs/serials/LISA/LISA_mm_metrics_GRCh38.csv'

# Segmented LISA
file_seg2  = '/gpfs/scratch/bsc18/bsc665560/Docs/serials/LISA/LISA_SEQ_mm_metrics_GRCh38.csv'
file_seg4  = '/gpfs/scratch/bsc18/bsc665560/Docs/serials/LISA/LISA_4SEQ_mm_metrics_GRCh38.csv'
file_seg8  = '/gpfs/scratch/bsc18/bsc665560/Docs/serials/LISA/LISA_8SEQ_mm_metrics_GRCh38.csv'
file_seg16 = '/gpfs/scratch/bsc18/bsc665560/Docs/serials/LISA/LISA_16SEQ_mm_metrics_GRCh38.csv'


files_config = [
    (file_baseline, "Baseline-LISA",        '#1f77b4', 'o'),  # Azul
    (file_seg2,     "2 Segments-LISA",      '#ffbf00', 's'),  # Amarillo
    (file_seg4,     "4 Segments-LISA",      '#2ca02c', '^'),  # Verde
    (file_seg8,     "8 Segments-LISA",      '#d62728', 'v'),  # Rojo
    (file_seg16,    "16 Segments-LISA",     '#9467bd', 'D'),  # Morado
]


def load_accesses_from_file(path):
    counts = {}

    pattern_exact = re.compile(r"Total number of accesses without range adjustments,(\d+)")
    pattern_adjust = re.compile(r"Total number of accesses with (\d+) range adjustment[s]?,(\d+)")
    pattern_more = re.compile(r"Total number of accesses with more than 8 range adjustments,(\d+)")
    pattern_total_accesses = re.compile(r"Total number of accesses to the index,(\d+)")

    total_accesses_real = None
    total_queries = 0

    try:
        with open(path, "r") as f:
            for line in f:
                line = line.strip()

                m = pattern_total_accesses.match(line)
                if m:
                    total_accesses_real = int(m.group(1))
                    continue

                m = pattern_exact.match(line)
                if m:
                    cnt = int(m.group(1))
                    counts[0] = cnt
                    total_queries += cnt
                    continue

                m = pattern_adjust.match(line)
                if m:
                    adj = int(m.group(1))
                    cnt = int(m.group(2))
                    counts[adj] = cnt
                    total_queries += cnt
                    continue

                m = pattern_more.match(line)
                if m:
                    cnt = int(m.group(1))
                    counts[9] = cnt   # bucket >8
                    total_queries += cnt
                    continue

    except FileNotFoundError:
        print(f"Error: Archivo no encontrado -> {path}")
        return [], [], 0

    if total_queries == 0:
        print(f"Advertencia: No se encontraron datos de accesos en {path}")
        return [], [], 0

    xs = sorted(counts.keys())

    # Accesos estimados
    accesses_estimated = []
    for k in xs:
        if k == 9:
            accesos = counts[k] * 10
        else:
            accesos = counts[k] * (k + 1)
        accesses_estimated.append(accesos)

    A_est = sum(accesses_estimated)
    A_real = total_accesses_real if total_accesses_real else A_est
    correction_factor = A_real / A_est if A_est > 0 else 1.0

    accesses_corrected = [a * correction_factor for a in accesses_estimated]

    # Acumulado
    cumulative = []
    running = 0
    for a in accesses_corrected:
        running += a
        cumulative.append(running)

    avg_access = A_real / total_queries if total_queries > 0 else 0.0

    print(f"\n[{path}]")
    print(f"  Accesos reales reportados  : {A_real:,}")
    print(f"  Queries totales            : {total_queries:,}")
    print(f"  Average accesses/query     : {avg_access:.6f}")

    return xs, cumulative, avg_access


plt.figure(figsize=(16, 9), dpi=300)

xs_all_union = set()

for path, label, color, marker in files_config:
    xs, ys, avg = load_accesses_from_file(path)

    if len(xs) > 0:
        xs_all_union.update(xs)
        plt.plot(
            xs,
            ys,
            marker=marker,
            linewidth=2.5,
            color=color,
            label=f"{label} (avg={avg:.3f})"
        )

def human_format(y, pos):
    if y >= 1e9:
        return f'{y/1e9:.1f}B'
    elif y >= 1e6:
        return f'{y/1e6:.1f}M'
    elif y >= 1e3:
        return f'{y/1e3:.1f}K'
    else:
        return f'{y:.0f}'


plt.title(
    "Baseline-LISA vs Segment=ed-LISA: Cumulative Index Accesses (GRCh38)",
    fontsize=16
)
plt.xlabel("Number of accesses / range adjustments", fontsize=16)
plt.ylabel("Cumulative accesses to index", fontsize=16)
plt.grid(True, linestyle='--', alpha=0.6)

ax = plt.gca()
ax.yaxis.set_major_formatter(FuncFormatter(human_format))

if xs_all_union:
    xs_all = sorted(xs_all_union)
    plt.xticks(xs_all, [x + 1 for x in xs_all])

plt.legend(fontsize=16)
plt.tight_layout()

os.makedirs('graphs', exist_ok=True)
output_file = "graphs/Baseline_vs_Segmented_Accesses.png"
plt.savefig(output_file, dpi=300)

print(f"\nGráfico guardado en: {output_file}")
