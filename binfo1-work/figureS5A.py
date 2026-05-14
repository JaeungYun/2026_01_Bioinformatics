import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator

START_LEFT, START_RIGHT = -50, 50
STOP_LEFT, STOP_RIGHT = -50, 20

def accumulate_counts(path, codon_col, left, right):
    counts = np.zeros(right - left + 1, dtype=np.int64)
    seen = set()
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            pos = int(parts[1])
            cnt = int(parts[3])
            tx = parts[7]
            codon = int(parts[codon_col])
            rel = pos - codon
            if rel < left or rel > right:
                continue
            key = (tx, pos)
            if key in seen:
                continue
            seen.add(key)
            counts[rel - left] += cnt
    return counts

samples = [('RPF-siLuc', 'siLuc'), ('RPF-siLin28a', 'siLin28a')]

data = {}
for sample, _ in samples:
    data[(sample, 'start')] = accumulate_counts(
        f'fivepcounts-startfiltered-{sample}.txt', 8, START_LEFT, START_RIGHT) / 1000.0
    data[(sample, 'stop')] = accumulate_counts(
        f'fivepcounts-stopfiltered-{sample}.txt', 8, STOP_LEFT, STOP_RIGHT) / 1000.0

ymax = max(arr.max() for arr in data.values())
ymax *= 1.05

x_start = np.arange(START_LEFT, START_RIGHT + 1)
x_stop = np.arange(STOP_LEFT, STOP_RIGHT + 1)

fig, axes = plt.subplots(2, 2, figsize=(12, 6),
                         gridspec_kw={'width_ratios': [(START_RIGHT - START_LEFT + 1),
                                                       (STOP_RIGHT - STOP_LEFT + 1)]})

for row, (sample, label) in enumerate(samples):
    ax_left = axes[row, 0]
    ax_right = axes[row, 1]

    ax_left.bar(x_start, data[(sample, 'start')], width=1.0, color='black')
    ax_left.axvline(0, color='red', linewidth=1)
    ax_left.set_xlim(START_LEFT, START_RIGHT)
    ax_left.set_ylim(0, ymax)
    ax_left.xaxis.set_major_locator(MultipleLocator(10))
    ax_left.set_xlabel("Relative position to start codon of 5'-end of reads")
    ax_left.set_ylabel(f"$\\bf{{{label}}}$\nRaw read count\n(x1000)")
    ax_left.text(0, ymax, ' start codon', color='red',
                 ha='left', va='top', fontsize=10)

    ax_right.bar(x_stop, data[(sample, 'stop')], width=1.0, color='black')
    ax_right.axvline(0, color='red', linewidth=1)
    ax_right.set_xlim(STOP_LEFT, STOP_RIGHT)
    ax_right.set_ylim(0, ymax)
    ax_right.xaxis.set_major_locator(MultipleLocator(10))
    ax_right.set_xlabel("Relative position to stop codon of 5'-end of reads")
    ax_right.text(0, ymax, ' stop codon', color='red',
                  ha='left', va='top', fontsize=10)

for ax in axes.ravel():
    ax.grid(True, alpha=0.3)

plt.tight_layout()
fig.savefig('figureS5A.png', dpi=300)
plt.close()
