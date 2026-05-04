import ssl
ssl._create_default_https_context = ssl._create_unverified_context

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

cnts = pd.read_csv('./read-counts.txt', sep='\t', comment='#', index_col=0)
cnts['clip_enrichment'] = cnts['./binfo1-datapack1/CLIP-35L33G.bam'] / cnts['./binfo1-datapack1/RNA-control.bam']
cnts['rden_change'] = (cnts['./binfo1-datapack1/RPF-siLin28a.bam'] / cnts['./binfo1-datapack1/RNA-siLin28a.bam']) / (cnts['./binfo1-datapack1/RPF-siLuc.bam'] / cnts['./binfo1-datapack1/RNA-siLuc.bam'])

mouselocal = pd.read_csv('https://hyeshik.qbio.io/binfo/mouselocalization-20210507.txt', sep='\t')

cnts = cnts.copy()
cnts['gene_id'] = cnts.index.str.split('.').str[0]
merged = cnts.merge(mouselocal, on='gene_id', how='inner')

nonzero_idx = (
    ~np.isnan(merged['clip_enrichment'].values)
    & ~np.isnan(merged['rden_change'].values)
    & (merged['clip_enrichment'].values > 0)
    & (merged['rden_change'].values > 0)
)
merged = merged[nonzero_idx]

color_map = {'nucleus': 'blue', 'integral membrane': 'red', 'cytoplasm': 'green'}

fig, ax = plt.subplots(1, 1, figsize=(5, 5))
for loc_type, color in color_map.items():
    subset = merged[merged['type'] == loc_type]
    ax.scatter(np.log2(subset['clip_enrichment'].values),
               np.log2(subset['rden_change'].values),
               s=2, alpha=0.6, color=color, label=loc_type)

ax.set_xlabel('LIN28A CLIP enrichment (log$_2$)')
ax.set_ylabel('Ribosome density change\nupon $\\mathit{Lin28a}$ knockdown (log$_2$)')
ax.grid(True)
ax.legend(loc='best', markerscale=3)
plt.tight_layout()
fig.savefig('./figure5B.png', dpi=300)
plt.close()
