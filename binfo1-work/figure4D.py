import pandas as pd
cnts = pd.read_csv('./read-counts.txt', sep='\t', comment='#', index_col=0)
cnts.head()

cnts['clip_enrichment'] = cnts['./binfo1-datapack1/CLIP-35L33G.bam'] / cnts['./binfo1-datapack1/RNA-control.bam']
cnts['rden_change'] = (cnts['./binfo1-datapack1/RPF-siLin28a.bam'] / cnts['./binfo1-datapack1/RNA-siLin28a.bam']) / (cnts['./binfo1-datapack1/RPF-siLuc.bam'] / cnts['./binfo1-datapack1/RNA-siLuc.bam'])
cnts.head()

from matplotlib import pyplot as plt
import numpy as np

fig, ax = plt.subplots(1, 1, figsize=(5, 5))
ax.scatter(np.log2(cnts['clip_enrichment']),
           np.log2(cnts['rden_change']))
plt.savefig("./figure4D_deprecated.png")
plt.close()

fig, ax = plt.subplots(1, 1, figsize=(5, 5))
nonzero_clip_idx = ~np.isnan(cnts['clip_enrichment'].values) & (cnts['rden_change'].values > 0)
nonzero_rden_idx = ~np.isnan(cnts['rden_change'].values) & (cnts['clip_enrichment'].values > 0)
nonzero_idx = nonzero_clip_idx & nonzero_rden_idx
ax.scatter(np.log2(cnts['clip_enrichment'].values[nonzero_idx]),np.log2(cnts['rden_change'].values[nonzero_idx]), s=1, alpha=0.5, color='black')

ax.set_xlabel('LIN28A CLIP enrichment (log$_2$)')
ax.set_ylabel('Ribosome density change\nupon $\\mathit{Lin28a}$ knockdown (log$_2$)')
ax.grid(True)
plt.tight_layout()
fig.savefig('./figure4D_revised.png', dpi=300)
plt.close()