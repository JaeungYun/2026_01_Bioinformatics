import os
import re
import subprocess
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

BAM = '../binfo1-datapack1/CLIP-35L33G.bam'

GENES = [
    ('Mirlet7g',   'chr9',  106056039, 106056126, '+'),
    ('Mirlet7d',   'chr13',  48689488,  48689590, '-'),
    ('Mirlet7f-1', 'chr13',  48691305,  48691393, '-'),
]

BASE_RE = re.compile(r'\^.|\$|[<>*#]|[+-](\d+)')

def clean_pileup_bases(s):
    out = []
    i = 0
    while i < len(s):
        c = s[i]
        if c == '^':
            i += 2
            continue
        if c == '$':
            i += 1
            continue
        if c in '<>*#':
            i += 1
            continue
        if c in '+-':
            m = re.match(r'(\d+)', s[i+1:])
            if m:
                n = int(m.group(1))
                i += 1 + len(m.group(1)) + n
                continue
        out.append(c)
        i += 1
    return ''.join(out)

def shannon_entropy(bases):
    if not bases:
        return 0.0
    counts = {}
    for b in bases.upper():
        if b in '.,':
            counts['M'] = counts.get('M', 0) + 1
        elif b in 'ACGTN':
            counts[b] = counts.get(b, 0) + 1
    total = sum(counts.values())
    if total == 0:
        return 0.0
    H = 0.0
    for v in counts.values():
        p = v / total
        if p > 0:
            H -= p * np.log2(p)
    return H

def run_pileup(gene, chrom, start, end):
    bam_out = f'CLIP-{gene}.bam'
    pileup_out = f'CLIP-{gene}.pileup'
    gene_pileup = f'CLIP-{gene}-gene.pileup'
    region = f'{chrom}:{start}-{end}'
    subprocess.run(['samtools', 'view', '-b', '-o', bam_out, BAM, region], check=True)
    with open(pileup_out, 'w') as fh:
        subprocess.run(['samtools', 'mpileup', bam_out], check=True, stdout=fh,
                       stderr=subprocess.DEVNULL)
    df = pd.read_csv(pileup_out, sep='\t',
                     names=['chrom', 'pos', '_ref', 'count', 'basereads', 'quals'],
                     dtype={'chrom': str, 'basereads': str, 'quals': str})
    df = df[(df['pos'] >= start) & (df['pos'] <= end)].copy()
    df.to_csv(gene_pileup, sep='\t', header=False, index=False)
    return df

def make_entropy_table(df):
    df = df.copy()
    df['matches'] = df['basereads'].fillna('').apply(clean_pileup_bases)
    df['entropy'] = df['matches'].apply(shannon_entropy)
    return df[['chrom', 'pos', 'count', 'matches', 'entropy']]

def write_bedgraph(df, gene, chrom, start, end):
    path = f'CLIP-{gene}.bedgraph'
    with open(path, 'w') as fh:
        fh.write(f'track type=bedGraph name="CIMS-entropy-{gene}" '
                 f'description="Shannon entropy of CLIP-35L33G at {gene}" '
                 f'visibility=full color=0,0,0\n')
        positions = dict(zip(df['pos'].astype(int), df['entropy']))
        for p in range(start, end + 1):
            H = positions.get(p, 0.0)
            fh.write(f'{chrom}\t{p-1}\t{p}\t{H:.6f}\n')
    return path

def main():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    results = {}
    for gene, chrom, start, end, strand in GENES:
        print(f'[{gene}] {chrom}:{start}-{end} ({strand})')
        pileup = run_pileup(gene, chrom, start, end)
        entropy_df = make_entropy_table(pileup)
        write_bedgraph(entropy_df, gene, chrom, start, end)
        results[gene] = (chrom, start, end, strand, entropy_df)

    fig, axes = plt.subplots(len(GENES), 1, figsize=(10, 2.2 * len(GENES)), sharey=True)
    if len(GENES) == 1:
        axes = [axes]
    for ax, (gene, chrom, start, end, strand) in zip(axes, GENES):
        _, _, _, _, df = results[gene]
        positions = np.arange(start, end + 1)
        entropy = np.zeros_like(positions, dtype=float)
        pos_map = dict(zip(df['pos'].astype(int), df['entropy']))
        for i, p in enumerate(positions):
            entropy[i] = pos_map.get(p, 0.0)
        ax.bar(positions, entropy, width=1.0, color='black')
        ax.set_xlim(start, end)
        ax.set_ylabel(f'$\\bf{{{gene}}}$\nShannon entropy')
        ax.set_title(f'{gene} — {chrom}:{start}-{end} ({strand})', fontsize=9, loc='left')
        ax.grid(True, alpha=0.3)
    axes[-1].set_xlabel('Genomic position (mm39)')
    plt.tight_layout()
    fig.savefig('figureS2A.png', dpi=300)
    plt.close()
    print('Wrote figureS2A.png')

if __name__ == '__main__':
    main()
