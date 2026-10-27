import os
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

import matplotlib
matplotlib.use("Agg")

import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from scipy.ndimage import gaussian_filter1d
import scanpy as sc

python_data_dir = rf"C:\Users\Nikla\Documents\MoBi_Master_FS.3\Internship AG Höfer\DynaVelo_testing\data"
python_plot_dir = rf"C:\Users\Nikla\Documents\MoBi_Master_FS.3\Internship AG Höfer\DynaVelo_testing\plots"

FRAGMENTS = rf"{python_data_dir}\atac_fragments.tsv.gz"
ATAC_H5AD = rf"{python_data_dir}\all_atac_peaks_corr.h5ad"
RNA_H5AD = rf"{python_data_dir}\total_counts_raw.h5ad"
GTF = rf"{python_data_dir}\genes.gtf.gz"
BINS_CSV = rf"{python_data_dir}\cell_bins.csv"
PEAKS_CSV = rf"{python_data_dir}\Exp4_ctr_peak_types.csv"

N_GENES, N_BINS, WINDOW, SMOOTH = 5, 10, 100, 50   # SMOOTH: gaussian sigma in bp
UP, DOWN = 2000, 100
SEED = 0

gene_data = pd.read_csv(rf"{python_data_dir}\pmtr_peak_sub\chrom_velo_all_genes.csv", index_col=0, header=0)
gene_subset = gene_data[(gene_data["monotonicity_r"] <= 0.1) 
          & (gene_data["monotonicity_a"] <= 0.1) 
          & (gene_data["rna_dir"] == 1)
          & (gene_data["atac_dir"] == -1)]
          #& (gene_data["rna_dir"] == gene_data["atac_dir"])]

if len(gene_subset) < N_GENES:
    N_GENES = len(gene_subset)
    print(f"Set N_GENES to {len(gene_subset)} as no more available.")

GENE_LIST = gene_subset.index.tolist() # rf"{python_data_dir}\genes.txt"

print("Loading all data...")
bins = pd.read_csv(BINS_CSV, header=0, index_col=0)["bin"]

atac = sc.read_h5ad(ATAC_H5AD)
rna = sc.read_h5ad(RNA_H5AD)

rna.obs_names = [i[9:-2] for i in rna.obs_names.values]

bins = bins[bins.index.isin(atac.obs_names)]
rna = rna[bins.index.intersection(rna.obs_names)]

bin_of = {f"{bc}-1": b for bc, b in bins.items()}   # fragment-style barcodes -> bin
barcodes = set(bin_of)

cells_per_bin = bins.value_counts().reindex(range(N_BINS))

first_chrom = str(pd.read_csv(FRAGMENTS, sep="\t", comment="#", header=None, nrows=1)[0].iloc[0])
prefix = "chr" if first_chrom.startswith("chr") else ""

gtf = pd.read_csv(GTF, sep="\t", comment="#", header=None, usecols=[0, 2, 3, 4, 6, 8],
                  names=["chrom", "feature", "start", "end", "strand", "attr"])
gtf = gtf[gtf.feature == "gene"].copy()
gtf["gene"] = gtf.attr.str.extract(r'gene_name "([^"]+)"')
gtf["chrom"] = prefix + gtf.chrom.astype(str).str.replace("^chr", "", regex=True)
genes = gtf.dropna(subset=["gene"]).drop_duplicates("gene").set_index("gene")

peaks = pd.read_csv(PEAKS_CSV, index_col=0)
peaks["chromosome"] = prefix + peaks.chromosome.astype(str).str.replace("^chr", "", regex=True)
peaks["group"] = np.where(peaks.peakType == "Promoter", "Promoter", "other")

print("Loading functions...")

FRAC = 0.05   # fraction of fragments to keep; 1.0 = everything

def load_fragments(path, regions):
    by_chrom = {}
    for gene, (c, s, e) in regions.items():
        by_chrom.setdefault(c, []).append((gene, s, e))
    out = {gene: [] for gene in regions}
    rng = np.random.default_rng(SEED)

    reader = pd.read_csv(path, sep="\t", comment="#", header=None, usecols=[0, 1, 2, 3],
                         names=["chrom", "start", "end", "bc"], chunksize=1_000_000)
    for chunk in reader:
        chunk = chunk[chunk.chrom.isin(by_chrom) & chunk.bc.isin(barcodes)]
        if FRAC < 1:
            chunk = chunk[rng.random(len(chunk)) < FRAC]
        for c, items in by_chrom.items():
            sub = chunk[chunk.chrom == c]
            for gene, s, e in items:
                out[gene].append(sub[(sub.end >= s) & (sub.start <= e)])
    return {gene: pd.concat(v) for gene, v in out.items()}


def coverage(frag, start, end):
    # Signac: 1-bp cut counts per group -> count / n_cells * median(n_cells) -> centered rolling mean -> gaussian smoothing
    pos = np.arange(start, end + 1)
    edges = np.arange(start, end + 2)
    cov = np.zeros((N_BINS, len(pos)))
    for b, grp in frag.assign(bin=frag.bc.map(bin_of)).groupby("bin"):
        cov[int(b)] = np.histogram(np.r_[grp.start, grp.end], bins=edges)[0]
    cov = cov / cells_per_bin.values[:, None] * np.median(cells_per_bin.values) / FRAC

    roll = pd.DataFrame(cov.T, index=pos).rolling(WINDOW, center=True).mean().dropna()
    return roll.index.values, gaussian_filter1d(roll.values.T, SMOOTH, axis=1, mode="nearest")


def expr_per_bin(gene):
    x = rna[:, gene].X
    x = x.toarray().ravel() if sp.issparse(x) else np.asarray(x).ravel()
    s, g = pd.Series(x, index=rna.obs_names), bins.reindex(rna.obs_names)
    return [s[g == b].values for b in range(N_BINS)]

print("Defining gene candidates...")

cmap, norm = plt.get_cmap("plasma"), colors.Normalize(0, N_BINS - 1)
rng = np.random.default_rng(SEED)
candidates = GENE_LIST
selected = rng.choice(candidates, min(N_GENES, len(candidates)), replace=False)

regions = {}
for gene in selected:
    g = genes.loc[gene]
    left, right = (UP, DOWN) if g.strand == "+" else (DOWN, UP)
    regions[gene] = (g.chrom, max(int(g.start) - left, 0), int(g.end) + right)

print("Loading ATAC fragments...")

fragments = load_fragments(FRAGMENTS, regions)

print("Visualizing results...")
for gene in selected:
    g = genes.loc[gene]
    chrom, start, end = regions[gene]
    x, cov = coverage(fragments[gene], start, end)
    expr = expr_per_bin(gene)
    ymax = float(f"{cov.max():.2g}")

    fig, axs = plt.subplots(3, 2, figsize=(11, 5), layout="constrained",
                            gridspec_kw={"width_ratios": [3, 1], "height_ratios": [8, 1, 0.6]})
    (ax1, ax2), (ax3, ax4), (ax5, ax6) = axs
    for b in range(N_BINS):
        ax1.plot(x, cov[b], color=cmap(norm(b)), lw=1.5, alpha=0.55)
    top = cov.argmax(axis=1)
    # ax1.plot(x[top], cov[range(N_BINS), top], color="grey", lw=1, zorder=2)
    ax1.scatter(x[top], cov[range(N_BINS), top], c=range(N_BINS), cmap=cmap, norm=norm,
                 s=30, edgecolor="k", linewidths=0.5, zorder=3)
    ax1.set(title=f"{gene} ({g.strand})", xlim=(start, end), ylim=(0, ymax+0.1), yticks=[],
            ylabel=f"Accessibility \n(range 0 - {ymax:g})")
    ax1.tick_params(labelbottom=False)

    gcol = "darkblue" if g.strand == "+" else "darkgreen"
    tail, head = (g.start, g.end) if g.strand == "+" else (g.end, g.start)
    ax3.annotate("", xy=(head, 0), xytext=(tail, 0), arrowprops=dict(arrowstyle="-|>", color=gcol, lw=2))
    ax3.text((g.start + g.end) / 2, 0.3, gene, ha="center", va="bottom", fontsize=8, color=gcol)
    ax3.set(xlim=(start, end), ylim=(-1, 1), yticks=[])
    ax3.tick_params(labelbottom=False)

    pk = peaks[(peaks.chromosome == chrom) & (peaks.stop >= start) & (peaks.start <= end)]
    for t, c in {"Promoter": "darkred", "other": "grey"}.items():
        sub = pk[pk.group == t]
        ax5.hlines(np.zeros(len(sub)), sub.start, sub.stop, color=c, lw=4, label=t)
    ax5.set(xlim=(start, end), ylim=(-1, 1), yticks=[], ylabel="Peaks", xlabel=f"{chrom} position (bp)")
    ax6.legend(*ax5.get_legend_handles_labels(), loc="center left", frameon=False)

    vp = ax2.violinplot(expr, positions=range(N_BINS), vert=False, widths=0.9, showextrema=False)
    for b, body in enumerate(vp["bodies"]):
        body.set_facecolor(cmap(norm(b)))
        body.set_alpha(0.8)
    means = [e.mean() for e in expr]
    ax2.plot(means, range(N_BINS), color="grey", lw=1, zorder=2)
    ax2.scatter(means, range(N_BINS), c=range(N_BINS), cmap=cmap, norm=norm, s=30,
                edgecolor="k", linewidths=0.5, zorder=3)
    # plot percent change
    maxes = np.array([e.max() for e in expr])
    pad = 0.03 * maxes.max()
    with np.errstate(divide="ignore", invalid="ignore"):
        pct = (np.array(means[1:]) - np.array(means[:-1])) / np.array(means[:-1]) * 100
    for i, p in enumerate(pct):
        ax2.text(max(maxes[i], maxes[i + 1]) + pad, i + 0.5,
                 f"{p:+.0f}%" if np.isfinite(p) else "n/a",
                 ha="left", va="center", fontsize=7)
        
    ax2.set_xlim(right=maxes.max() * 1.35)
    ax2.invert_yaxis()
    ax2.set(xlabel="mRNA", yticks=[])
    ax4.axis("off")
    ax6.axis("off")

    fig.colorbar(cm.ScalarMappable(norm, cmap), ax=axs, label="Cell bin along PC1")
    print(rf"Saving plot to: {python_data_dir}\{gene}_atac_rna.png !")

    fig.savefig(rf"{python_plot_dir}\{gene}_atac_rna.png")
    fig.savefig(rf"{python_plot_dir}\{gene}_atac_rna.svg")
    plt.close(fig)