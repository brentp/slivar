import sys

from matplotlib import pyplot as plt
import seaborn as sns
sns.set_style("white")
import pandas as pd

colors = sns.color_palette("RdYlBu", 13)
colors = sns.cubehelix_palette(7, start=2, rot=0, dark=0.05, light=.95, reverse=False)
colors = sns.color_palette("Set1", 13)

df_exome = pd.read_csv(sys.argv[1], sep="\t")
df_genome = pd.read_csv(sys.argv[2], sep="\t")

for df in [df_exome, df_genome]:
    cols = list(df.columns)
    cols[0] = cols[0].lstrip('#')
    df.columns = cols


fig, axes = plt.subplots(len(df.GQ.unique()), 2, figsize=(8, 12))
assert len(df_exome.GQ.unique()) == len(df_genome.GQ.unique())

for ci, df in enumerate((df_exome, df_genome)):
    tmax = float(df.total.max())
    gqs = sorted(df.GQ.unique())

    for i, gq in enumerate(gqs):
        df_gq = df[df.GQ == gq]
        ax = axes[i, ci]
        dfs = df_gq
        if len(dfs) == 0: continue
        dfs.tp = dfs.tp * (dfs.total / tmax)

        cutoffs = [0.1, 0.15, 0.2, 0.25, 0.3]
        ax.plot(dfs.fp, dfs.tp, color=(0.6, 0.6, 0.6))
        for j, cutoff in enumerate(cutoffs):
            if j < len(cutoffs) - 1:
                c1 = cutoffs[j + 1]
            else:
                c1 = 1
            cut = dfs[(dfs.cutoff >= cutoff) & (dfs.cutoff < c1)].iloc[0, :]
            ax.plot([cut.fp], [cut.tp], color=colors[j+1], ls='none',
                    marker='o', label="%.2g-%.2g" % (cutoff, 1 - cutoff))

            if ci == 0 and cutoff == 0.2 and gq == 10:
                ax.plot([cut.fp], [cut.tp], color='gray', ls='none',
                        markersize=12,
                        marker='o', zorder=-1)
            elif ci == 1 and cutoff == 0.25 and gq == 10:
                ax.plot([cut.fp], [cut.tp], color='gray', ls='none',
                        markersize=12,
                        marker='o', zorder=-1)

        if ci == 0:
            ax.set_xlim(0.0, 0.0015)
            ax.set_ylim(0.85, 1)
        else:
            ax.set_xlim(0.0, 0.006)
            ax.set_ylim(0.6, 1)

        sns.despine()
        if i == 0 and ci == 0:
            ax.legend(title="Allele-balance cut-offs")
        ax.set_title("Genotype-quality cutoff: %d" % gqs[i])
        ax.set_ylabel("Transmission rate")

        if i != len(gqs) - 1:
            ax.set_xticklabels([])

axes[len(axes)-1, 0].set_xlabel("Mendelian-violation rate")
axes[len(axes)-1, 1].set_xlabel("Mendelian-violation rate")


axes[0, 0].set_title("Exome", fontsize=15)
axes[0, 1].set_title("Genome", fontsize=15)

#plt.ylim(0.8, 1)
#sns.despine()
ax = axes[0]
plt.tight_layout()
plt.show()


