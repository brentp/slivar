import sys

from matplotlib import pyplot as plt
import seaborn as sns
sns.set_style("white")
import pandas as pd

colors = sns.color_palette("RdYlBu", 13)
colors = sns.cubehelix_palette(7, start=2, rot=0, dark=0.05, light=.95, reverse=False)
colors = sns.color_palette("Set1", 13)



df_exome = pd.read_csv(sys.argv[1], sep="\t")
#df_genome = pd.read_csv(sys.argv[2], sep="\t")

for df in [df_exome]:
    cols = list(df.columns)
    cols[0] = cols[0].lstrip('#')
    df.columns = cols

gqs = [x for x in df_exome.GQ.unique() if x != 1]

fig, axes = plt.subplots(1, len(gqs), figsize=(8, 3), sharey=True, sharex=True)
#assert len(df_exome.GQ.unique()) == len(df_genome.GQ.unique())

for ci, df in enumerate((df_exome, )):
    tmax = float(df.total.max())

    for i, gq in enumerate(gqs):
        df_gq = df[df.GQ == gq]
        try:
            ax = axes[i, ci]
        except IndexError:
            ax = axes[i]
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

            if cutoff == 0.2 and gq == 20:
                print(sys.argv[ci + 1], cut)
                ax.plot([cut.fp], [cut.tp], color='gray', ls='none',
                        markersize=12,
                        marker='o', zorder=-1)

            if cutoff == 0.2:
                ax.annotate(("FNR:%.4f" % (1 - float(cut.tp))),
                            [cut.fp, cut.tp], [cut.fp-0.0002, cut.tp-0.25],
                            arrowprops=dict(arrowstyle="simple", facecolor='black'),
                            )

        ax.set_xlim(0, 0.0015)

        sns.despine()
        if i == 1 and ci == 0:
            ax.legend(title="Allele-balance cut-offs", fontsize=9)
        ax.set_title("Genotype-quality cutoff: %d" % gqs[i])
        ax.set_ylabel("Transmission rate")

        #if i != len(gqs) - 1:
        #    ax.set_xticklabels([])

try:
    axes[len(axes)-2, 0].set_xlabel("Mendelian-violation rate")
    axes[len(axes)-2, 1].set_xlabel("Mendelian-violation rate")
except IndexError:
    axes[len(axes)-2].set_xlabel("Mendelian-violation rate")


#axes[0, 0].set_title("Exome", fontsize=15)
#axes[0, 1].set_title("Genome", fontsize=15)
#axes[0, 0].set_title("Deep Variant", fontsize=15)
#axes[0, 1].set_title("GATK", fontsize=15)

#plt.ylim(0.8, 1)
sns.despine()
ax = axes[0]
plt.tight_layout(rect=(0, 0.005, 1, 0.995))
plt.savefig("figure1-exome-GQ-AB-cutoffs.eps")
plt.savefig("figure1-exome-GQ-AB-cutoffs.png")
plt.show()


