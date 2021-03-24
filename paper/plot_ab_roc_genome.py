import sys
import re
import matplotlib
matplotlib.use("Agg")

from matplotlib import pyplot as plt
import seaborn as sns
sns.set_style("white")
import pandas as pd

from matplotlib.ticker import FormatStrFormatter


colors = sns.color_palette("RdYlBu", 13)
colors = sns.cubehelix_palette(7, start=2, rot=0, dark=0.05, light=.95, reverse=False)
colors = sns.color_palette("Set1", 13)

title = sys.argv[1]

args = sys.argv[2:]
args.sort(key=lambda f: int(re.search("-d(\d+)-", f).groups()[0]))

df_genomes = [pd.read_csv(f, sep="\t") for f in args]
#df_genome = pd.read_csv(sys.argv[2], sep="\t")

axis_label_size=13

for df in df_genomes:
    cols = list(df.columns)
    cols[0] = cols[0].lstrip('#')
    df.columns = cols

gqs = list(sorted([x for x in df_genomes[0].GQ.unique() if x > 0]))
print(args)
print(gqs)

fig, axes = plt.subplots(len(df_genomes), len(gqs), figsize=(12, 8), sharey=True, sharex=True)
#assert len(df_genome.GQ.unique()) == len(df_genome.GQ.unique())
tmax = 0

for ci, df in enumerate(df_genomes):
    tmax = max(tmax, float(df.total.max()))

for ci, df in enumerate(df_genomes):

    depth = int(re.search("-d(\d+)-", args[ci]).groups()[0])

    for i, gq in enumerate(gqs):
        df_gq = df[df.GQ == gq]
        ax = axes[ci, i]
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.3f'))


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
            if cutoff == 0.2:
                ax.text(0.15, 0.1, "FP: %d TP: %d" % (cut.n_fps, cut.n_tps),
                        transform=ax.transAxes,
                        color=colors[j+1])

            if cutoff == 0.2:
                ax.plot([cut.fp], [cut.tp], color='gray', ls='none',
                        markersize=12,
                        marker='o', zorder=-1)

                ax.annotate(("FNR:%.4f" % (1 - float(cut.tp))),
                            [cut.fp, cut.tp], [cut.fp-0.0002, cut.tp-0.25],
                            arrowprops=dict(arrowstyle="simple", facecolor='black'))


        sns.despine()
        if i == 1 and ci == 0:
            ax.legend(title="Allele-balance cut-offs", fontsize=9)

        ax.set_title("GQ cutoff: %d depth: %s" % (gqs[i], depth),
                size=axis_label_size)
        ax.set_ylabel("Transmission rate", size=axis_label_size)

        #if i != len(gqs) - 1:
        #    ax.set_xticklabels([])

try:
    for ax in axes[len(axes)-1]:
      ax.set_xlabel("Mendelian-violation rate", size=axis_label_size)
except IndexError:
    axes[len(axes)-1].set_xlabel("Mendelian-violation rate",
            size=axis_label_size)


#axes[0, 0].set_title("genome", fontsize=15)
#axes[0, 1].set_title("Genome", fontsize=15)
#axes[0, 0].set_title("Deep Variant", fontsize=15)
#axes[0, 1].set_title("GATK", fontsize=15)

#plt.ylim(0.8, 1)
sns.despine()
ax = axes[0]
plt.tight_layout(rect=(0, 0.005, 1, 0.995), w_pad=1.8, h_pad=1.8)
plt.savefig("supp-figure1-%s.png" % title.lower().replace(" ", "-"))
plt.savefig("supp-figure1-%s.eps" % title.lower().replace(" ", "-"))
#plt.show()


