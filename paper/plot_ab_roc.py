import sys

from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd

colors = sns.color_palette("RdYlBu", 13)
colors = sns.cubehelix_palette(7, start=2, rot=0, dark=0.05, light=.95, reverse=False)
colors = sns.color_palette("Set1", 13)

df = pd.read_csv(sys.argv[1], sep="\t")
cols = list(df.columns)
cols[0] = cols[0].lstrip('#')
df.columns = cols

if "genome" in sys.argv[1]:
    df = df[df.exclude == "exclude"]
else:
    df = df[df.exclude == "."]
tmax = float(df.total.max())

fig, axes = plt.subplots(len(df.GQ.unique()), 1, sharex=True, sharey=True, figsize=(5, 12))

gqs = sorted(df.GQ.unique())

for i, gq in enumerate(gqs):
    df_gq = df[df.GQ == gq]
    ax = axes[i]

    # each df has exclude and '' (for no exclusion)
    for k, sel in enumerate(("exclude", ".")):

        dfs = df_gq[df_gq.exclude == sel]
        if len(dfs) == 0: continue
        dfs.tp = dfs.tp * (dfs.total / tmax)

        #if k == 0: continue

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

        """
        if not "genome" in sys.argv[1] and gqs[i] == 5:
            cutoff = dfs[dfs.cutoff > 0.2].iloc[0, :]

            ax.annotate('0.2 <= AB < 0.8', xy=(cutoff.fp, cutoff.tp),
                    xytext=(cutoff.fp+0.0003, cutoff.tp-0.03),
                arrowprops=dict(facecolor='gray', shrink=0.02),
                )
        elif "genome" in sys.argv[1] and gqs[i] == 10:
            cutoff = dfs[dfs.cutoff > 0.22].iloc[0, :]

            ax.annotate('0.25 <= AB < 0.75', xy=(cutoff.fp, cutoff.tp),
                    xytext=(cutoff.fp+0.0003, cutoff.tp-0.03),
                arrowprops=dict(facecolor='gray', shrink=0.02),
                )
        """


        if not "genome" in sys.argv[1]:
            ax.set_xlim(0.0, 0.0024)
            ax.set_ylim(0.8, 1)
        else:
            ax.set_xlim(0.0, 0.006)
            ax.set_ylim(0.5, 1)

        sns.despine()
        if i == 0:
            ax.legend(title="Allele-balance cut-offs")
        ax.set_title("Genotype-quality cutoff: %d" % gqs[i])
        ax.set_ylabel("Transmission rate")

axes[len(axes)-1].set_xlabel("Mendelian-violation rate")
#plt.ylim(0.8, 1)
#sns.despine()
ax = axes[0]
plt.tight_layout()
plt.show()


