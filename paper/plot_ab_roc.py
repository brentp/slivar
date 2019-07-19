import sys

from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd

colors = sns.color_palette("RdYlBu", 13)
colors = sns.cubehelix_palette(7, start=2, rot=0, dark=0.05, light=.95, reverse=False)
colors = sns.color_palette("Set1", 13)

gqs = [int(sys.argv[i].split("gq")[-1].split(".")[0]) for i in range(1, len(sys.argv))]

dfs = [pd.read_csv(sys.argv[i], sep="\t") for i in range(1, len(sys.argv))]
tmax = 0
for df in dfs:
    cols = list(df.columns)
    cols[0] = cols[0].lstrip('#')
    df.columns = cols
    tmax = max(tmax, df.total.max())

tmax = float(tmax)
fig, axes = plt.subplots(len(dfs), 1, sharex=True, sharey=True, figsize=(5, 12))
print(axes)

for i, df in enumerate(dfs):
    df.tp = df.tp * (df.total / tmax)

    ax = axes[i]

    cutoffs = [0, 0.1, 0.15, 0.2, 0.25, 0.3]
    for j, cutoff in enumerate(cutoffs):
        if j < len(cutoffs) - 1:
            c1 = cutoffs[j + 1]
        else:
            c1 = 1
        cut = df[(df.cutoff >= cutoff) & (df.cutoff < c1)]
        ax.scatter(cut.fp, cut.tp, s=13, ec=colors[j+1], color=colors[j+1], label="%.2g-%.2g" % (cutoff, 1-cutoff))

    if gqs[i] == 5:
        cutoff = df[df.cutoff > 0.2].iloc[0, :]

        ax.annotate('0.2 <= AB < 0.8', xy=(cutoff.fp, cutoff.tp),
                xytext=(cutoff.fp+0.0003, cutoff.tp-0.03),
            arrowprops=dict(facecolor='gray', shrink=0.02),
            )


    ax.set_xlim(0.0, 0.003)
    ax.set_ylim(0.8, 1)
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


