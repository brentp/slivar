import matplotlib
matplotlib.use("Agg")
import seaborn as sns
from matplotlib import pyplot as plt

import sys
import pandas as pd
import numpy as np

col_names = ["impact", "count"]

ex = pd.read_csv(sys.argv[1], sep="\t", names=col_names)
ex["tech"] = "exome"
wg = pd.read_csv(sys.argv[2], sep="\t", names=col_names)
wg["tech"] = "genome"

bo = pd.concat([ex, wg])

bosum = bo.groupby("impact").sum()
sel = bosum["count"] > 20
print(sel.shape)

bosum = bosum[sel]

order = np.array(bosum.sort_values("count", ascending=False).index)


bop = pd.DataFrame(bo.pivot(index="impact", columns="tech", values="count"))

bop.fillna(0, inplace=True)

# now make values scaled since we don't care about actual number

bop["exome"] /= bop.exome.sum()
bop["genome"] /= bop.genome.sum()

bop = bop.loc[order, :]
print(bop)
bop["impact"] = bop.index

mlt = pd.melt(bop, id_vars=["impact"], value_vars=["exome", "genome"])

ax = sns.barplot(data=mlt, x="impact", y="value", hue="tech")
ax.set_ylabel("Proportion of impactful variants")
ax.set_xlabel("Highest impact")

for label in ax.get_xticklabels():
    label.set_ha("right")
    label.set_rotation(25)

plt.tight_layout()

# https://stackoverflow.com/a/51535326
# (Thanks Sharon)
def show_values_on_bars(axs):
    def _show_on_single_plot(ax):
        for p in ax.patches:
            _x = p.get_x() + p.get_width() / 2
            _y = p.get_y() + p.get_height() + 0.005
            value = '{:.2f}'.format(p.get_height())
            ax.text(_x, _y, value, ha="center", size=8, rotation=25)

    if isinstance(axs, np.ndarray):
        for idx, ax in np.ndenumerate(axs):
            _show_on_single_plot(ax)
    else:
        _show_on_single_plot(axs)

show_values_on_bars(ax)
#for index, row in mlt.iterrows():
#    print(row)
#    ax.text(row.name, row.value, "%.2f" % row.value, ha="left")
plt.savefig("rev2q3.png")
