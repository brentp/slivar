import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

import sys
import numpy as np

colors = sns.color_palette("Set1", 12)

c = 0.5
colors = [(c, c, c), (c, c, c)]

dfo = pd.read_csv(sys.argv[1], sep="\t", index_col=0)

df = dfo.melt(value_name="number_of_variants")


df.number_of_variants += np.random.uniform(low=-0.24, high=0.24, size=len(df.variable))


fig, axes = plt.subplots(1, figsize=(9, 7))
plt.subplots_adjust(wspace=0.5)
axes = [axes]
sns.swarmplot(x="variable", y="number_of_variants",
       data=df, ax=axes[0], palette=colors)

axes[0].set_xlabel("Filtering strategy", size=15)
axes[0].set_ylabel("Number of candidate $\it{de}$ $\it{novo}$ variants")
labels = ["0.2 <= AB < 0.8\n& GQ >= 20",
          "gnomAD popmax AF < 0.001",
                         "impactful"]


for i, l in enumerate(labels[1:], start=1):
    labels[i] = labels[i - 1] + "\n& " + l



axes[0].annotate("additive filtering", [2, 7], [0.4, 7],
        ha="center",
        arrowprops=dict(arrowstyle="simple"))


axes[0].set_xticklabels(labels)

axes[0].set_yticks(range(11))

print(axes[0].get_xlim())

ax = axes[0]
ax.set_ylim(-0.5, 10.5)
#ax.axhline(o, 0.1, 0.4, lw=4)
#ax.axhline(f, 1-0.4, 1 - 0.1, lw=4)

sns.despine()
plt.tight_layout()
plt.savefig("figure2-exome-denovos.png")
plt.savefig("figure2-exome-denovos.eps")
plt.show()




