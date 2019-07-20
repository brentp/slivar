import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

import sys

colors = sns.color_palette("Set1", 12)

df = pd.read_csv(sys.argv[1], sep="\t", index_col=0)
dfo = df.join(pd.read_csv(sys.argv[2], sep="\t", index_col=0))
dfo.drop("comphet_side", inplace=True, axis="columns")

print(dfo.head())

df = dfo.melt(value_name="number_of_variants")



fig, axes = plt.subplots(1, figsize=(7, 7))
axes = [axes]
sns.swarmplot(x="variable", y="number_of_variants",
       data=df, ax=axes[0], palette=colors)

axes[0].set_xlabel("Inheritance mode")
axes[0].set_ylabel("Candidate variants")
#axes[0].set_xticklabels(["0.2 <= AB < 0.8\n& GQ >= 5",
#                         "gnomAD popmax AF < 0.001\n& PASS in gnomAD\n & topMed AF < 0.01\n"])

print(axes[0].get_xlim())

ax = axes[0]
#ax.set_ylim(-0.5, 14.5)
#ax.axhline(o, 0.1, 0.4, lw=4)
#ax.axhline(f, 1-0.4, 1 - 0.1, lw=4)

plt.tight_layout()
sns.despine()
plt.show()




