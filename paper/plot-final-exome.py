import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt

import sys

colors = sns.color_palette("Set1", 12)

colors = [(0.6, 0.6, 0.6)] * 12

df = pd.read_csv(sys.argv[1], sep="\t", index_col=0)
dfo = df.join(pd.read_csv(sys.argv[2], sep="\t", index_col=0))
dfo.drop("comphet_side", inplace=True, axis="columns")

ad = dfo.auto_dom
ad = pd.DataFrame({"auto_dom": ad})

dfo.drop("auto_dom", inplace=True, axis="columns")

df = dfo.melt(value_name="number_of_variants")

df.number_of_variants += np.random.random(size=df.shape[0]) / 2.5
size=3


fig, axes = plt.subplots(1, 2, figsize=(7, 7), gridspec_kw={'width_ratios': [3, 1]})
#axes = [axes]
sns.swarmplot(x="variable", y="number_of_variants",
       data=df, ax=axes[0], palette=colors, size=size)

print(ad.head())

sns.swarmplot(y="auto_dom", data=ad, ax=axes[1],
        palette=colors, size=size)

axes[1].set_ylabel("")
axes[0].set_xlabel("Inheritance mode", horizontalalignment='left')
axes[1].set_xlabel("")
#plt.xlabel("Inheritance mode")
axes[0].set_ylabel("Candidate variants")

ax = axes[0]

labels = ax.get_xticklabels()
lookups = {"denovo":"de novo", "compound-het":"compound-heterozygote",
        "x_recessive":"x-linked recessive", "recessive": "autosomal recessive",
        "auto_dom":"autosomal dominant",
        "x_denovo": "x-linked de novo"}
labels = [lookups[l.get_text()] for l in labels]

ax.set_xticklabels(labels, rotation=12)
axes[1].set_xticklabels(["autosomal dominant"], rotation=12)

plt.tight_layout()
sns.despine()
plt.show()
