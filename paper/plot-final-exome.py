import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt

import sys

colors2 = sns.color_palette("Set1", 12)[3:]

colors = [(0.6, 0.6, 0.6)] * 12


def read_dfs(i, j):

    df = pd.read_csv(sys.argv[i], sep="\t", index_col=0)
    dfo = df.join(pd.read_csv(sys.argv[j], sep="\t", index_col=0))
    dfo.drop("comphet_side", inplace=True, axis="columns")
    ad = dfo.auto_dom
    ad = pd.DataFrame({"number_of_variants": ad})
    dfo.drop("auto_dom", inplace=True, axis="columns")
    df = dfo.melt(value_name="number_of_variants")
    df.number_of_variants += np.random.random(size=df.shape[0]) / 2.0
    return df, ad

def get_mean_df(df_means, order):

    df_means = df_means.to_dict(orient="dict")["number_of_variants"]

    df_means_o = [df_means[(l, "all")] for l in order]
    df_means_i = [df_means[(l, "impactful")] for l in order]

    df_means_b = [[labels[i], o, "all"] for i, o in enumerate(df_means_o)] + \
                 [[labels[i], o, "impactful"] for i, o in enumerate(df_means_i)]

    df_means_bar = pd.DataFrame(df_means_b, columns=["inheritance mode", "mean number of variants",
        "variant class"])
    return df_means_bar


df, ad = read_dfs(1, 2)



dfi, adi = read_dfs(3, 4)

df["cat"] = "all"
dfi["cat"] = "impactful"

ad["cat"] = "all"
adi["cat"] = "impactful"

df = pd.concat((df, dfi))
ad = pd.concat((ad, adi), ignore_index=True)
ad["variable"] = "autosomal dominant"

print(df.head())
#print(ad.head())


df_means = df.groupby(["variable", "cat"]).mean()
ad_means = ad.groupby(["variable", "cat"]).mean()
print(df_means, ad_means)


size = 2

fig, axes = plt.subplots(2, 2, figsize=(8, 5), gridspec_kw={'width_ratios': [5,
    1], 'height_ratios': [1, 3]})

sns.swarmplot(x="variable", y="number_of_variants", hue="cat", dodge=True,
       data=df, ax=axes[1,0], palette=colors2, size=size)

sns.swarmplot(x="variable", y="number_of_variants", data=ad, ax=axes[1,1], hue="cat", dodge=True,
        palette=colors2, size=size)

axes[1,1].set_ylabel("")
axes[1,0].set_xlabel("Inheritance mode", horizontalalignment='left')
axes[1,1].set_xlabel("")
#plt.xlabel("Inheritance mode")
axes[1,0].set_ylabel("Candidate variants")

ax = axes[1,0]

labels = ax.get_xticklabels()
lookups = {"denovo":"de novo", "compound-het":"compound-heterozygote",
        "x_recessive":"x-linked recessive", "recessive": "autosomal recessive",
        "auto_dom":"autosomal dominant",
        "x_denovo": "x-linked de novo"}
order = [l.get_text() for l in labels]
labels = [lookups[l.get_text()] for l in labels]

ax.set_xticklabels(labels, rotation=12)
axes[1,1].set_xticklabels(["autosomal dominant"], rotation=12)

axes[1,1].get_legend().remove()
axes[1,0].get_legend().set_title("variant class")


df_means_bar = get_mean_df(df_means, order)


sns.barplot(data=df_means_bar, hue="variant class", x="inheritance mode",
        y="mean number of variants", ax=axes[0, 0], palette=colors2)
axes[0, 0].set_ylabel(None)
axes[0, 0].set_xlabel(None)
axes[0, 0].set_xticklabels([])
axes[0, 0].set_title("mean number of candidate variants per trio")
axes[0, 0].get_legend().remove()


order = ["autosomal dominant"]
ad_means_bar = get_mean_df(ad_means, ["autosomal dominant"])


sns.barplot(data=ad_means_bar, hue="variant class", x="inheritance mode",
        y="mean number of variants", ax=axes[0, 1], palette=colors2)
axes[0, 1].set_ylabel(None)
axes[0, 1].set_xlabel(None)
axes[0, 1].set_xticklabels([])
axes[0, 1].get_legend().remove()


plt.tight_layout()
sns.despine()
plt.savefig("figure3-exome-counts.png")
plt.show()
