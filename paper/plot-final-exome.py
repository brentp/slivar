import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt

import sys

colors2 = sns.color_palette("Set1", 12)[3:]

colors = [(0.6, 0.6, 0.6)] * 12


def read_dfs(i, j, sample_set):

    df = pd.read_csv(sys.argv[i], sep="\t", index_col=0)
    dfo = df.join(pd.read_csv(sys.argv[j], sep="\t", index_col=0))
    dfo.drop("comphet_side", inplace=True, axis="columns")
    dfo = dfo.loc[dfo.index.intersection(sample_set)]
    print(dfo[dfo.auto_dom == 0])

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


ped = sys.argv[5]
# some affected parents
affected_kids = {l.split("\t")[1] for l in open(ped) if l.split("\t")[5] == "2" and l.split("\t")[3] not in "0."}
df, ad = read_dfs(1, 2, affected_kids)

dfi, adi = read_dfs(3, 4, affected_kids)



df["cat"] = "all"
dfi["cat"] = "impactful"

ad["cat"] = "all"
adi["cat"] = "impactful"

df = pd.concat((df, dfi))
ad = pd.concat((ad, adi), ignore_index=True)
ad["variable"] = "autosomal dominant"

df_means = df.groupby(["variable", "cat"]).mean()
ad_means = ad.groupby(["variable", "cat"]).mean()
print(df_means, ad_means)


size = 2

fig, axes = plt.subplots(1, 2, figsize=(8, 5), gridspec_kw={'width_ratios': [5,
    1]})

sns.swarmplot(x="variable", y="number_of_variants", hue="cat", dodge=True,
       data=df, ax=axes[0], palette=colors2, size=size, alpha=0.8)

sns.swarmplot(x="variable", y="number_of_variants", data=ad, ax=axes[1], hue="cat", dodge=True,
        palette=colors2, size=size, alpha=0.8)

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
order = [l.get_text() for l in labels]
labels = [lookups[l.get_text()] for l in labels]

ax.set_xticklabels(labels, rotation=12)
axes[1].set_xticklabels(["autosomal dominant"], rotation=12)

axes[1].get_legend().remove()
axes[0].get_legend().set_title("variant class")


df_means_bar = get_mean_df(df_means, order)
ad_means_bar = get_mean_df(ad_means, ["autosomal dominant"])

xs = axes[0].get_xticks()
for i, y in enumerate(df_means_bar.loc[df_means_bar["variant class"] == "all", "mean number of variants"]):
   axes[0].plot([xs[i] - 0.4, xs[i] - 0.05], [y, y], ls='-', color="gray", lw=2, zorder=-10)

for i, y in enumerate(df_means_bar.loc[df_means_bar["variant class"] == "impactful", "mean number of variants"]):
   axes[0].plot([xs[i] + 0.05, xs[i] + 0.4], [y, y], ls='-', color="gray", lw=2, zorder=-10)

xs = axes[1].get_xticks()
for i, y in enumerate(ad_means_bar.loc[ad_means_bar["variant class"] == "all", "mean number of variants"]):
   axes[1].plot([xs[i] - 0.4, xs[i] - 0.05], [y, y], ls='-', color="gray", lw=2, zorder=-10)

for i, y in enumerate(ad_means_bar.loc[ad_means_bar["variant class"] == "impactful", "mean number of variants"]):
   axes[1].plot([xs[i] + 0.05, xs[i] + 0.4], [y, y], ls='-', color="gray", lw=2, zorder=-10)


print(df_means_bar.loc[df_means_bar["variant class"] == "all", "mean number of variants"])
axes[0].set_ylim(0, 15)

plt.tight_layout()
sns.despine()
plt.savefig("figure3-exome-counts.png")
plt.savefig("figure3-exome-counts.eps")
plt.show()
