import pandas as pd
import numpy as np
import matplotlib
#matplotlib.use("Agg")
import seaborn as sns
from matplotlib import pyplot as plt

import sys

colors2 = sns.color_palette("Set1", 12)[6:]

colors = [(0.6, 0.6, 0.6)] * 12


def read_dfs(i, j):

    df = pd.read_csv(sys.argv[i], sep="\t", index_col="sample")
    df2 = pd.read_csv(sys.argv[j], sep="\t", index_col="sample")
    try:
        df = df.drop(["RGP_426_3", "RGP_430_3", "RGP_462_3"])
        df2 = df2.drop(["RGP_426_3", "RGP_430_3", "RGP_462_3"])
    except:
        pass

    dfo = df.join(df2)

    dfo.drop("comphet_side", inplace=True, axis="columns")
    ad = dfo.auto_dom.copy()
    ad = pd.DataFrame({"number_of_variants": ad})
    dfo.drop("auto_dom", inplace=True, axis="columns")
    df = dfo.melt(value_name="number_of_variants")
    df["number_of_variants_r"] = df.number_of_variants + np.random.random(size=df.shape[0]) / 2.0
    ad["number_of_variants_r"] = ad.number_of_variants + np.random.random(size=ad.shape[0]) / 2.0
    return df, ad

def get_mean_df(df_means, order, use_rand):
    field = "number_of_variants_r" if use_rand else "number_of_variants"

    df_means = df_means.to_dict(orient="dict")[field]

    df_means_o = [df_means[(l, "GATK")] for l in order]
    df_means_i = [df_means[(l, "DeepVariant")] for l in order]

    df_means_b = [[labels[i], o, "GATK"] for i, o in enumerate(df_means_o)] + \
                 [[labels[i], o, "DeepVariant"] for i, o in enumerate(df_means_i)]

    df_means_bar = pd.DataFrame(df_means_b, columns=["inheritance mode", "mean number of variants",
        "caller"])
    return df_means_bar


df, ad = read_dfs(1, 2)

dfi, adi = read_dfs(3, 4)

df["cat"] = "GATK"
dfi["cat"] = "DeepVariant"

ad["cat"] = "GATK"
adi["cat"] = "DeepVariant"

df = pd.concat((df, dfi))
ad = pd.concat((ad, adi), ignore_index=True)
ad["variable"] = "autosomal dominant"


df_means = df.groupby(["variable", "cat"]).mean()
ad_means = ad.groupby(["variable", "cat"]).mean()
#print(df_means, ad_means)
print(ad_means)


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
axes[0].get_legend().set_title("Caller")


#print(df_means.head())
df_means_bar = get_mean_df(df_means, order, True)
ad_means_bar = get_mean_df(ad_means, ["autosomal dominant"], False)

xs = axes[0].get_xticks()
for i, y in enumerate(df_means_bar.loc[df_means_bar["caller"] == "GATK", "mean number of variants"]):
   axes[0].plot([xs[i] - 0.4, xs[i] - 0.05], [y, y], ls='-', color="gray", lw=2, zorder=-10)

for i, y in enumerate(df_means_bar.loc[df_means_bar["caller"] == "DeepVariant", "mean number of variants"]):
   axes[0].plot([xs[i] + 0.05, xs[i] + 0.4], [y, y], ls='-', color="gray", lw=2, zorder=-10)

xs = axes[1].get_xticks()
for i, y in enumerate(ad_means_bar.loc[ad_means_bar["caller"] == "GATK", "mean number of variants"]):
   axes[1].plot([xs[i] - 0.4, xs[i] - 0.05], [y, y], ls='-', color="gray", lw=2, zorder=-10)

for i, y in enumerate(ad_means_bar.loc[ad_means_bar["caller"] == "DeepVariant", "mean number of variants"]):
   axes[1].plot([xs[i] + 0.05, xs[i] + 0.4], [y, y], ls='-', color="gray", lw=2, zorder=-10)


axes[0].set_ylim(0, 15)
#print(np.arange(len(labels)))

#print(df_means_bar.loc[df_means_bar["caller"] == "all", "mean number of variants"])
#print(axes[0].get_xticks())

df_means_bar = get_mean_df(df_means, order, False)
print(df_means_bar)
print("AD", ad_means_bar)

plt.suptitle("Impactful Candidate Variants")
#axes[1].set_title("Genic Candidate Variants")

plt.tight_layout()
sns.despine()
plt.savefig("figure5-genome-counts.png")
plt.savefig("figure5-genome-counts.eps")
plt.show()
