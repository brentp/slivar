import pandas as pd
import seaborn as sns

import matplotlib
matplotlib.use("Agg")
sns.set_style("whitegrid")

from matplotlib import rc
#from brokenaxes import brokenaxes

rc('font',**{'family':'sans-serif','sans-serif':['Arial']})


from matplotlib import pyplot as plt

import sys

colors = sns.color_palette("Set1", 12)

c = 0.5
colors = [(c, c, c), (c, c, c)]

args = sys.argv[1:][:]
#assert len(sys.argv) == 5

# order is [GATK, DV, GATK-noLCR, DV-noLCR]
#args[0] = [x for x in sys.argv[1:] if not "noLCR" in x and not 'dv' in x][0]
#args[1] = [x for x in sys.argv[1:] if not "noLCR" in x and 'dv' in x][0]

# sample  dn_ab_gq_dp     dn_gnomad_af_01 dn_gnomad_af_001        dn_gnomad_af_001_parent_AB      dn_impactful    dn_impactful_p01

dfoxs = [pd.read_csv(arg, sep="\t", index_col=0) for arg in args]
dfos = [df["dn_ab_gq_dp dn_gnomad_af_01	dn_gnomad_af_001 dn_gnomad_af_001_parent_AB".split()] for df in dfoxs]

#dfimp = [df[["dn_impactful"]] for df in dfoxs]

for d in dfos:
    d["sample"] = d.index

#for d in dfimp:
#    d["sample"] = d.index


dfs = [dfo.melt(value_name="number_of_variants", id_vars="sample") for dfo in dfos]

#dfimps = [dfi.melt(value_name="number_of_variants", id_vars="sample") for dfi in dfimp]
#dfimps = [dfi.loc[~dfi["sample"].isin(("RGP_426_3", "RGP_430_3")), :] for dfi in dfimps]


dfs = [d.loc[~d["sample"].isin(("RGP_426_3", "RGP_430_3")), :] for d in dfs]

#print("GATK all impactful:", dfimps[0]["number_of_variants"].median(), dfimps[0]["number_of_variants"].mean())
#print("DeepVariant all impactful:", dfimps[1]["number_of_variants"].median(), dfimps[1]["number_of_variants"].mean())
#print("GATK no LCR impactful:", dfimps[2]["number_of_variants"].median(), dfimps[2]["number_of_variants"].mean())
#print("DeepVariant no LCR impactful:", dfimps[3]["number_of_variants"].median(), dfimps[3]["number_of_variants"].mean())


dfs = [[dfs[0], dfs[1]]]#, [dfs[2], dfs[3]]]


fig, axes = plt.subplots(1, 2, figsize=(7, 4), sharey=True)


for i in range(0, 2):
    df = dfs[0][i]
    sns.boxplot(x="variable", y="number_of_variants",
             data=df, ax=axes[i], palette=colors) #, s=3)

axes[0].set_title("GATK", size=15)
axes[1].set_title("DeepVariant", size=15)

axes[0].set_xlabel(None, size=15)
axes[1].set_xlabel(None, size=15)

axes[0].set_ylabel("Candidate $\it{de novo}$ variants")
axes[1].set_ylabel(None)

axes[0].set_xlabel(None, size=15)
axes[1].set_xlabel(None, size=15)

axes[0].set_ylabel("Candidate $\it{de}$ $\it{novo}$ variants")
axes[1].set_ylabel(None)

xlabels = ["AB in 0.2..0.8\nGQ >= 20\nDP >= 10",
           "AF < 0.01",
           "AF < 0.001",
           "parent AB < 0.02"]

lbls = [x for x in axes[0].get_xticklabels()]
for i, l in enumerate(lbls):
    print(i)
    l.set_text(xlabels[i])
axes[0].set_xticklabels(lbls, rotation=12, ha="center")
axes[1].set_xticklabels(lbls, rotation=12, ha="center")

ax = axes[0]
ax.set_ylim(ymin=-0.5)
#ax.axhline(o, 0.1, 0.4, lw=4)
#ax.axhline(f, 1-0.4, 1 - 0.1, lw=4)

for ax in axes:
    ax.annotate("additive filtering", [3, 1600], [1.4, 1600],
        ha="center", va="center",
        arrowprops=dict(arrowstyle="simple"))


sns.despine()
plt.tight_layout()
plt.savefig("figure4-genome-denovos.png")
plt.savefig("figure4-genome-denovos.eps")
plt.show()
