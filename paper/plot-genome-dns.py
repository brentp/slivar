import pandas as pd
import seaborn as sns

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
assert len(sys.argv) == 5

# order is [GATK, DV, GATK-noLCR, DV-noLCR]
args[0] = [x for x in sys.argv[1:] if not "noLCR" in x and not 'dv' in x][0]
args[1] = [x for x in sys.argv[1:] if not "noLCR" in x and 'dv' in x][0]
args[2] = [x for x in sys.argv[1:] if "noLCR" in x and not 'dv' in x][0]
args[3] = [x for x in sys.argv[1:] if "noLCR" in x and 'dv' in x][0]

dfoxs = [pd.read_csv(arg, sep="\t", index_col=0) for arg in args]
dfos = [df["dn_gnomad_af_01	dn_gnomad_af_001 dn_gnomad_af_filter".split()] for df in dfoxs]

dfimp = [df[["dn_impactful"]] for df in dfoxs]

for d in dfos:
    d["sample"] = d.index

for d in dfimp:
    d["sample"] = d.index


dfs = [dfo.melt(value_name="number_of_variants", id_vars="sample") for dfo in dfos]

dfimps = [dfi.melt(value_name="number_of_variants", id_vars="sample") for dfi in dfimp]

dfimps = [dfi.loc[~dfi["sample"].isin(("RGP_426_3", "RGP_430_3")), :] for dfi in dfimps]


dfs = [d.loc[~d["sample"].isin(("RGP_426_3", "RGP_430_3")), :] for d in dfs]

print("GATK all impactful:", dfimps[0]["number_of_variants"].median(), dfimps[0]["number_of_variants"].mean())
print("DeepVariant all impactful:", dfimps[1]["number_of_variants"].median(), dfimps[1]["number_of_variants"].mean())
print("GATK no LCR impactful:", dfimps[2]["number_of_variants"].median(), dfimps[2]["number_of_variants"].mean())
print("DeepVariant no LCR impactful:", dfimps[3]["number_of_variants"].median(), dfimps[3]["number_of_variants"].mean())


dfs = [[dfs[0], dfs[1]], [dfs[2], dfs[3]]]


fig, axes = plt.subplots(2, 2, figsize=(7, 7), sharey=True)


for i in range(0, 2):
    for j in range(0, 2):
        df = dfs[i][j]

        sns.boxplot(x="variable", y="number_of_variants",
               data=df, ax=axes[i, j], palette=colors) #, s=3)


axes[0,0].set_title("GATK", size=15)
axes[0,1].set_title("DeepVariant", size=15)
axes[1,0].set_title("GATK excluding LCR", size=15)
axes[1,1].set_title("DeepVariant excluding LCR", size=15)

axes[0,0].set_xlabel(None, size=15)
axes[0,1].set_xlabel(None, size=15)
axes[1,0].set_xlabel(None, size=15)
axes[1,1].set_xlabel(None, size=15)

axes[0,0].set_ylabel("Candidate $\it{de novo}$ variants")
axes[1,0].set_ylabel("Candidate $\it{de novo}$ variants")
axes[0,1].set_ylabel(None)
axes[1,1].set_ylabel(None)

xlabels = ["AB in 0.2..0.8\nGQ >= {cutoff}\nDP >= 10\nAF < 0.01",
           "AF < 0.001",
           "not FILTER'ed\nin gnomAD"]

axes[0, 0].set_xticklabels([])


axes[1, 0].set_xticklabels([x.format(cutoff=20) for x in xlabels], rotation=12, ha="center")
axes[1, 1].set_xticklabels([x.format(cutoff=15) for x in xlabels], rotation=12, ha="center")
axes[0, 0].set_xticklabels([])
axes[0, 1].set_xticklabels([])



#axes[0,0].set_xticklabels(["0.2 <= AB < 0.8\n& GQ >= 20",
#                         "gnomAD popmax AF < 0.001\n& not non-PASS in gnomAD",
#                         "impactful"])

#axes[0].set_yticks(range(11))

#print(axes[0].get_xlim())

ax = axes[0, 0]
ax.set_ylim(-0.5, 1350)
#ax.axhline(o, 0.1, 0.4, lw=4)
#ax.axhline(f, 1-0.4, 1 - 0.1, lw=4)

sns.despine()
plt.tight_layout()
plt.savefig("figure4-genome-denovos.png")
plt.show()




