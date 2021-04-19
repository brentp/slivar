import sys
from cyvcf2 import VCF
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_curve

np.random.seed(42)

colors = sns.color_palette()
dpi = 120

# python ab-truth.py eval/fp.vcf.gz eval/tp.vcf.gz eval/fn.vcf.gz

fp_path = sys.argv[1]
tp_path = sys.argv[2]
fn_path = sys.argv[3]


def read(path, tpi):
    vcf = VCF(path)

    for v in vcf:
        if v.FILTER is not None: continue
        if len(v.ALT) > 1: continue

        gs = v.genotypes[0]
        if gs != [0, 1, False]: continue

        try:
            ad = v.format("AD")[0]
            ab = ad[1] / float(ad[0] + ad[1])
        except:
            ab = 0.0

        if ab > 0.5: ab = 1 - ab

        try:
            DP = v.format("DP")[0][0]
        #if gq > 40: continue
            GQ = float(v.format("GQ")[0][0])
        except:
            DP = 0
            GQ = 0
        yield (tpi, GQ, ab, int(DP))



fps = list(read(fp_path, 0))
tps = list(read(tp_path, 1))

ps = np.array(fps + tps)

fns = np.array(list(read(fn_path, 1)))

fig, ax = plt.subplots(3, 1, figsize=(4, 7), dpi=dpi)

trs = []
tprs = []
fprs = []

for i, n in enumerate(["GQ", "AB", "DP"], start=1):
    ps[:, i] += (np.random.random(ps.shape[0]) - 0.5) / 151.0
    fpr, tpr, tr = roc_curve(ps[:, 0], ps[:, i], drop_intermediate=False)
    print(tr.shape, fpr.shape, tpr.shape)
    fprs.append(fpr); tprs.append(tpr); trs.append(np.array(tr))

    ax[i-1].plot(fpr, tpr, color=colors[0], label="all %s values" % n)
    if n == "GQ":
        x = fpr[tr >= 20][-1]
        y = tpr[tr >= 20][-1]
        for j, cutoff in enumerate((20, 30, 40)):
            ax[i-1].plot(fpr[tr >= cutoff], tpr[tr >= cutoff], label="GQ >= %d" %
                    cutoff, c=colors[j+1])

    if n == "AB":
        x = fpr[tr >= 0.2][-1]
        y = tpr[tr >= 0.2][-1]
        for j, cutoff in enumerate((0.1, 0.2, 0.3)):
            ax[i-1].plot(fpr[tr >= cutoff], tpr[tr >= cutoff], label="AB >= %.2f" %
                    cutoff, c=colors[j+1])
    if n == "DP":
        x = fpr[tr >= 10][-1]
        y = tpr[tr >= 10][-1]
        for j, cutoff in enumerate((10, 20, 30)):
            ax[i-1].plot(fpr[tr >= cutoff], tpr[tr >= cutoff], label="DP >= %d" %
                    cutoff, c=colors[j+1])

    ax[i-1].plot([x], [y], 'o', c=colors[0], label="selected %s cutoff\nTPR: %.3f FPR: %.3f" % (n, y, x))


    if i == len(ax):
        ax[i-1].set_xlabel("FPR")
    ax[i-1].set_ylabel("TPR")
    ax[i-1].text(0.95, 0.05, "varying %s" % n, transform=ax[i-1].transAxes,
            ha="right", size=17)
    #plt.plot(fpr[tr >= 40], tpr[tr >= 40], label="GQ >= 40")

print("fns:", fns.shape)

#plt.suptitle("ROC Curve for DeepVariant on GIABv4.1", size=15)
# now show combined TPR, FPR
tprs = np.array(tprs)
fprs = np.array(fprs)
print(trs[0].shape)
trs = np.array(trs)
print("trs.shape:", trs.shape)
print(trs)

sel = (ps[:, 1] >= 20) & (ps[:, 2] >= 0.2) & (ps[:, 3] >= 10)

m = ps[:, 3].mean()
s = ps[:, 3].std()
#sel = (ps[:, 1] >= 20) & (ps[:, 2] >= 0.2) & (ps[:, 3] >= 8) # & (ps[:, 3] <= (m + 3*s))

from sklearn.metrics import confusion_matrix
from sklearn.metrics import plot_confusion_matrix
from sklearn.metrics import ConfusionMatrixDisplay

print(ps.shape, sel.shape)

cm = confusion_matrix(ps[:, 0].astype(int), sel.astype(int))
tn, fp, fn, tp = cm.ravel()

tot_fn = fn + len(fns)

filtered_fdr = fp / (fp + tp)
unfiltered_fdr = len(fps) / (len(fps) + len(tps))
print(f"filtered    tn:{tn} fp:{fp} fn:{tot_fn} tp:{tp}  TPR:{tp / (tot_fn + tp):.3f} FDR:{filtered_fdr:.5f}")
print(f"unfiltered          fp:{len(fps)} fn:{len(fns)} tp:{len(tps)} TPR: {len(tps) / (len(tps) + len(fns)):.3f} FDR:{unfiltered_fdr:.5f}")
print(f"improvement in fdr:{unfiltered_fdr/filtered_fdr:.2f}X")

plt.tight_layout() #h_pad=4, pad=4)
plt.savefig("supp-figure-12.png", dpi=dpi)
plt.savefig("supp-figure-12.eps", dpi=dpi)
plt.show()
