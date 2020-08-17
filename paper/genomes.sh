set -euo pipefail
gff=/uufs/chpc.utah.edu/common/HIPAA/u6000771/Data/Homo_sapiens.GRCh38.95.chr.prefix.gff3.gz
fasta=/scratch/ucgd/lustre/common/data/Reference/GRCh38/human_g1k_v38_decoy_phix.fasta
LCR=/uufs/chpc.utah.edu/common/HIPAA/u6000771/Data/LCR-hs38.bed.gz
zip=/uufs/chpc.utah.edu/common/HIPAA/u6000771/Data/gnomad.hg38.genomes.v3.fix.zip

SRC=~/Projects/src/slivar/paper/
js=~/Projects/src/slivar/js/slivar-functions.js

here="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$here"

mkdir -p vcfs/
cohort=FAKE
d=0
<<DONE

for d in 0 6 8 10 12; do
	for cohort in genome genome-dv; do
	  dn_roc -d $d --gq 1 --gq 5 --gq 10 --gq 15 --gq 20 data-links/$cohort.ped data-links/$cohort.bcf > genome-rocs/$cohort-d$d-roc.txt &
	  dn_roc -d $d -x $LCR --gq 1 --gq 5 --gq 10 --gq 15 --gq 20 data-links/$cohort.ped data-links/$cohort.bcf > genome-rocs/$cohort-d$d-noLCR-roc.txt &
	done
done
wait

python $SRC/plot_ab_roc_genome.py "GATK" genome-rocs/genome-d{0,6,8,10,12}-roc.txt
python $SRC/plot_ab_roc_genome.py "GATK excluding LCR" genome-rocs/genome-d{0,6,8,10,12}-noLCR-roc.txt
python $SRC/plot_ab_roc_genome.py "DeepVariant excluding LCR" genome-rocs/genome-dv-d{0,6,8,10,12}-noLCR-roc.txt
python $SRC/plot_ab_roc_genome.py "DeepVariant" genome-rocs/genome-dv-d{0,6,8,10,12}-roc.txt

exit
DONE

GQ=20
DP=10
parent_AB=0.02
out=xx
exclude=xx

set +u
<<DENOVO

for exclude in "" "--exclude $LCR"; do
	for cohort in genome genome-dv; do
		
		export SLIVAR_SUMMARY_FILE=pab$cohort.dn.summary.tsv
		out=vcfs/$cohort.dn.vcf
		if [[ "$exclude" != "" ]]; then
			export SLIVAR_SUMMARY_FILE=pab$cohort.dn.noLCR.summary.tsv
		    out=vcfs/$cohort.dn.noLCR.vcf
		fi;

		#--info "variant.FILTER == 'PASS' && variant.ALT[0] != '*' & INFO.VQSLOD > 0.0 && variant.CHROM != 'chrX' && variant.CHROM != 'X' && variant.ALT[0].length < 5 && variant.REF.length < 9" \
		pslivar expr --vcf data-links/$cohort.bcf --ped data-links/$cohort.ped \
			$exclude \
			--fasta $fasta \
			-g $zip \
			-g /uufs/chpc.utah.edu/common/HIPAA/u6000771/Data/1kg-dv/slivar.thousand_genomes.deep_variant.hg38.zip \
			--pass-only \
			--info "variant.FILTER == 'PASS' && variant.ALT[0] != '*' && variant.CHROM != 'chrX' && variant.CHROM != 'X' && variant.CHROM != 'chrY'" \
			--trio "dn_ab_gq_dp:kid.het && mom.hom_ref && dad.hom_ref && kid.GQ >= $GQ && mom.GQ >= $GQ && dad.GQ >= $GQ && kid.AB >= 0.20 && kid.AB <= 0.80 && kid.DP >= $DP && mom.DP >= $DP && dad.DP >= $DP" \
			--trio "dn_gnomad_af_01:kid.het && mom.hom_ref && dad.hom_ref && kid.GQ >= $GQ && mom.GQ >= $GQ && dad.GQ >= $GQ && kid.AB >= 0.20 && kid.AB <= 0.80 && INFO.gnomad_popmax_af < 0.01 && kid.DP >= $DP && mom.DP >= $DP && dad.DP >= $DP"  \
			--trio "dn_gnomad_af_001:kid.het && mom.hom_ref && dad.hom_ref && dad.AB < $parent_AB && mom.AB < $parent_AB && kid.GQ >= $GQ && mom.GQ >= $GQ && dad.GQ >= $GQ && kid.AB >= 0.20 && kid.AB <= 0.80 && INFO.gnomad_popmax_af < 0.001"  \
			--trio "dn_gnomad_af_001_parent_AB:kid.het && mom.hom_ref && dad.hom_ref && dad.AB < $parent_AB && mom.AB < $parent_AB && kid.GQ >= $GQ && mom.GQ >= $GQ && dad.GQ >= $GQ && kid.AB >= 0.20 && kid.AB <= 0.80 && INFO.gnomad_popmax_af < 0.001 && kid.DP >= $DP && mom.DP >= $DP && dad.DP >= $DP"  \
			--trio "dn_impactful:kid.het && mom.hom_ref && dad.hom_ref && dad.AB < $parent_AB && mom.AB < $parent_AB && kid.GQ >= $GQ && mom.GQ >= $GQ && dad.GQ >= $GQ && kid.AB >= 0.20 && kid.AB <= 0.80 && INFO.gnomad_popmax_af < 0.001 && INFO.impactful && kid.DP >= $DP && mom.DP >= $DP && dad.DP >= $DP" \
			--trio "dn_impactful_p01:kid.het && mom.hom_ref && dad.hom_ref && dad.AB < $parent_AB && mom.AB < $parent_AB && kid.GQ >= $GQ && mom.GQ >= $GQ && dad.GQ >= $GQ && kid.AB >= 0.20 && kid.AB <= 0.80 && INFO.gnomad_popmax_af < 0.01 && INFO.impactful && kid.DP >= $DP && mom.DP >= $DP && dad.DP >= $DP" \
			> $out
	done
done

wait
exit

python ~/Projects/src/slivar/paper/plot-genome-dns.py pabgenome.dn.summary.tsv pabgenome-dv.dn.summary.tsv
mv figure4-genome-denovos.png supp-figure-genome-denovos-including-LCR.png
python ~/Projects/src/slivar/paper/plot-genome-dns.py pabgenome.dn.noLCR.summary.tsv pabgenome-dv.dn.noLCR.summary.tsv

# total denovos:
#awk 'NF > 1 && $1 != "RGP_426_3" && $1 != "RGP_430_3" {s += $5; t += 1 } END{ print s/t }' pabgenome.dn.noLCR.summary.ts

# impactful 0.01
#awk 'NF > 1 && $1 != "RGP_426_3" && $1 != "RGP_430_3" {s += $7; t += 1 } END{ print s/t }' pabgenome.dn.noLCR.summary.tsv


DENOVO
set -u


echo "################################################################"
echo "done with denovos"
echo "################################################################"

for cohort in genome genome-dv; do
    for impactful in "INFO.genic" "INFO.impactful" "true"; do
        name="";
        if [[ "$impactful" == "INFO.impactful" ]]; then name=".impactful"; fi
        if [[ "$impactful" == "INFO.genic" ]]; then name=".genic"; fi

        export SLIVAR_SUMMARY_FILE=dev.$cohort$name.summary.tsv
# now get plot of total variant counts:

        pslivar expr --vcf data-links/$cohort.bcf --ped data-links/$cohort.ped \
			--fasta $fasta \
			-g $zip \
			--exclude /uufs/chpc.utah.edu/common/HIPAA/u6000771/Data/LCR-hs38.bed.gz \
            --pass-only \
            --js $js \
            --info "$impactful && INFO.gnomad_popmax_af < 0.01 && variant.FILTER == 'PASS' && variant.ALT[0] != '*'" \
            --trio "denovo:denovo(kid, mom, dad) && INFO.gnomad_popmax_af < 0.001" \
            --trio "recessive:recessive(kid, mom, dad)" \
            --trio "x_denovo:x_denovo(kid, mom, dad) && (variant.CHROM == 'chrX' || variant.CHROM == 'X')" \
            --trio "x_recessive:x_recessive(kid, mom, dad) && (variant.CHROM == 'chrX' || variant.CHROM == 'X')" \
            --trio "auto_dom:fake_auto_dom(kid, mom, dad) && variant.CHROM != 'chrX' && variant.CHROM != 'X' && INFO.gnomad_popmax_af < 0.001 && INFO.gnomad_nhomalt < 10" \
            --trio "comphet_side:comphet_side(kid, mom, dad)" \
            > vcfs/$cohort$name.vcf
# from compount het
# && INFO.gnomad_nhomalt_controls < 10 && INFO.gnomad_popmax_af < 0.01
        export SLIVAR_SUMMARY_FILE=dev.$cohort$name.ch.summary.tsv
        slivar compound-hets --sample-field comphet_side --sample-field denovo -p data-links/$cohort.ped -v vcfs/$cohort$name.vcf > vcfs/$cohort$name.ch.vcf
    done
done

# figure 5a
python $SRC/plot-final-genome.py dev.genome.impactful.summary.tsv dev.genome.impactful.ch.summary.tsv dev.genome-dv.impactful.summary.tsv dev.genome-dv.impactful.ch.summary.tsv
mv figure5-genome-counts.png figure5-genome-counts.impactful.png
mv figure5-genome-counts.eps figure5-genome-counts.impactful.eps
# figure 5b
python $SRC/plot-final-genome.py dev.genome.genic.summary.tsv dev.genome.genic.ch.summary.tsv dev.genome-dv.genic.summary.tsv dev.genome-dv.genic.ch.summary.tsv
mv figure5-genome-counts.png figure5-genome-counts.genic.png
mv figure5-genome-counts.eps figure5-genome-counts.genic.eps


## for checking that we don't remove "true" positives:
slivar tsv --sample-field denovo --sample-field recessive --sample-field x_denovo --sample-field x_recessive -o genomes.tsv -p data-links/genome.ped --csq-column CSQ vcfs/genome.vcf
slivar tsv --sample-field slivar_comphet -p data-links/genome.ped --csq-column CSQ vcfs/genome.ch.vcf | grep -v ^# >> genomes.tsv

# solved cases are from: https://app.terra.bio/#workspaces/broad-genomics-delivery/RGP_Rehm_RareDisease_WGS/data downloaed 8/15/2020

