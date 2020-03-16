set -euo pipefail
gff=/uufs/chpc.utah.edu/common/HIPAA/u6000771/Data/Homo_sapiens.GRCh38.95.chr.prefix.gff3.gz
fasta=/scratch/ucgd/lustre/common/data/Reference/GRCh38/human_g1k_v38_decoy_phix.fasta
LCR=/uufs/chpc.utah.edu/common/HIPAA/u6000771/Data/LCR-hs38.bed.gz
NOLCR=/uufs/chpc.utah.edu/common/HIPAA/u6000771/Data/hg38.LCR-complement.bed.gz
zip=/uufs/chpc.utah.edu/common/HIPAA/u6000771/Data/gnomad.hg38.genomes.v3.fix.zip

js=~/Projects/src/slivar/paper/paper.js

here="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$here"

mkdir -p vcfs/
cohort=FAKE
d=0
<<DONE

for d in 0 6 8 10 12; do
	for cohort in genome genome-dv; do
	  dn_roc -d $d --gq 1 --gq 5 --gq 10 --gq 15 --gq 20 data-links/$cohort.ped data-links/$cohort.bcf > $cohort-d$d-roc.txt &
	  dn_roc -d $d -x $LCR --gq 1 --gq 5 --gq 10 --gq 15 --gq 20 data-links/$cohort.ped data-links/$cohort.bcf > $cohort-d$d-noLCR-roc.txt &
	done
done
wait

python plot_ab_roc_genome.py "GATK" genome-rocs/genome-d{0,6,8,10,12}-roc.txt
python plot_ab_roc_genome.py "GATK excluding LCR" genome-rocs/genome-d{0,6,8,10,12}-noLCR-roc.txt
python plot_ab_roc_genome.py "DeepVariant excluding LCR" genome-rocs/genome-dv-d{0,6,8,10,12}-noLCR-roc.txt
python plot_ab_roc_genome.py "DeepVariant" genome-rocs/genome-dv-d{0,6,8,10,12}-roc.txt

exit
DONE

region=""
GQ=20
DP=10
out=xx
exclude=xx
<<DENOVO

for exclude in "" "--exclude $LCR"; do
	for cohort in genome genome-dv; do
		
		export SLIVAR_SUMMARY_FILE=y$cohort.dn.summary.tsv
		out=vcfs/$cohort.dn.vcf
		if [[ "$exclude" != "" ]]; then
			export SLIVAR_SUMMARY_FILE=y$cohort.dn.noLCR.summary.tsv
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
			--trio "dn_ab_gq_dp:kid.het && mom.hom_ref && dad.hom_ref && dad.AB < 0.03 && mom.AB < 0.03 && kid.GQ >= $GQ && mom.GQ >= $GQ && dad.GQ >= $GQ && kid.AB >= 0.20 && kid.AB <= 0.80 && kid.DP > 9 && mom.DP > 9 && dad.DP > 9" \
			--trio "dn_gnomad_af_01:kid.het && mom.hom_ref && dad.AB < 0.03 && mom.AB < 0.03 && dad.hom_ref && kid.GQ >= $GQ && mom.GQ >= $GQ && dad.GQ >= $GQ && kid.AB >= 0.20 && kid.AB <= 0.80 && INFO.gnomad_popmax_af < 0.01 && kid.DP > 9 && mom.DP > 9 && dad.DP > 9"  \
			--trio "dn_gnomad_af_001:kid.het && mom.hom_ref && dad.hom_ref && dad.AB < 0.03 && mom.AB < 0.03 && kid.GQ >= $GQ && mom.GQ >= $GQ && dad.GQ >= $GQ && kid.AB >= 0.20 && kid.AB <= 0.80 && INFO.gnomad_popmax_af < 0.001 && kid.DP > 9 && mom.DP > 9 && dad.DP > 9"  \
		--trio "dn_gnomad_af_filter:kid.het && mom.hom_ref && dad.hom_ref && dad.AB < 0.03 && mom.AB < 0.03 && kid.GQ >= $GQ && mom.GQ >= $GQ && dad.GQ >= $GQ && kid.AB >= 0.20 && kid.AB <= 0.80 && INFO.gnomad_popmax_af < 0.001 && !('gnomad_popmax_af_filter' in INFO) && kid.DP >= $DP && mom.DP >= $DP && dad.DP >= $DP"  \
			--trio "dn_impactful:kid.het && mom.hom_ref && dad.hom_ref && dad.AB < 0.03 && mom.AB < 0.03 && kid.GQ >= $GQ && mom.GQ >= $GQ && dad.GQ >= $GQ && kid.AB >= 0.20 && kid.AB <= 0.80 && INFO.gnomad_popmax_af < 0.001 && !('gnomad_popmax_af_filter' in INFO) && INFO.impactful && kid.DP >= $DP && mom.DP >= $DP && dad.DP >= $DP" \
			--trio "dn_impactful_p01:kid.het && mom.hom_ref && dad.hom_ref && dad.AB < 0.03 && mom.AB < 0.03 && kid.GQ >= $GQ && mom.GQ >= $GQ && dad.GQ >= $GQ && kid.AB >= 0.20 && kid.AB <= 0.80 && INFO.gnomad_popmax_af < 0.01 && !('gnomad_popmax_af_filter' in INFO) && INFO.impactful && kid.DP >= $DP && mom.DP >= $DP && dad.DP >= $DP" \
			> $out
#			--trio "dn_gnomad_af_filter_tight_ac:mom.hom_ref && dad.hom_ref && kid.het && dad.AB < 0.03 && mom.AB < 0.03 && kid.GQ >= $GQ && mom.GQ >= $GQ && dad.GQ >= $GQ && kid.AB >= 0.25 && kid.AB <= 0.75 && INFO.gnomad_popmax_af < 0.001 && !('gnomad_popmax_af_filter' in INFO)  && kid.DP >= $DP && mom.DP >= $DP && dad.DP >= $DP && INFO.AC == 1"  &
	done
done

wait
exit


DENOVO


echo "################################################################"
echo "done with denovos"
echo "################################################################"

for cohort in genome genome-dv; do
    for impactful in "INFO.genic" "INFO.impactful" "true"; do
        name=".impactful";
        if [[ "$impactful" == "true" ]]; then name=""; fi
        if [[ "$impactful" == "INFO.genic" ]]; then name=".genic"; fi

        export SLIVAR_SUMMARY_FILE=devx.$cohort$name.summary.tsv
# now get plot of total variant counts:
	        #--region $NOLCR \

        pslivar expr --vcf data-links/$cohort.bcf --ped data-links/$cohort.ped \
			--fasta $fasta \
			-g $zip \
			--exclude /uufs/chpc.utah.edu/common/HIPAA/u6000771/Data/LCR-hs38.bed.gz \
            --pass-only \
            --js $js \
            --info "$impactful && INFO.gnomad_popmax_af < 0.01 && variant.FILTER == 'PASS' && variant.ALT[0] != '*' && !('gnomad_popmax_af_filter' in INFO)" \
            --trio "denovo:denovo(kid, mom, dad) && INFO.gnomad_popmax_af < 0.001" \
            --trio "recessive:recessive(kid, mom, dad)" \
            --trio "x_denovo:x_denovo(kid, mom, dad) && (variant.CHROM == 'chrX' || variant.CHROM == 'X')" \
            --trio "x_recessive:x_recessive(kid, mom, dad) && (variant.CHROM == 'chrX' || variant.CHROM == 'X')" \
            --trio "auto_dom:fake_auto_dom(kid, mom, dad) && variant.CHROM != 'chrX' && variant.CHROM != 'X' && INFO.gnomad_popmax_af < 0.001 && INFO.gnomad_nhomalt < 10" \
            --trio "comphet_side:comphet_side(kid, mom, dad)" \
            > vcfs/$cohort$name.vcf
# from compount het
# && INFO.gnomad_nhomalt_controls < 10 && INFO.gnomad_popmax_af < 0.01
        export SLIVAR_SUMMARY_FILE=devx.$cohort$name.ch.summary.tsv
        slivar compound-hets --sample-field comphet_side --sample-field denovo -p data-links/$cohort.ped -v vcfs/$cohort$name.vcf > vcfs/$cohort$name.ch.vcf
    done
done

#python plot-final-exome.py exome.summary.tsv exome.ch.summary.tsv
