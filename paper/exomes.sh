set -euo pipefail
gff=/data/human/Homo_sapiens.GRCh38.96.chr.gff3.gz
bcf=/home/brentp/src/varx/var.bcf
ped=/home/brentp/src/varx/var.ped
fasta=/data/human/GRCh38_full_analysis_set_plus_decoy_hla.fa
LCR=/data/human/LCR-hs38.bed.gz
mkdir -p vcfs/

cohort=exome
export SLIVAR_SUMMARY_FILE=$cohort.summary.tsv PATH=.:$PATH

here="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$here"
mkdir -p vcfs

<<DONE
./dn_roc --gq 1 --gq 5 --gq 10 --gq 20 $ped $bcf > $cohort-roc.txt
python plot_ab_roc.py $cohort-roc.txt

export SLIVAR_SUMMARY_FILE=$cohort.dn.summary.tsv

../slivar expr --vcf $bcf --ped $ped \
    --pass-only \
    -o vcfs/$cohort.dn.bcf \
    -g /home/brentp/src/slivar/gnomad.hg38.genomes.v3.fix.zip \
    --info "variant.FILTER == 'PASS' && variant.ALT[0] != '*'" \
    --trio "dn_ab_gq:mom.hom_ref && dad.hom_ref && kid.het && kid.GQ >= 20 && mom.GQ >= 20 && dad.GQ >= 20 && kid.AB >= 0.2 && kid.AB <= 0.8" \
    --trio "dn:mom.hom_ref && dad.hom_ref && kid.het && kid.GQ >= 20 && mom.GQ >= 20 && dad.GQ >= 20 && kid.AB >= 0.2 && kid.AB <= 0.8 && INFO.gnomad_popmax_af < 0.001"  \
    --trio "impactful_dn:INFO.impactful && mom.hom_ref && dad.hom_ref && kid.het && kid.GQ >= 20 && mom.GQ >= 20 && dad.GQ >= 20 && kid.AB >= 0.2 && kid.AB <= 0.8 && INFO.gnomad_popmax_af < 0.001"  \

DONE
python plot-exome-dn-summary.py $cohort.dn.summary.tsv


DONE

export SLIVAR_SUMMARY_FILE=$cohort.summary.tsv
echo $(which slivar)

pslivar expr --vcf $bcf \
    --fasta $fasta \
    --ped $ped \
    --pass-only \
    --js /home/brentp/src/slivar/js/slivar-functions.js \
    -g /home/brentp/src/slivar/gnomad.hg38.genomes.v3.fix.zip \
    --info 'INFO.gnomad_popmax_af < 0.01 && variant.FILTER == "PASS" && variant.ALT[0] != "*"' \
    --family-expr 'denovo:fam.every(segregating_denovo) && INFO.gnomad_popmax_af < 0.001' \
    --family-expr 'recessive:fam.every(segregating_recessive)' \
    --family-expr 'x_denovo:fam.every(segregating_denovo_x) && INFO.gnomad_popmax_af < 0.001 && variant.CHROM == "chrX"' \
    --family-expr 'x_recessive:fam.every(segregating_recessive_x) && variant.CHROM == "chrX"' \
    --trio "auto_dom:fake_auto_dom(kid, mom, dad) && variant.CHROM != 'chrX' && variant.CHROM != 'X' && INFO.gnomad_popmax_af < 0.001 && INFO.gnomad_nhomalt < 10" \
    --trio 'comphet_side:comphet_side(kid, mom, dad) && INFO.gnomad_nhomalt < 10 && INFO.gnomad_popmax_af < 0.005' \
    | bcftools csq -s - --ncsq 40 -g $gff -l -f $fasta - -o vcfs/$cohort.vcf

export SLIVAR_SUMMARY_FILE=$cohort.ch.summary.tsv
slivar compound-hets --sample-field comphet_side --sample-field denovo -p $ped -v vcfs/$cohort.vcf > vcfs/$cohort.ch.vcf

export SLIVAR_SUMMARY_FILE=$cohort.impactful.summary.tsv
pslivar expr --vcf $bcf \
    --fasta $fasta \
    --ped $ped \
    --pass-only \
    --js /home/brentp/src/slivar/js/slivar-functions.js \
    -g /home/brentp/src/slivar/gnomad.hg38.genomes.v3.fix.zip \
    --info 'INFO.impactful && INFO.gnomad_popmax_af < 0.01 && variant.FILTER == "PASS" && variant.ALT[0] != "*"' \
    --family-expr 'denovo:fam.every(segregating_denovo) && INFO.gnomad_popmax_af < 0.001' \
    --family-expr 'recessive:fam.every(segregating_recessive)' \
    --family-expr 'x_denovo:fam.every(segregating_denovo_x) && INFO.gnomad_popmax_af < 0.001 && variant.CHROM == "chrX"' \
    --family-expr 'x_recessive:fam.every(segregating_recessive_x) && variant.CHROM == "chrX"' \
    --trio "auto_dom:fake_auto_dom(kid, mom, dad) && variant.CHROM != 'chrX' && variant.CHROM != 'X' && INFO.gnomad_popmax_af < 0.001 && INFO.gnomad_nhomalt < 10" \
    --trio 'comphet_side:comphet_side(kid, mom, dad) && INFO.gnomad_nhomalt < 10 && INFO.gnomad_popmax_af < 0.005' \
    | bcftools csq -s - --ncsq 40 -g $gff -l -f $fasta - -o vcfs/$cohort.impactful.vcf

export SLIVAR_SUMMARY_FILE=$cohort.impactful.ch.summary.tsv
slivar compound-hets --sample-field comphet_side --sample-field denovo -p $ped -v vcfs/$cohort.impactful.vcf > vcfs/$cohort.impactful.ch.vcf

python plot-final-exome.py $cohort.summary.tsv $cohort.ch.summary.tsv $cohort.impactful.summary.tsv $cohort.impactful.ch.summary.tsv $ped
