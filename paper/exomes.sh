set -euo pipefail
gff=/data/human/Homo_sapiens.GRCh38.96.chr.gff3.gz
bcf=/home/brentp/src/varx/var.bcf
ped=/home/brentp/src/varx/var.ped
fasta=/data/human/GRCh38_full_analysis_set_plus_decoy_hla.fa
LCR=/home/brentp/src/slivar/LCR-hs38.bed.gz
mkdir -p vcfs/

cohort=exome
export SLIVAR_SUMMARY_FILE=$cohort.summary.tsv
PATH=.:$PATH

here="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$here"
mkdir -p vcfs

<<DONE
# get all passing variants and annotate with gnomad and topmed
slivar expr --vcf $bcf \
    --pass-only \
    -g /home/brentp/src/slivar/gnomad.hg38.v2.zip \
    -g /home/brentp/src/slivar/topmed.hg38.dbsnp.151.zip \
    -o vcfs/$cohort.annotated.bcf \
    --info "variant.FILTER == 'PASS'" 
DONE

./dn_roc --gq 1 --gq 5 --gq 10 --gq 15 $ped vcfs/$cohort.annotated.bcf > $cohort-roc.txt
python plot_ab_roc.py $cohort-roc.txt

exit

bcf=vcfs/$cohort.annotated.bcf

#<<DENOVO
export SLIVAR_SUMMARY_FILE=$cohort.dn.summary.tsv

slivar expr --vcf $bcf --ped $ped \
    --pass-only \
    -o vcfs/$cohort.dn.bcf \
    --info "variant.FILTER == 'PASS' && variant.ALT[0] != '*'" \
    --trio "dn_ab_gq:mom.hom_ref && dad.hom_ref && kid.het && kid.GQ >= 10 && mom.GQ >= 10 && dad.GQ >= 10 && kid.AB >= 0.2 && kid.AB <= 0.8" \
    --trio "dn:mom.hom_ref && dad.hom_ref && kid.het && kid.GQ >= 10 && mom.GQ >= 10 && dad.GQ >= 10 && kid.AB >= 0.2 && kid.AB <= 0.8 && INFO.gnomad_popmax_af < 0.001 && !('gnomad_popmax_af_filter' in INFO) && INFO.topmed_af < 0.01"  \

python plot-exome-dn-summary.py $cohort.dn.summary.tsv

exit
DENOVO

export SLIVAR_SUMMARY_FILE=$cohort.summary.tsv
# now get plot of total variant counts:
slivar expr --vcf $bcf --ped $ped \
    --pass-only \
    --js paper.js \
    --info "variant.FILTER == 'PASS' && variant.ALT[0] != '*' && INFO.gnomad_popmax_af < 0.01 && INFO.topmed_af < 0.02" \
    --trio "denovo:mom.hom_ref && dad.hom_ref && kid.het && kid.GQ >= 5 && mom.GQ >= 5 && dad.GQ >= 5 && kid.AB >= 0.2 && kid.AB <= 0.8 && INFO.gnomad_popmax_af < 0.001 && !('gnomad_popmax_af_filter' in INFO) && INFO.topmed_af < 0.01" \
    --trio "recessive:recessive(kid, mom, dad)" \
    --trio "x_denovo:x_denovo(kid, mom, dad) && (variant.CHROM == 'chrX' || variant.CHROM == 'X')" \
    --trio "x_recessive:x_recessive(kid, mom, dad) && (variant.CHROM == 'chrX' || variant.CHROM == 'X')" \
    --trio "auto_dom:fake_auto_dom(kid, mom, dad) && variant.CHROM != 'chrX' && variant.CHROM != 'X' && INFO.gnomad_popmax_af < 0.001 && INFO.gnomad_nhomalt < 4 && INFO.topmed_af < 0.1" \
    --trio "comphet_side:comphet_side(kid, mom, dad) && INFO.gnomad_nhomalt_controls < 10 && INFO.gnomad_popmax_af < 0.005" \
    | bcftools csq -s - --ncsq 40 -g $gff -l -f $fasta - -o vcfs/$cohort.vcf

export SLIVAR_SUMMARY_FILE=$cohort.ch.summary.tsv
slivar compound-hets -f BCSQ  --sample-field comphet_side --sample-field denovo -p $ped -v vcfs/$cohort.vcf > vcfs/$cohort.ch.vcf

python plot-final-exome.py exome.summary.tsv exome.ch.summary.tsv
