set -euo pipefail
gff=/data/human/Homo_sapiens.GRCh38.96.chr.gff3.gz
bcf=/home/brentp/src/varx/var.bcf
ped=/home/brentp/src/varx/var.ped
fasta=/data/human/GRCh38_full_analysis_set_plus_decoy_hla.fa
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

./paper/dn_roc ../varx/var.ped paper/vcfs/exome.annotated.bcf 0 > roc-gq0.txt
./paper/dn_roc ../varx/var.ped paper/vcfs/exome.annotated.bcf 5 > roc-gq5.txt
./paper/dn_roc ../varx/var.ped paper/vcfs/exome.annotated.bcf 10 > roc-gq10.txt

python paper/plot_ab_roc.py roc-gq{0,5,10}.txt

exit
DONE

bcf=vcfs/$cohort.annotated.bcf

<<DENOVO
export SLIVAR_SUMMARY_FILE=$cohort.dn.summary.tsv

slivar expr --vcf $bcf --ped $ped \
    --pass-only \
    -o vcfs/$cohort.dn.bcf \
    --info "variant.FILTER == 'PASS' && variant.ALT[0] != '*'" \
    --trio "dn_ab_gq:mom.hom_ref && dad.hom_ref && kid.het && kid.GQ >= 5 && mom.GQ >= 5 && dad.GQ >= 5 && kid.AB >= 0.2 && kid.AB <= 0.8" \
    --trio "dn:mom.hom_ref && dad.hom_ref && kid.het && kid.GQ >= 5 && mom.GQ >= 5 && dad.GQ >= 5 && kid.AB >= 0.2 && kid.AB <= 0.8 && INFO.gnomad_popmax_af < 0.001 && !('gnomad_popmax_af_filter' in INFO) && INFO.topmed_af < 0.01"  \

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
    --trio "comphet_side:comphet_side(kid, mom, dad) && INFO.gnomad_nhomalt_controls < 10 && INFO.gnomad_popmax_af < 0.005" \
    | bcftools csq -s - --ncsq 40 -g $gff -l -f $fasta - -o vcfs/$cohort.vcf

export SLIVAR_SUMMARY_FILE=$cohort.ch.summary.tsv
slivar compound-hets -f BCSQ  --sample-field comphet_side --sample-field denovo -p $ped -v vcfs/$cohort.vcf > vcfs/$cohort.ch.vcf

python plot-final-exome.py exome.summary.tsv exome.ch.summary.tsv
