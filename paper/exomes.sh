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

slivar expr --vcf $bcf --ped $ped \
    --pass-only \
    -o vcfs/$cohort.dn.bcf \
    --info "variant.FILTER == 'PASS'" \
    --trio "dn_naive:mom.alts == 0 && dad.alts == 0 && kid.alts == 1" \
    --trio "dn_ab_gq:mom.alts == 0 && dad.alts == 0 && kid.alts == 1 && kid.GQ >= 5 && mom.GQ >= 5 && dad.GQ >= 5 && kid.AB >= 0.2 && kid.AB <= 0.8" \
    --trio "dn_gnomad:mom.alts == 0 && dad.alts == 0 && kid.alts == 1 && kid.GQ >= 5 && mom.GQ >= 5 && dad.GQ >= 5 && kid.AB >= 0.2 && kid.AB <= 0.8 && INFO.gnomad_popmax_af < 0.001 && !('gnomad_popmax_af_filter' in INFO)" \
    --trio "dn:mom.alts == 0 && dad.alts == 0 && kid.alts == 1 && kid.GQ >= 5 && mom.GQ >= 5 && dad.GQ >= 5 && kid.AB >= 0.2 && kid.AB <= 0.8 && INFO.gnomad_popmax_af < 0.001 && !('gnomad_popmax_af_filter' in INFO) && INFO.topmed_af < 0.02" 


#bcftools csq -s - --ncsq 40 -g $gff -l -f $fasta vcfs/$cohort.bcf -O u \
#  | slivar compound-hets -f BCSQ -i 2 --sample-field AR_filtered --sample-field denovo -p $ped > vcfs/$cohort.ch.vcf
