wget -qc ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.1/GRCh37/HG002_GRCh37_1_22_v4.1_draft_benchmark.vcf.gz
wget -qc ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.1/GRCh37/HG002_GRCh37_1_22_v4.1_draft_benchmark.vcf.gz.tbi
wget -qc ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.1/GRCh37/HG002_GRCh37_1_22_v4.1_draft_benchmark.bed.gz

set -euo pipefail

ref=/data/human/Homo_sapiens_assembly38.fasta
PATH=$PATH:~/src/rtg-tools/dist/rtg-tools-3.11-39691f9f/
sdf=hg38.sdf
#see: https://precision.fda.gov/challenges/10/view

if [[ -z $sdf ]]; then
	rtg format -o $sdf $ref
fi

vcf=DeepVariant.Illumina.HG002.vcf.gz
bed=HG002_GRCh38_1_22_v4.1_draft_benchmark.bed.gz
truth=HG002_GRCh38_1_22_v4.1_draft_benchmark.vcf.gz

zcat $vcf | bgzip -c > t.vcf.gz
mv t.vcf.gz $vcf
bcftools index -f --threads 4 --tbi $vcf

rm -rf eval


rtg vcfeval --bed-regions $bed -b $truth -o eval -t $sdf -c $vcf

python ab-truth.py eval/fp.vcf.gz eval/tp.vcf.gz eval/fn.vcf.gz
<<DONE
bcftools view --threads 3 -g "het" $vcf -O z -o het.vcf.gz
bcftools index --threads 3 --tbi het.vcf.gz

bcftools view --threads 3 -g "het" -i "FMT/GQ >= 20 && FMT/DP >= 10"  $vcf -O z -o gq20dp10.vcf.gz
bcftools index --threads 3 --tbi gq20dp10.vcf.gz

bcftools view --threads 3 -g "het" $vcf -O z -o bench.vcf.gz
bcftools index --threads 3 --tbi bench.vcf.gz

rm -rf eval-het
rm -rf eval-gq20dp10

rtg vcfeval --bed-regions $bed -b bench.vcf.gz -o eval-het -t $sdf -c het.vcf.gz
rtg vcfeval --bed-regions $bed -b bench.vcf.gz -o eval-gq20dp10 -t $sdf -c gq20dp10.vcf.gz


DONE
