wget -cq https://ftp.ncbi.nih.gov/snp/pre_build152/organisms/human_9606_b151_GRCh38p7/VCF/00-All.vcf.gz
wget -cq https://ftp.ncbi.nih.gov/snp/pre_build152/organisms/human_9606_b151_GRCh38p7/VCF/00-All.vcf.gz.tbi

./scripts/topmed 00-All.vcf.gz dbsnp.pre.151.bcf
bcftools index --threads 3 dbsnp.pre.151.bcf

fasta=/data/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa
bcftools norm -m - -f $fasta -w 10000 dbsnp.pre.151.bcf -O b -o dbsnp.151.bcf

./src/slivar make-gnotate --field topmed_af -m "dbsnp 151 from https://ftp.ncbi.nih.gov/snp/pre_build152/organisms/human_9606_b151_GRCh38p7/VCF/00-All.vcf.gz " \
	--prefix topmed.hg38.dbsnp.151 \
	dbsnp.151.bcf
