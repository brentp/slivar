wget -O - ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/OsloUniversityHospital_Exome_GATK_jointVC_11242015/HG002-HG003-HG004.jointVC.filter.vcf | grep -v "##GATKCommandLine" | grep -v "GVCF" | head -10000 | perl -pe 's/(Sample_Diag-excap51-|-EEogPU)//g' | bgzip -l 9 -c > tests/ashk-trio.vcf.gz
echo "ash	HG002	HG003	HG004	1	0
ash	HG003	0	0	1	0
ash	HG004	0	0	2	0" > tests/ashk-trio.ped
