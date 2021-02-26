nextflow.enable.dsl=2

js = "~/Projects/src/slivar/js/slivar-functions.js"

outdir = params.outdir

process slivar {
  publishDir "$outdir", mode: "copy" 
  cpus 2

  input: tuple(val(vcf), val(csi), val(tech), val(outdir), val(ped), path(zip), path(exclude), path(fasta), path(fai))
  output: tuple(path(slivar_vcf), path(slivar_csi), path(slivar_summary), path(slivar_ch_summary), path(slivar_ch_vcf), path(slivar_dn_summary), path(slivar_dn_vcf),
           path(slivar_impactful_vcf), path(slivar_impactful_csi), path(slivar_impactful_ch_vcf), path(slivar_impactful_ch_summary), path(slivar_impactful_summary),
           path(slivar_genic_vcf), path(slivar_genic_csi), path(slivar_genic_ch_vcf), path(slivar_genic_ch_summary), path(slivar_genic_summary)
            )

  script:
  slivar_dn_summary = "${tech}.dn.slivar.summary.txt"
  slivar_dn_vcf = "${tech}.dn.slivar.vcf.gz"

  slivar_vcf = "${tech}.slivar.vcf.gz"
  slivar_csi = "${tech}.slivar.vcf.gz.csi"
  slivar_ch_vcf = "${tech}.ch.slivar.vcf"
  slivar_summary = "${tech}.slivar.summary.txt"
  slivar_ch_summary = "${tech}.ch.slivar.summary.txt"

  slivar_impactful_ch_summary = "${tech}.ch.slivar.impactful.summary.txt"
  slivar_impactful_summary = "${tech}.slivar.impactful.summary.txt"
  slivar_impactful_vcf = "${tech}.slivar.impactful.vcf.gz"
  slivar_impactful_csi = "${tech}.slivar.impactful.vcf.gz.csi"
  slivar_impactful_ch_vcf = "${tech}.ch.slivar.impactful.vcf"

  slivar_genic_ch_summary = "${tech}.ch.slivar.genic.summary.txt"
  slivar_genic_summary = "${tech}.slivar.genic.summary.txt"
  slivar_genic_vcf = "${tech}.slivar.genic.vcf.gz"
  slivar_genic_csi = "${tech}.slivar.genic.vcf.gz.csi"
  slivar_genic_ch_vcf = "${tech}.ch.slivar.genic.vcf"

  """

export SLIVAR_SUMMARY_FILE=$slivar_dn_summary
pslivar expr --vcf $vcf --ped $ped \
        --fasta $fasta \
        -g $zip \
        --exclude $exclude \
        --pass-only \
        --js $js \
        --info "INFO.gnomad_popmax_af < 0.01 && variant.FILTER == 'PASS' && variant.ALT[0] != '*'" \
        --trio "denovo:denovo(kid, mom, dad)" \
        | bgzip -c -@ 4 > $slivar_dn_vcf


### ONLY IMPACTFUL

export SLIVAR_SUMMARY_FILE=$slivar_impactful_summary
pslivar expr --vcf $vcf --ped $ped \
        --fasta $fasta \
        -g $zip \
        --exclude $exclude \
        --pass-only \
        --js $js \
        --info "INFO.impactful && INFO.gnomad_popmax_af < 0.01 && variant.FILTER == 'PASS' && variant.ALT[0] != '*'" \
        --trio "denovo:denovo(kid, mom, dad) && INFO.gnomad_popmax_af < 0.001" \
        --trio "recessive:recessive(kid, mom, dad)" \
        --trio "x_denovo:x_denovo(kid, mom, dad) && (variant.CHROM == 'chrX' || variant.CHROM == 'X')" \
        --trio "x_recessive:x_recessive(kid, mom, dad) && (variant.CHROM == 'chrX' || variant.CHROM == 'X')" \
        --trio "auto_dom:fake_auto_dom(kid, mom, dad) && variant.CHROM != 'chrX' && variant.CHROM != 'X' && INFO.gnomad_popmax_af < 0.001 && INFO.gnomad_nhomalt < 10" \
        --trio "comphet_side:comphet_side(kid, mom, dad)" \
        | bgzip -c -@ 4 > $slivar_impactful_vcf

bcftools index --threads 3 $slivar_impactful_vcf;

export SLIVAR_SUMMARY_FILE=$slivar_impactful_ch_summary
slivar compound-hets --sample-field comphet_side --sample-field denovo -p $ped -v $slivar_impactful_vcf > $slivar_impactful_ch_vcf

### ONLY GENIC

export SLIVAR_SUMMARY_FILE=$slivar_genic_summary
pslivar expr --vcf $vcf --ped $ped \
        --fasta $fasta \
        -g $zip \
        --exclude $exclude \
        --pass-only \
        --js $js \
        --info "INFO.genic && INFO.gnomad_popmax_af < 0.01 && variant.FILTER == 'PASS' && variant.ALT[0] != '*'" \
        --trio "denovo:denovo(kid, mom, dad) && INFO.gnomad_popmax_af < 0.001" \
        --trio "recessive:recessive(kid, mom, dad)" \
        --trio "x_denovo:x_denovo(kid, mom, dad) && (variant.CHROM == 'chrX' || variant.CHROM == 'X')" \
        --trio "x_recessive:x_recessive(kid, mom, dad) && (variant.CHROM == 'chrX' || variant.CHROM == 'X')" \
        --trio "auto_dom:fake_auto_dom(kid, mom, dad) && variant.CHROM != 'chrX' && variant.CHROM != 'X' && INFO.gnomad_popmax_af < 0.001 && INFO.gnomad_nhomalt < 10" \
        --trio "comphet_side:comphet_side(kid, mom, dad)" \
        | bgzip -c -@ 4 > $slivar_genic_vcf

bcftools index --threads 3 $slivar_genic_vcf;

export SLIVAR_SUMMARY_FILE=$slivar_genic_ch_summary
slivar compound-hets --sample-field comphet_side --sample-field denovo -p $ped -v $slivar_genic_vcf > $slivar_genic_ch_vcf

### INCLUDES ALL VARIANTS, NOT JUST IMPACTFUL
export SLIVAR_SUMMARY_FILE=$slivar_summary
pslivar expr --vcf $vcf --ped $ped \
        --fasta $fasta \
        -g $zip \
        --exclude $exclude \
        --pass-only \
        --js $js \
        --info "INFO.gnomad_popmax_af < 0.01 && variant.FILTER == 'PASS' && variant.ALT[0] != '*'" \
        --trio "denovo:denovo(kid, mom, dad) && INFO.gnomad_popmax_af < 0.001" \
        --trio "recessive:recessive(kid, mom, dad)" \
        --trio "x_denovo:x_denovo(kid, mom, dad) && (variant.CHROM == 'chrX' || variant.CHROM == 'X')" \
        --trio "x_recessive:x_recessive(kid, mom, dad) && (variant.CHROM == 'chrX' || variant.CHROM == 'X')" \
        --trio "auto_dom:fake_auto_dom(kid, mom, dad) && variant.CHROM != 'chrX' && variant.CHROM != 'X' && INFO.gnomad_popmax_af < 0.001 && INFO.gnomad_nhomalt < 10" \
        --trio "comphet_side:comphet_side(kid, mom, dad)" \
        | bgzip -c -@ 4 > $slivar_vcf

bcftools index --threads 3 $slivar_vcf;

export SLIVAR_SUMMARY_FILE=$slivar_ch_summary
slivar compound-hets --sample-field comphet_side --sample-field denovo -p $ped -v $slivar_vcf > $slivar_ch_vcf


    """

}

workflow {
      slivar(tuple(params.vcf, params.vcf + ".csi", params.name, params.outdir, params.ped, params.zip, params.exclude, params.fasta, params.fasta + ".fai"))
}
