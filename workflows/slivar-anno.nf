nextflow.enable.dsl=2

params.help = false
params.genome_build = "38"
params.outdir = false
if (params.help) {
    log.info("""

Usage
-----
    nextflow run -resume \
        -profile local \
        -entry annotate \
        -with-singularity docker://brentp/slivar-anno:v0.2.1 \
        --genome_build 37 \
        --bcf $vcf \
        --name $cohort \
        --fasta $ref \
        slivar-anno.nf

Required arguments:
-------------------

    --bcf      path to bcf to annotate
    --name     name for project
    --genome_build    genome-build. must be "38" or "37"
    --fasta    path to reference fasta file
	""")
}
   

bcf = file(params.bcf, checkIfExists: true)
name = params.name
fasta = file(params.fasta, checkIfExists: true)
outdir = "."
if (params.outdir) {
  outdir = params.outdir
}

if (params.genome_build == 38) {
  snpeffbuild = "GRCh38.99"
  gff = "/opt/slivar/Homo_sapiens.GRCh38.102.chr.gff3.gz"
  zip = "/opt/slivar/gnomad.hg38.genomes.v3.fix.zip"
} else if (params.genome_build == 37) {
  snpeffbuild = "GRCh37.75"
  gff = "/opt/slivar/Homo_sapiens.GRCh37.87.gff3.gz"
  zip = "/opt/slivar/gnomad.hg37.zip"
} else {
  log.info("${params.genome_build}")
   exit 1, "--genome_build must be 37 or 38"
}

log.info("$outdir $gff $zip $snpeffbuild $bcf $fasta ")


process annotate_vcf {

  publishDir "$outdir", mode: "copy"
  maxRetries 0

  cpus 5

  //container 'docker://brentp/slivar-anno:v0.2.1'

  input: tuple(val(bcf), val(name))
  output: tuple(path(output_file), path(output_csi), val(name))
  script:
  output_file = "${name}.anno.bcf"
  output_csi = "${name}.anno.bcf.csi"
  """
  set -euo pipefail
  bcftools norm --threads 3 -m - -w 10000 -f $fasta -O v $bcf \
    | snpEff -Xmx4G eff -noStats ${snpeffbuild} \
    | bcftools csq --threads 3 -s - --ncsq 40 -g $gff -l -f $fasta - -O b \
    | slivar expr -g $zip -o $output_file -v -

  bcftools index --threads 8 $output_file
  """                                                                                                                                                                                                                   
}

workflow annotate {
    annotate_vcf(tuple(bcf, name)) | view
}
