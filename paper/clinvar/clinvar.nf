nextflow.enable.dsl=2

params.fasta = false
if(!params.fasta) {
 exit 1, "--fasta argument required"
}

params.genome_build = "GRCh38.99"

process get_clinvar {
    publishDir "results/"
    input: 
        val(url)
    output: path(output_file)
    script:
    output_file = file(url).baseName + ".gz"

    """
    wget -qc $url
    """
}

process get_gff {
    publishDir "results/"
    input: 
        val(url)
    output: path(output_file)
    script:
    output_file = file(url).baseName + ".gz"

    """
    wget -qc $url
    """
}

process csq {
    publishDir "results/"
    input:
       val("vcf")
       val("gff")
       val("fasta")
    output: 
       path(output_file)
       path(output_csi)
    script:
    output_file = file(file(vcf).baseName).baseName + ".norm.csq.vcf.gz"
    output_csi =  file(file(vcf).baseName).baseName + ".norm.csq.vcf.gz.csi"
    """
    bcftools index --threads 3 $vcf || true
    bcftools norm --threads 3 -m - -w 10000 -f $fasta -O v $vcf \
        | snpEff -Xmx4G eff -noStats ${params.genome_build} \
        | bcftools csq --threads 2 -O z -o ${output_file} -s - --ncsq 40 -g $gff -l -f $fasta -
    bcftools index --threads 3 $output_file
    """
}

process slivar {
    publishDir "results/"
    input: 
      val("vcf")
      val("csi")
    output: path(output_file)

    script:
    output_file = file(file(vcf).baseName).baseName + ".slivar.vcf.gz"
    """
    slivar expr -g ~/src/slivar/gnomad.hg38.genomes.v3.fix.zip  -v $vcf -o $output_file
    """
}

process pathoimpact {
    publishDir "results/"
    input:
        val("slivar_vcf")
    output: path(output_file)
    script:
    output_file = file(slivar_vcf).baseName + ".txt"
    """
    pathoimpact ${slivar_vcf} > ${output_file}
    """
}


workflow {
    get_clinvar("https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2020/clinvar_20201031.vcf.gz")
    //get_gff("http://ftp.ensembl.org/pub/current_gff3/homo_sapiens/Homo_sapiens.GRCh38.102.gff3.gz")
    get_gff("http://ftp.ensembl.org/pub/release-99/gff3/homo_sapiens/Homo_sapiens.GRCh38.99.gff3.gz")
    
    csq(get_clinvar.out, get_gff.out, params.fasta)

    slivar(csq.out)
    pathoimpact(slivar.out)

}
