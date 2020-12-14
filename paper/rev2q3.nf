#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process reslivar {
    // new version of slivar has highest_impact_order
    publishDir 'rev2q3'
    input:
        tuple(path(vcf), val(name))
      
    output: tuple(path(output_file), val(name))
    script:
      output_file = "${name}-rev2q3.vcf"
      """
      slivar expr -v $vcf -o $output_file --alias <(echo "#")
      """
}

process count_highest_impact {
	publishDir 'rev2q3', mode:"copy"
    cache false
    input:
      tuple(val(vcf), val(name))
      val(order)
    output:
       path("*.txt")
    script:
        """
        python $HERE/rev2q3-count-highest-impact.py $order $vcf > ${name}.counts.txt
        """
}

process plot {
	publishDir 'rev2q3', mode:"copy"
    cache false
    input: val(txts)
    output: path("rev2q3.png")
    script:
        """
        python $HERE/rev2q3-plot.py ${txts.join(" ")}
        """
}


HERE = "/uufs/chpc.utah.edu/common/HIPAA/u6000771/Projects/2019/slivar-paper-analyses"
ORDER = "/uufs/chpc.utah.edu/common/HIPAA/u6000771/Projects/src/slivar/src/slivarpkg/default-order.txt"

workflow {

    vcfs = channel.from([["$HERE/vcfs/exome.impactful.vcf", "exome"], ["$HERE/vcfs/genome.impactful.vcf", "genome"]])
    
	reslivar(vcfs)
	count_highest_impact(reslivar.out, ORDER) 
    plot(count_highest_impact.out.collect())
    

}
