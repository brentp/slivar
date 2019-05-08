import argparse
import strutils
import os
import hts/vcf
import strformat
import bpbiopkg/pedfile
import tables

proc toIndexLookup(samples:seq[Sample]): TableRef[string,Sample] =
  result = newTable[string,Sample]()
  for s in samples:
    result[s.id] = s

proc getDP(ads:var seq[int32], sample:Sample): array[3, string] =
  result = [".", ".", "."]
  if ads.len == 0: return
  result[0] = &"{ads[2*sample.i] + ads[2*sample.i+1]}"
  if sample.dad != nil:
    result[1] = &"{ads[2*sample.dad.i] + ads[2*sample.dad.i+1]}"
  if sample.mom != nil:
    result[2] = &"{ads[2*sample.mom.i] + ads[2*sample.mom.i+1]}"

proc getAB(ads:var seq[int32], sample:Sample): array[3, string] =
  result = [".", ".", "."]
  if ads.len == 0: return
  result[0] = &"{ads[2*sample.i+1].float32 / (ads[2*sample.i] + ads[2*sample.i+1]).float32:g}"
  if sample.dad != nil:
    result[1] = &"{ads[2*sample.dad.i+1].float32 / (ads[2*sample.dad.i] + ads[2*sample.dad.i+1]).float32:g}"
  if sample.mom != nil:
    result[2] = &"{ads[2*sample.mom.i+1].float32 / (ads[2*sample.mom.i] + ads[2*sample.mom.i+1]).float32:g}"

template lookup(a:int8): string =
  $a
  #if a == -1: return "unknown"
  #const lookup = ["hom-ref", "het", "hom-alt"]
  #return lookup[a]

proc getGenotype(alts:seq[int8], sample:Sample): array[3, string] =
  result = [".", ".", "."]
  result[0] = lookup(alts[sample.i])
  if sample.dad != nil:
    result[1] = lookup(alts[sample.dad.i])
  if sample.mom != nil:
    result[2] = lookup(alts[sample.mom.i])

proc getField(v:Variant, field:string, ivcf:VCF): string =
  case ivcf.header.get(field, BCF_HEADER_TYPE.BCF_HL_INFO)["Type"]
  of "String":
    discard v.info.get(field, result)
    return
  of "Integer":
    var x: seq[int32]
    if v.info.get(field, x) != Status.OK:
      return "."
    var xs = newSeqOfCap[string](x.len)
    for val in x:
      xs.add($val)
    return join(xs, ",")
  of "Float":
    var x: seq[float32]
    if v.info.get(field, x) != Status.OK:
      return "."
    var xs = newSeqOfCap[string](x.len)
    for val in x:
      xs.add(&"{val:0g}")
    return join(xs, ",")
  of "Flag":
    if v.info.has_flag(field): return field
    else: return ""
  else:
    stderr.write_line "[slivar] type unknown for field:" & field
    return "."

type GeneIndexes = object
  gene: int
  consequence: int
  transcript: int

proc set_csq_fields(ivcf:VCF, field:string, gene_fields: var GeneIndexes) =
  gene_fields.gene = -1
  gene_fields.consequence = -1
  gene_fields.transcript = -1
  var desc = ivcf.header.get(field, BCF_HEADER_TYPE.BCF_HL_INFO)["Description"]
  var adesc = desc.split("Format:")[1].strip().strip(chars={'"', '\''}).multiReplace(("[", ""), ("]", ""), ("'", ""), ("*", "")).split("|")
  for v in adesc.mitems: v = v.toUpperAscii
  for check in ["SYMBOL", "GENE"]:
    gene_fields.gene = adesc.find(check)
    if gene_fields.gene != -1: break
  for check in ["CONSEQUENCE"]:
    gene_fields.consequence = adesc.find(check)
    if gene_fields.consequence != -1: break
  for check in ["FEATURE", "TRANSCRIPT"]:
    gene_fields.transcript = adesc.find(check)
    if gene_fields.transcript != -1: break

  if gene_fields.gene == -1:
    quit &"[slivar] unable to find gene field in {field}"
  if gene_fields.consequence == -1:
    quit &"[slivar] unable to find consequence field in {field}"
  if gene_fields.transcript == -1:
    quit &"[slivar] unable to find transcript field in {field}"


proc get_gene_info(v:Variant, csq_field_name:string, gene_fields:GeneIndexes, just_gene:bool=false): seq[string] =
  ## get the gene_names and consequences for each transcript.
  var s = ""
  if v.info.get(csq_field_name, s) != Status.OK:
    return

  for tr in s.split(','):
    var toks = tr.split('|')
    var key = toks[gene_fields.gene].strip()
    if not just_gene:
      key &= "/" & toks[gene_fields.consequence] & "/" & toks[gene_fields.transcript]
    if key.strip().len == 0 or key in result: continue
    result.add(key)

proc gene2description(fname:string): TableRef[string,string] =
  result = newTable[string, string](1024)
  for line in fname.lines:
    if line[0] == '#': continue
    var toks = line.strip(chars={'\r', '\n', ' '}).split("\t")
    if toks[0] in  result:
      result[toks[0]] &= ';' & toks[1]
    else:
      result[toks[0]] = toks[1]

proc main*(dropfirst:bool=false) =
  var p = newParser("slivar tsv"):
    help("""convert filtered VCF to spreadsheet for final filtering

create a --gene-description file e.g. from human phenotype ontology with:
  wget -qO - http://compbio.charite.de/jenkins/job/hpo.annotations.monthly/lastSuccessfulBuild/artifact/annotation/ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt | cut -f 2,3 > gene_desc.txt
or from clinvar with:
  wget -qO - ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/gene_condition_source_id | cut -f 2,5 | head | grep -v ^$'\t' > clinvar_gene_desc.txt
    """)
    option("-p", "--ped", default="", help="required ped file describing the trios in the VCF")
    option("-c", "--csq-field", help="INFO field containing the gene name and impact. Usually CSQ or BCSQ")
    option("-s", "--sample-field", multiple=true, help="INFO field(s) that contains list of samples (kids) that have passed previous filters.\ncan be specified multiple times.")
    option("-g", "--gene-description", help="tab-separated lookup of gene (column 1) to description (column 2) to add to output. the gene is case-sensitive")
    option("-i", "--info-field", multiple=true, help="INFO field(s) that should be added to output (e.g. gnomad_popmax_af)")
    option("-o", "--out-vcf", default="/dev/stdout", help="path to output tab-separated file")
    arg("vcf", default="/dev/stdin", help="input VCF")

  var argv = commandLineParams()
  if len(argv) > 0 and argv[0] == "tsv":
    argv = argv[1..argv.high]

  var tsv_header = @["mode", "family_id", "sample_id", "chr:pos:ref:alt", "genotype(kid,dad,mom)"]

  let opts = p.parse(argv)
  if opts.help:
    quit 0
  if opts.ped == "":
    echo p.help
    quit "--ped is required"
  if opts.sample_field.len == 0:
    echo p.help
    quit "must specify at least one --sample-field"


  var ivcf:VCF
  if not open(ivcf, opts.vcf):
    quit "[slivar tsv] couldn't open vcf:" & opts.vcf

  var g2d:TableRef[string, string]
  if opts.gene_description != "":
    g2d = opts.gene_description.gene2description

  for f in opts.info_field:
    doAssert ivcf.header.get(f, BCF_HEADER_TYPE.BCF_HL_INFO)["Type"] != ""
    tsv_header.add(f)

  var gene_fields :GeneIndexes

  if opts.csq_field != "":
    set_csq_fields(ivcf, opts.csq_field, gene_fields)
    tsv_header.add("gene")

  for f in opts.sample_field:
    doAssert ivcf.header.get(f, BCF_HEADER_TYPE.BCF_HL_INFO)["Type"] == "String"

  var samples = parse_ped(opts.ped)
  samples = samples.match(ivcf)
  var sampleId2Obj = samples.toIndexLookup

  tsv_header.add("depths(kid,mom,dad)")
  tsv_header.add("allele_balance(kid,mom,dad)")
  if opts.csq_field != "":
    tsv_header.add("gene_impact_transcript")
    if g2d != nil:
      tsv_header.add("gene_description")
  echo "#" & join(tsv_header, "\t")

  var str:string
  var xg:seq[int32]
  var ad:seq[int32]
  for v in ivcf:
    for f in opts.sample_field:
      if v.info.get(f, str) != Status.OK: continue
      if v.format.get("AD", ad) != Status.OK:
        ad.setLen(0)
      var alts = v.format.genotypes(xg).alts
      for sample_id in str.split(seps={','}):
        if sample_id notin sampleId2Obj: continue
        var sample = sampleId2Obj[sample_id]
        var line = @[f,sample.family_id, sample.id, &"""{v.CHROM}:{v.start+1}:{v.REF}:{join(v.ALT, ",")}"""]
        line.add(join(getGenotype(alts, sample), ","))

        for f in opts.info_field:
          line.add(v.getField(f, ivcf))

        var genes: seq[string]
        if gene_fields.gene != -1:
          genes = v.get_gene_info(opts.csq_field, gene_fields, just_gene=true)
          line.add(join(genes, ";"))

        line.add(join(getDP(ad, sample), ","))
        line.add(join(getAB(ad, sample), ","))

        if gene_fields.gene != -1:
          line.add(join(v.get_gene_info(opts.csq_field, gene_fields), ";"))
          if g2d != nil:
            var ds = ""
            for gene in genes:
              var dss = g2d.getOrDefault(gene)
              if dss == "": continue
              ds &= g2d.getOrDefault(gene) & ";;"
            line.add(ds)

        echo join(line, "\t")


when isMainModule:
  main()
