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
  result = [".", "", ""]
  if ads.len == 0: return
  result[0] = &"{ads[2*sample.i] + ads[2*sample.i+1]}"
  if sample.dad != nil:
    result[1] = &"{ads[2*sample.dad.i] + ads[2*sample.dad.i+1]}"
  if sample.mom != nil:
    result[2] = &"{ads[2*sample.mom.i] + ads[2*sample.mom.i+1]}"

proc getAB(ads:var seq[int32], sample:Sample): array[3, string] =
  result = [".", "", ""]
  if ads.len == 0: return
  result[0] = &"{ads[2*sample.i+1].float32 / max(1, ads[2*sample.i] + ads[2*sample.i+1]).float32:g}"
  if sample.dad != nil:
    result[1] = &"{ads[2*sample.dad.i+1].float32 / max(1, ads[2*sample.dad.i] + ads[2*sample.dad.i+1]).float32:g}"
  if sample.mom != nil:
    result[2] = &"{ads[2*sample.mom.i+1].float32 / max(1, ads[2*sample.mom.i] + ads[2*sample.mom.i+1]).float32:g}"

template lookup(a:int8): string =
  $a
  #if a == -1: return "unknown"
  #const lookup = ["hom-ref", "het", "hom-alt"]
  #return lookup[a]

proc getGenotype(alts:seq[int8], sample:Sample): array[3, string] =
  result = [".", "", ""]
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

  columns: OrderedTableRef[string, int]

proc set_csq_fields(ivcf:VCF, field:string, gene_fields: var GeneIndexes, csq_columns: seq[string]) =
  gene_fields.gene = -1
  gene_fields.consequence = -1
  gene_fields.transcript = -1
  gene_fields.columns = newOrderedTable[string, int]()

  var desc = ivcf.header.get(field, BCF_HEADER_TYPE.BCF_HL_INFO)["Description"]
  var spl = (if "Format: '" in desc: "Format: '" else: "Format: ")
  var adesc = desc.split(spl)[1].split("'")[0].strip().strip(chars={'"', '\''}).multiReplace(("[", ""), ("]", ""), ("'", ""), ("*", "")).split("|")

  for v in adesc.mitems: v = v.toUpperAscii

  for cq in csq_columns:
    gene_fields.columns[cq] = adesc.find(cq.toUpperAscii)
    if gene_fields.columns[cq] == -1:
      raise newException(KeyError, &"[slivar] requested csq column '{cq}' not found in {adesc}")

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
      for c, ci in gene_fields.columns:
        key &= "/" & (if ci < toks.len: toks[ci] else: "")
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
    help("""convert filtered VCF to tab-separated-value spreadsheet for final filtering

create a --gene-description file e.g. from human phenotype ontology with:
  wget -qO - http://compbio.charite.de/jenkins/job/hpo.annotations.monthly/lastSuccessfulBuild/artifact/annotation/ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt | cut -f 2,3 > gene_desc.txt
or from clinvar with:
  wget -qO - ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/gene_condition_source_id | cut -f 2,5 | grep -v ^$'\t' > clinvar_gene_desc.txt
or gene->pLI with:
   wget -qO - https://storage.googleapis.com/gnomad-public/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz | zcat \
       | cut -f 1,21,24 | tail -n+2 | awk '{ printf("%s\tpLI=%.3g;oe_lof=%.5g\n", $1, $2, $3)}'

    """)
    option("-p", "--ped", default="", help="required ped file describing the trios in the VCF")
    option("-c", "--csq-field", help="INFO field containing the gene name and impact. Usually CSQ or BCSQ")
    option("--csq-column", multiple=true, help="CSQ sub-field(s) to extract (in addition to gene, impact, transcript). may be specified multiple times.")
    option("-s", "--sample-field", multiple=true, help="INFO field(s) that contains list of samples that have passed previous filters.\ncan be specified multiple times.")
    option("-g", "--gene-description", multiple=true, help="tab-separated lookup of gene (column 1) to description (column 2) to add to output. the gene is case-sensitive. can be specified multiple times")
    option("-i", "--info-field", multiple=true, help="INFO field(s) that should be added to output (e.g. gnomad_popmax_af)")
    option("-o", "--out", default="/dev/stdout", help="path to output tab-separated file")
    arg("vcf", help="input VCF", default="/dev/stdin")

  var argv = commandLineParams()
  if len(argv) > 0 and argv[0] == "tsv":
    argv = argv[1..argv.high]


  var opts = p.parse(argv)
  if opts.help:
    quit 0

  if opts.vcf.len == 0:
    opts.vcf.add("/dev/stdin")

  var extra: string
  if opts.ped == "":
    stderr.write_line "[slivar] no ped file specified, not able to determine family relationships."
  else:
    extra = "(sample,dad,mom)"

  var tsv_header = @["mode", "family_id", "sample_id", "chr:pos:ref:alt", "genotype" & extra]

  if opts.sample_field.len == 0:
    echo p.help
    quit "must specify at least one --sample-field"

  var ivcf:VCF
  if not open(ivcf, opts.vcf):
    quit "[slivar tsv] couldn't open vcf:" & opts.vcf

  var g2ds:seq[TableRef[string, string]]
  for gd in opts.gene_description:
    g2ds.add(gd.gene2description)

  for f in opts.info_field:
    doAssert ivcf.header.get(f, BCF_HEADER_TYPE.BCF_HL_INFO)["Type"] != ""
    tsv_header.add(f)

  var gene_fields :GeneIndexes
  var outfh: File
  var has_comphet = false
  if opts.out == "/dev/stdout":
    outfh = stdout
  else:
    if not open(outfh, opts.out, fmWrite):
      quit "couldn't open output file:" & opts.out

  if opts.csq_field != "":
    set_csq_fields(ivcf, opts.csq_field, gene_fields, opts.csq_column)
    tsv_header.add("gene")

  for f in opts.sample_field:
    doAssert ivcf.header.get(f, BCF_HEADER_TYPE.BCF_HL_INFO)["Type"] == "String"
    if f == "slivar_comphet":
      has_comphet = true

  var samples: seq[Sample]
  if opts.ped != "":
    samples = parse_ped(opts.ped)
    samples = samples.match(ivcf)
  else:
    samples = newSeq[Sample](ivcf.n_samples)
    for i, sid in ivcf.samples:
      samples[i] = Sample(id:sid, i:i)

  var sampleId2Obj = samples.toIndexLookup

  tsv_header.add("depths" & extra)
  tsv_header.add("allele_balance" & extra)
  if opts.csq_field != "":
    tsv_header.add("gene_impact_transcript")

    if opts.csq_column.len > 0:
      tsv_header[tsv_header.high] &= "_" & join(opts.csq_column,  "_")
    for i, g in g2ds:
      tsv_header.add(&"gene_description_{i+1}")
  outfh.write_line "#" & join(tsv_header, "\t")

  var str:string
  var xg:seq[int32]
  var ad:seq[int32]
  for v in ivcf:
    for f in opts.sample_field:
      if v.info.get(f, str) != Status.OK: continue
      if v.format.get("AD", ad) != Status.OK:
        ad.setLen(0)
      var alts = v.format.genotypes(xg).alts
      for osample_id in str.split(seps={','}):
        var sample_id = osample_id
        if has_comphet and f == "slivar_comphet":
          var tmp = sample_id.split("/", 1)
          sample_id = tmp[0]
        if sample_id notin sampleId2Obj: continue
        var sample = sampleId2Obj[sample_id]
        var line = @[f,sample.family_id, sample.id, &"""{v.CHROM}:{v.start+1}:{v.REF}:{join(v.ALT, ",")}"""]
        line.add(join(getGenotype(alts, sample), ",").replace(",,", ""))

        for f in opts.info_field:
          line.add(v.getField(f, ivcf))

        var genes: seq[string]
        if gene_fields.gene != -1:
          genes = v.get_gene_info(opts.csq_field, gene_fields, just_gene=true)
          line.add(join(genes, ";"))

        if has_comphet and f == "slivar_comphet":
          # add the compound het id to the "mode" so we can tell which variants
          # travel together in the spreadsheet.
          line[0] &= "_" & osample_id.split("/")[2]

        line.add(join(getDP(ad, sample), ",").replace(",,", ""))
        line.add(join(getAB(ad, sample), ",").replace(",,", ""))

        if gene_fields.gene != -1:
          line.add(join(v.get_gene_info(opts.csq_field, gene_fields), ";"))
          for g2d in g2ds:
            var ds = ""
            for gi, gene in genes:
              var dss = g2d.getOrDefault(gene)
              if dss == "": continue
              ds &= g2d.getOrDefault(gene)
              if gi < genes.high: ds &= ";;"
            line.add(ds.strip(chars={';'}))

        out_fh.write_line join(line, "\t")

when isMainModule:
  main()
