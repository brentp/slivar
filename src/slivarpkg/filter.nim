import docopt
import os
import tables
import hts/vcf
import ./evaluator

proc main*(dropfirst:bool=false) =
  let doc = """
slivar -- variant expression for great good

Usage: slivar filter [options --pass-only --out-vcf <path> --vcf <path> --ped <path> --expr=<expressions>... --info=<expression>]

About:

    <expressions>...  as many name:expression pairs as required. an example would be:
    'low_quality:variant.FILTER != "PASS" || INFO.DP < 10 || variant.call_rate < 20'

    if a --info expression is specified, it is excuted first with access only to the variant and INFO objects. the boolean
    returned from this expression indicates whether the other expressions (in --expr) should be excecuted. This is an optimization
    to allow slivar to avoid loading all the sample data unless necessary.

Options:

  -v --vcf <path>       VCF/BCF
  --region <string>     optional region to limit evaluation. e.g. chr1 or 1:222-333
  -j --js <path>        path to javascript functions to expose to user
  -p --ped <path>       pedigree file with trio relations
  -o --out-vcf <path>   VCF/BCF
  --pass-only           only output variants that pass at least one of the filters [default: false]
  --expr <string>...    each trio expressions is applied to each trio where "mom", "dad", "kid" labels are available
  --info <string>       apply a filter using only variables from  the info field and variant attributes. If this filter
                        does not pass, the trio filters will not be applied.
  """

  var args: Table[string, docopt.Value]
  if dropfirst:
    var argv = commandLineParams()
    args = docopt(doc, argv=argv[1..argv.high])
  else:
    args = docopt(doc)

  if $args["--vcf"] == "nil":
    stderr.write_line "must specify the --vcf"
    quit doc
  if $args["--out-vcf"] == "nil":
    stderr.write_line "must specify the --out-vcf"
    quit doc
  var
    ivcf:VCF
    ovcf:VCF

  if not open(ivcf, $args["--vcf"], threads=1):
    quit "couldn't open:" & $args["--vcf"]

  var pass_only = bool(args["--pass-only"])

  ovcf.copy_header(ivcf.header)

  var tbl = ovcf.getExpressionTable(@(args["--trio"]), $args["--vcf"])
  doAssert ovcf.write_header


