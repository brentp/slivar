import hts/vcf
import strutils
import os

var vcf_path = paramStr(1)

var ivcf:VCF

if not ivcf.open(vcf_path):
  quit "could not open vcf:" & vcf_path

var clnsig:string

var patho_impactful = 0
var patho_notimpactful = 0

for v in ivcf:
  let st = v.info.get("CLNSIG", clnsig)
  if st == Status.NotFound:
    continue
  if st != Status.OK:
    quit "bad info?" & $st

  clnsig = toLowerAscii(clnsig)
  if "conflicting_interpretations_of_pathogenicity" in clnsig: continue

  var is_pathogenic = "pathogenic" in clnsig
  if not is_pathogenic: continue
  var impactful = v.info.has_flag("impactful")
  if impactful:
    patho_impactful += 1
  else:
    patho_notimpactful += 1
    #echo v.tostring()


echo "not impactful:", patho_notimpactful
echo "impactful:", patho_impactful
