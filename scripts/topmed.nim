import hts/vcf

import strutils
import os

var ivcf:VCF
var ovcf:VCF

if not open(ivcf, paramStr(1)):
  quit "couldnot open input VCF"

if not open(ovcf, paramStr(2), mode="wb"):
  quit "couldnot open input VCF"

doAssert ivcf.header.add_info("topmed_af", "A", "Float", "topmed alternate allele frequency") == Status.OK
ovcf.copy_header(ivcf.header)
doAssert ovcf.write_header

var tm: string
var tmf: seq[float32]
for v in ivcf:

  if v.info.get("TOPMED", tm) != Status.OK: continue

  var tms = tm.split(",")
  tmf.setLen(tms.len - 1)
  for i, s in tms:
    if i == 0: continue
    tmf[i - 1] = parseFloat(s).float32

  if v.info.set("topmed_af", tmf) != Status.OK:
    quit "couldn't set variant"

  doAssert ovcf.write_variant(v)

ovcf.close()
