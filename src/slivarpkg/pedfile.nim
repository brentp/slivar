import strutils
import math
import binaryheap
import sequtils
import algorithm
import strformat
import deques
import hashes
import sets
import tables
import hts/vcf

type
  Sample* = ref object
    family_id*: string
    id*: string
    paternal_id*: string
    maternal_id*: string
    mom*: Sample
    dad*: Sample
    sex*: int
    affected*: bool
    kids*: seq[Sample]
    i*:int

proc `$`*(s:Sample): string =
  return format(&"Sample(id:{s.id}, i:{s.i}, affected:{$s.affected})")

iterator siblings*(s:Sample): Sample =
  ## report the samples with same mom or dad.
  ## TODO: look at both mom and dad.
  ## TODO: handle cases where mom and dad ids
  var kids = if s.mom != nil: s.mom.kids elif s.dad != nil: s.dad.kids else: newSeq[Sample]()
  for kid in kids:
    if kid != s: yield kid
  if s.mom == nil and s.dad == nil and s.maternal_id != "" and s.paternal_id != "":
    # TODO: send in seq[Samples] and find samples with same ids
    discard

proc spouse*(s:Sample): Sample {.inline.} =
  if s.kids.len == 0: return nil
  var k = s.kids[0]
  if k.dad == s: return k.mom
  return k.dad

proc hash*(s:Sample): Hash =
  var h: Hash = 0

  h = h !& hash(s.family_id)
  h = h !& hash(s.id)
  result = !$h


type sample_path = tuple[sample:Sample, a_path:seq[Sample], b_path:seq[Sample]]

proc LCA(a:Sample, b:Sample, lcas: var HashSet[sample_path], a_path:seq[Sample], b_path:seq[Sample]) =

  var a_path = a_path
  var b_path = b_path

  if a == nil: return
  if b == nil: return

  if a == b:
    lcas.incl((a, a_path, b_path))

  block:
    var a_path = a_path & a.mom
    LCA(a.mom, b, lcas, a_path, b_path)
  block:
    var a_path = a_path & a.dad
    LCA(a.dad, b, lcas, a_path, b_path)
  block:
    var b_path = b_path & b.mom
    LCA(a, b.mom, lcas, a_path, b_path)
  block:
    var b_path = b_path & b.dad
    LCA(a, b.dad, lcas, a_path, b_path)

proc double_use(a:seq[Sample], b:seq[Sample]): bool {.inline.} =
  # a is a chain of samples up to common ancestor
  # b is chain of samples up to common ancestor
  # if a and b have the same tails, then we have double-use
  # and we don't count this chain
  # http://genetic-genealogy.co.uk/Toc115570135.html
  var ai = a.high
  var bi = b.high
  doAssert a[ai] == b[bi]
  ai -= 1
  bi -= 1
  if ai < 0 or bi < 0: return false

  # a:@[Sample(id:8451, i:-1, affected:false), Sample(id:8448, i:-1, affected:false)]
  # b:@[Sample(id:8683, i:-1, affected:false), Sample(id:8451, i:-1, affected:false), Sample(id:8448, i:-1, affected:false)]

  # go backwards through the chains
  result = false
  while ai > 0 or bi > 0:
    #echo "testing:", ai, ",", bi
    if a[ai] != b[bi]:
      #echo "bailing. ai:", ai, " bi:", bi
      return false
    return true
    #result = true
    #ai -= 1
    #bi -= 1
    #if ai <= 0 or bi <= 0: break

  when defined(superdebug):
    echo "ai:", ai, " bi:", bi, "dup?", result

proc is_full_sib_with(a:Sample, b:Sample): bool {.inline.} =
  return a.dad != nil and b.dad != nil and a.dad == b.dad and a.mom != nil and b.mom != nil and a.mom == b.mom

proc relatedness*(a:Sample, b:Sample): float =
  if a.family_id != b.family_id: return -1.0
  var r = initHashSet[sample_path](8)
  LCA(a, b, r, @[a], @[b])
  if len(r) == 0: return 0.0

  # now get a seq and sort by shortest totlal path-length
  var rs = newSeq[sample_path]()
  for x in r:
    when defined(superdebug):
      echo ">>>>>>>>>>>>>>>>>>>>>>> "
      echo "a:", x.a_path
      echo "b:", x.b_path
    if double_use(x.a_path, x.b_path):
      when defined(superdebug):
        echo "double"
        echo "<<<<<<<<<<<<<<<<<<<<"
      continue
    rs.add(x)
    when defined(superdebug):
      echo "not double"
      echo "<<<<<<<<<<<<<<<<<<<<"
    #echo " "
  for r in rs:
    result += pow(2'f64, -float64(r.a_path.len - 1 + r.b_path.len - 1))#*pow(1 + prel, 0.5)

  if result == 0.5 and a.is_full_sib_with(b):
    result = 0.49

proc parse_ped*(path: string, verbose:bool=true): seq[Sample] =
  result = new_seq_of_cap[Sample](10)

  var look = newTable[string,Sample]()

  var i = -1
  for line in lines(path):
    i += 1
    if line.len > 0 and line[0] == '#': continue
    if line.strip().len == 0: continue
    var toks = line.strip().split('\t')
    var likely_header = i == 0 and toks[0].toLowerAscii in ["family_id", "familyid", "famid", "kindred_id", "kindredid", "family id", "kindred id"]
    if toks.len < 6:
      stderr.write_line "[pedfile] error: expected at least 5 tab-delimited columns in ped file: " & path
      stderr.write_line "[pedfile] error: line was:" & $toks


    var s = Sample(family_id: toks[0], id: toks[1], kids:new_seq[Sample](), paternal_id: toks[2], maternal_id:toks[3], i: -1)
    s.affected = toks[5].toLowerAscii in  ["2", "affected", "yes"]
    if toks[4].toLowerAscii in ["XXXXXX", "unknown", "male", "female"]:
      s.sex = @["XXXXXX", "unknown", "male", "female"].find(toks[4].toLowerAscii) - 1
    else:
      try:
        s.sex = if toks[4] == ".": -9 else: parseInt(toks[4])
      except:
        if likely_header:
          stderr.write_line "skipping line as apparent header:" & line
          continue
        else:
          raise getCurrentException()
    result.add(s)
    look[s.id] = s

  for s in result:
    if s.paternal_id in look:
      s.dad = look[s.paternal_id]
      s.dad.kids.add(s)
    elif verbose and not (s.paternal_id in @[".", "-9", "", "0"]):
      stderr.write_line &"[pedfile] paternal_id: \"{s.paternal_id}\" referenced for sample \"{s.id}\" not found"
    if s.maternal_id in look:
      s.mom = look[s.maternal_id]
      s.mom.kids.add(s)
    elif verbose and not (s.maternal_id in @[".", "-9", "", "0"]):
      stderr.write_line &"[pedfile] maternal_id: \"{s.maternal_id}\" referenced for sample \"{s.id}\" not found"

proc match*(samples: seq[Sample], vcf:var VCF, verbose:bool=true): seq[Sample] =
  ## adjust the VCF samples and the samples to match
  var si = newTable[string,int]()
  for i, s in vcf.samples:
    si[s] = i
  # get samples in same order as VCF
  result = new_seqOfCap[Sample](len(si))
  var samplesById = newTable[string,Sample]()
  var missing = newSeq[string]()
  for sample in samples:
    sample.i = -1
    if not (sample.id in si):
      if verbose:
        missing.add(sample.id)
      continue
    sample.i = si[sample.id]
    samplesById[sample.id] = sample

  if missing.len > 0:
    sort(missing, system.cmp)
    stderr.write_line(&"[pedfile] ped samples: {join(missing, \",\")} not found in VCF samples")

  var sampleList = newSeq[string]()
  missing.setLen(0)
  for s in vcf.samples:
    if not (s in samplesById):
      if verbose:
        missing.add(s)
      continue
    sampleList.add(s)
  if missing.len > 0:
    sort(missing, system.cmp)
    stderr.write_line(&"[pedfile] VCF samples: {join(missing, \",\")} not found in pedigree file")

  vcf.set_samples(sampleList)
  for i, s in vcf.samples:
    var sample = samplesById[s]
    sample.i = i
    result.add(sample)

when isMainModule:
  import algorithm
  import unittest
  import times

  var samples = parse_ped("tests/testa.ped")
  var ovcf:VCF
  if not open(ovcf, "tests/test.vcf"):
    quit "bad"

  suite "pedfile test suite":

    #[
    test "that samples match those in vcf":
      var osamples = samples.match(ovcf)
      for i, s in osamples:
        check s.i == i
        check s.id == ovcf.samples[i]

    test "that reversed samples match those in vcf":
      var osamples = samples.reversed
      osamples = osamples.match(ovcf)
      for i, s in osamples:
        check s.i == i
        check s.id == ovcf.samples[i]

    test "siblings":
      check samples[0].id == "101976"
      for sib in samples[0].siblings:
        check sib.dad == samples[0].dad
        check sib.mom == samples[0].mom

    ]#

    test "dijkstra and relatedness":
      var k1 = Sample(family_id:"1", id:"kid1")
      var k2 = Sample(family_id:"1", id:"kid2")
      var mom = Sample(family_id:"1", id:"mom")
      var dad = Sample(family_id:"1", id:"dad")
      var uncle = Sample(family_id: "1", id:"uncle")
      var cousin = Sample(family_id: "1", id:"cousin")
      var gma = Sample(family_id:"1", id:"gma")
      var gpa = Sample(family_id:"1", id:"gma")
      var ggma = Sample(family_id:"1", id:"ggma")
      var ggpa = Sample(family_id:"1", id:"ggpa")
      var unrel = Sample(family_id:"1", id:"un")
      var extern = Sample(family_id:"xxx", id:"extern")
      k1.mom = mom
      k2.mom = mom
      k1.dad = dad
      k2.dad = dad
      dad.mom = gma
      dad.dad = gpa
      uncle.mom = gma
      uncle.dad = gpa
      uncle.kids.add(cousin)
      cousin.dad = uncle
      gma.kids.add(@[dad, uncle])
      mom.kids.add(@[k1, k2])
      dad.kids.add(@[k1, k2])
      ggma.kids.add(gma)
      gma.mom = ggma
      gma.dad = ggpa

      check relatedness(uncle, k1) == 0.25
      check relatedness(dad, k1) == 0.5
      check relatedness(k1, k2) == 0.49 # ^^^^^
      check relatedness(k2, mom) == 0.5
      check relatedness(k2, gma) == 0.25
      check relatedness(k2, ggma) == 0.125
      check relatedness(extern, gma) == -1'f64
      check relatedness(unrel, gma) == 0.0'f64
      check relatedness(k1, cousin) == 0.125

    test "relatedness":
      var k1 = Sample(family_id:"1", id:"kid1")
      var dad = Sample(family_id:"1", id:"kid1")
      k1.dad = dad
      dad.kids.add(k1)
      var fam = @[k1, dad]

      check 0.5 == relatedness(k1, dad)

    test "relatedness with no relations":
      echo "xxx"
      var a = Sample(family_id:"1", id:"a")
      var b = Sample(family_id:"2", id:"b")
      check relatedness(a, b) == -1'f64
      b.family_id = "1"
      check relatedness(a, b) == -0'f64


    test "relatedness with 603 ceph samples":
      var samples = parse_ped("tests/ceph.ped")
      var t = cpuTime()
      var n = 0

      for i, sampleA in samples[0..<samples.high]:
        for j, sampleB in samples[i + 1..samples.high]:
          var rel = relatedness(sampleA, sampleB)
          doAssert -1'f64 <= rel
          if rel > 0.5:
            echo &"{$sampleA} {$sampleB} {$rel}"
          n += 1
          if n mod 50000 == 0: echo "tested:", n

      echo &"time for {n} calculations: {cpuTime() - t:.1f}"

  test "ped file with kindred header":

    var fh:File
    check open(fh, "__k.ped", fmWrite)
 
    fh.write_line("Kindred_Id\tSample_ID\tPaternal_ID\tMaternal_ID\tSex\tPhenotype")
    fh.write_line("Kindred_Id\ta\tdad\tmom\t1\t1")

    fh.close()

    var samples = parse_ped("__k.ped")
    check samples.len == 1
