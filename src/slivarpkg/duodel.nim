import argparse
import math
import sets
import hts/vcf
import strformat
import strutils
import ./gnotate
import bpbiopkg/pedfile

type Duo = object
  kid: Sample
  parent: Sample
  parent_label: string # "mom" or "dad"
  i: int

template AB(sample_i:int, AD:seq[int32]): float32 =
  var a = AD[2*sample_i+1].float32
  a / max(1, a + AD[2*sample_i].float32)

template DP(sample_i: int, AD:seq[int32]): int32 =
  AD[2*sample_i+1] + AD[2*sample_i]

proc hq(sample_i: int, GQ: seq[int32], AD: seq[int32], alts: seq[int8], min_dp=9): bool =
  if sample_i.DP(AD) < min_dp: return false
  if GQ[sample_i] < 20: return false
  var ab = sample_i.AB(AD)
  case alts[sample_i]:
    of 0:
      return ab < 0.01
    of 1:
      return ab > 0.18 and ab < 0.82
    of 2:
      return ab > 0.99
    else:
      return false

type Transmit {.pure.} = enum
  LowQual = -1
  # 0, 0 or 1, 1
  Uninformative = 0
  No = 1
  KidHet = 2
  Maybe = 3

proc transmit(duo: Duo, GQ: seq[int32], AD: seq[int32], alts: seq[int8]): set[Transmit] =
  # mendelian error
  if not hq(duo.kid.i, GQ, AD, alts, 9): result.incl  Transmit.LowQual
  if not hq(duo.parent.i, GQ, AD, alts): result.incl Transmit.LowQual
  var k = alts[duo.kid.i]
  var p = alts[duo.parent.i]

  if k == 1:
    return result + {Transmit.KidHet}

  if (k == 0 and p == 2) or (k == 2 and p == 0):
    return result + {Transmit.No}

  if (k == 0 and p == 0) or (k == 2 and p == 2):
    return result + {Transmit.Uninformative}

  return result + {Transmit.Maybe}


proc duos(samples:seq[Sample], affected_only:bool): seq[Duo] =
  for s in samples:
    if affected_only and not s.affected: continue
    if s.mom != nil: result.add(Duo(kid:s, parent: s.mom, i: result.len, parent_label: "mom"))
    if s.dad != nil: result.add(Duo(kid:s, parent: s.dad, i: result.len, parent_label: "dad"))

type Candidate = object
  chrom: string
  start: int
  sample_i: int
  status: set[Transmit]
  i: int
  ads: array[4, int16]

proc near(group: seq[Candidate], other:Candidate, i_dist:int=3, genome_dist:int=10000): bool =
  for g in group:
    if abs(g.i - other.i) <= i_dist and abs(g.start - other.start) < genome_dist: return true
  return false

type CandidateStats = object
  # count of conflicting variants within the outer bounds
  hq_kid_hets: int
  kid_hets: int
  # excluding low-qual
  supporting: int
  maybe: int
  # all including low-qual
  supportingLowQ: int
  # nearby lowQ up/downstream
  nearbyLowQ: int
  nearbyN: int

proc stats(group: seq[Candidate], candidates: seq[Candidate]): CandidateStats =
  ## how many variants inside the event are bad.
  var used = initHashSet[int]()
  var imin = group[0].i
  var imax = group[group.high].i
  for site in group:
    used.incl(site.i)
    imin = min(imin, site.i)
    imax = max(imax, site.i)
  result.supporting = group.len

  var dist = 5000

  for i in imin..<imax:
    if i in used: continue
    if Transmit.KidHet in candidates[i].status:
      if Transmit.LowQual notin candidates[i].status:
        result.hq_kid_hets += 1
      result.kid_hets += 1
    if Transmit.LowQual in candidates[i].status and Transmit.No in candidates[i].status:
      result.supportingLowQ += 1
    if Transmit.Maybe in candidates[i].status:
      result.maybe += 1
    used.incl(i)

  var left = candidates[imin].start
  while imin > 0:
    imin -= 1
    if left - candidates[imin].start > dist: break
    if Transmit.LowQual in candidates[imin].status:
      result.nearbyLowQ += 1
      result.nearbyN += 1
    elif (Transmit.Maybe in candidates[imin].status) or (Transmit.KidHet in candidates[imin].status):
      result.nearbyN += 1

  var right = candidates[imax].start
  while imax < candidates.high:
    imax += 1
    if candidates[imax].start - right > dist: break
    if Transmit.LowQual in candidates[imax].status:
      result.nearbyLowQ += 1
      result.nearbyN += 1
    elif (Transmit.Maybe in candidates[imin].status) or (Transmit.KidHet in candidates[imin].status):
      result.nearbyN += 1

proc cull(duo: Duo, candidates: seq[Candidate], i_dist:int=20, genome_dist:int=30000, min_sites:int=2, min_size:int=20) =
   var seeds = newSeq[Candidate]()
   var nlq = 0
   for c in candidates:
    if Transmit.LowQual in c.status:
      nlq.inc
      continue
    if Transmit.No in c.status: seeds.add(c)

   var used = initHashSet[int]()
   for j, candidate in seeds:
    if candidate.i in used: continue
    used.incl(candidate.i)
    var group = @[candidate]
    for k, other in seeds:
      if other.i - seeds[seeds.high].i > i_dist: break
      if other.i in used: continue
      if group.near(other, i_dist, genome_dist):
        group.add(other)
        used.incl(other.i)
    if group.len < min_sites: continue
    if group[group.high].start - group[0].start < min_size: continue
    var stats = group.stats(candidates)
    if stats.hq_kid_hets / stats.supporting > 0.2: continue
    var parent = duo.parent_label
    echo &"{candidate.chrom}\t{group[0].start}\t{group[group.high].start + 1}\t{group.len}\t{parent}\t{group[group.high].i - group[0].i + 1}\t{duo.kid.id}\t{stats}"


template min16(i:int32):int16 =
  min(int16.high.int32, i).int16

include ./utils

proc main*(dropfirst:bool=false) =
  ## TODO: find homozygous deletions with ./.
  var p = newParser("slivar duodel"):
    help("""find denovo structural deletions in parent-child duos using non-transmission of SNPs
    see: https://github.com/brentp/slivar/wiki/finding-de-novo-deletions-in-parent-child-duos
    """)
    option("-p", "--ped", help="required ped file describing the duos in the VCF")
    option("-g", "--gnotate", help="optional gnotate file to check for flagged variants to exclude")
    option("-s", "--min-sites", default="3", help="minimum number of variants required to define a region (use 1 to output all putative deletions)")
    option("-S", "--min-size", default="50", help="minimum size in base-pairs of a region")
    option("-x", "--exclude", help="path to BED file of exclude regions e.g. (LCRs or self-chains)")
    flag("-a", "--affected-only", help="only output DEL calls for affected kids")
    arg("vcf", default="/dev/stdin", help="input SNP/indel VCF")

  var argv = commandLineParams()
  if len(argv) > 0 and argv[0] == "duo-del":
    argv = argv[1..argv.high]

  let opts = p.parse(argv)

  if opts.ped == "":
    echo p.help
    quit "--ped is required"


  var samples = parse_ped(opts.ped)
  var ivcf:VCF
  if not open(ivcf, opts.vcf, threads=2):
    quit "couldn't open VCF"
  samples = samples.match(ivcf)
  var duos = samples.duos(opts.affected_only)
  stderr.write_line &"[slivar] found {duos.len} duos"
  var gno:Gnotater
  if opts.gnotate != "":
    if not gno.open(opts.gnotate):
      quit "couldn't open gnotate file at:" & opts.gnotate
    gno.update_header(ivcf)

  var candidates = newSeq[seq[Candidate]](duos.len)

  var exclude: TableRef[string, Lapper[region]]
  var empty:seq[region]
  if opts.exclude != "":
    exclude = read_bed(opts.exclude)

  var GQ: seq[int32]
  var AD: seq[int32]
  var x: seq[int32]
  var MQ: seq[float32]
  var last_rid = -1
  var min_sites = parseInt(opts.min_sites)
  var min_size = parseInt(opts.min_size)
  for v in ivcf:

    if v.rid != last_rid:
      if last_rid != -1:
        for duo in duos:
            duo.cull(candidates[duo.i], min_sites=min_sites, min_size=min_size)
            candidates[duo.i].setLen(0)
      last_rid = v.rid

    if v.CHROM == "X" or v.CHROM == "chrX" or v.CHROM == "MT" or v.CHROM == "chrM" or v.CHROM == "chrMT": continue
    if v.FILTER notin ["PASS", "", "."]: continue
    if len(v.ALT) > 1: continue
    if exclude != nil and stripChr(v.CHROM) in exclude and exclude[stripChr(v.CHROM)].find(v.start, v.start + 1, empty):
      continue
    var LowQual = false
    if gno != nil and gno.annotate(v) and v.info.has_flag(gno.names[0] & "_filter"):
      LowQual = true

    if v.info.get("MQ", MQ) == Status.OK and MQ[0] < 10:
      LowQual = true

    if v.format.get("GQ", GQ) != Status.OK:
      quit "GQ: int32 required to use this tool"
    if v.format.get("AD", AD) != Status.OK:
      quit "AD: int32 required to use this tool"
    for v in AD.mitems:
      if v < 0: v = 0
    if AD.sum == 0: continue

    var alts = v.format.genotypes(x).alts

    for duo in duos:
      var m = transmit(duo, GQ, AD, alts)
      if LowQual:
        m.incl(Transmit.LowQual)
      var c = Candidate(chrom: $v.CHROM, start: v.start, status: m, i: candidates[duo.i].len)
      c.ads[0] = min16(AD[2*duo.kid.i])
      c.ads[1] = min16(AD[2*duo.kid.i+1])
      c.ads[2] = min16(AD[2*duo.parent.i])
      c.ads[3] = min16(AD[2*duo.parent.i+1])
      candidates[duo.i].add(c)

  for duo in duos:
    duo.cull(candidates[duo.i], min_sites=min_sites, min_size=min_size)
    #stderr.write_line "BREAKING EARLY"
    break

when isMainModule:
  main()
