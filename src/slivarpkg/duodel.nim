import argparse
import algorithm
import math
import stats
import sets
import hts/vcf
import strformat
import strutils
import ./gnotate
import pedfile

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

proc hq*(sample_i: int, GQ: seq[int32], AD: seq[int32], alts: seq[int8], min_dp=9): bool {.inline.} =
  if GQ[sample_i] < 20: return false
  if sample_i.DP(AD) < min_dp: return false
  var ab = sample_i.AB(AD)
  case alts[sample_i]:
    of 0:
      return ab < 0.01
    of 1:
      return ab > 0.3 and ab < 0.7
    of 2:
      return ab > 0.98
    else:
      return false

type Transmit {.pure.} = enum
  LowQual
  # 0, 0 or 1, 1
  Uninformative
  No
  KidHet
  ParentHet
  ParentHomAlt
  Maybe

proc transmit(duo: Duo, GQ: seq[int32], AD: seq[int32], alts: seq[int8]): set[Transmit] =
  # mendelian error
  var k = alts[duo.kid.i]
  var p = alts[duo.parent.i]
  if k == -1 or p == -1:
    result.incl Transmit.Uninformative


  if not hq(duo.kid.i, GQ, AD, alts, 9): result.incl  Transmit.LowQual
  if not hq(duo.parent.i, GQ, AD, alts): result.incl Transmit.LowQual
  if p == 1:
    result.incl Transmit.ParentHet
  elif p == 2:
    result.incl Transmit.ParentHomAlt

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

type Site = object
  chrom*: string
  start*: int64
  sample_i*: int
  status*: set[Transmit]
  i*: int # index of sites used in this duo
  dps_i*: int # index into DPs seq
  ads*: array[4, uint16]

proc near(group: seq[Site], other:Site, i_dist:int=3, genome_dist:int=10000): bool =
  for g in group:
    if abs(g.i - other.i) <= i_dist and abs(g.start - other.start) < genome_dist: return true
  return false

type SiteStats = object
  # count of conflicting variants within the outer bounds
  hq_kid_hets: int
  kid_hets: int
  # excluding low-qual
  hq_parent_hets: int
  hq_parent_hom_alts: int

  supporting: int
  maybe: int
  # all including low-qual
  supportingLowQ: int
  # nearby lowQ up/downstream
  nearbyLowQ: int
  nearbyN: int

proc stats(group: seq[Site], sites: seq[Site]): SiteStats =
  ## how many variants inside the event are bad.
  var used = initHashSet[int]()
  var imin = group[0].i
  var imax = group[group.high].i
  for site in group:
    #used.incl(site.i)
    imin = min(imin, site.i)
    imax = max(imax, site.i)
  result.supporting = group.len

  var dist = 5000

  for i in imin..<imax:
    if i in used: continue
    used.incl(i)
    if Transmit.KidHet in sites[i].status:
      if Transmit.LowQual notin sites[i].status:
        result.hq_kid_hets += 1
      result.kid_hets += 1
    if Transmit.LowQual notin sites[i].status:
      if Transmit.ParentHet in sites[i].status:
        result.hq_parent_hets += 1
      if Transmit.ParentHomAlt in sites[i].status:
        result.hq_parent_hom_alts += 1

    if Transmit.LowQual in sites[i].status and Transmit.No in sites[i].status:
      result.supportingLowQ += 1
    if Transmit.Maybe in sites[i].status:
      result.maybe += 1

  var left = sites[imin].start
  while imin > 0:
    imin -= 1
    if left - sites[imin].start > dist: break
    if Transmit.LowQual in sites[imin].status:
      result.nearbyLowQ += 1
      result.nearbyN += 1
    elif (Transmit.Maybe in sites[imin].status) or (Transmit.KidHet in sites[imin].status):
      result.nearbyN += 1

  var right = sites[imax].start
  while imax < sites.high:
    imax += 1
    if sites[imax].start - right > dist: break
    if Transmit.LowQual in sites[imax].status:
      result.nearbyLowQ += 1
      result.nearbyN += 1
    elif (Transmit.Maybe in sites[imin].status) or (Transmit.KidHet in sites[imin].status):
      result.nearbyN += 1

proc normalize*(DPs: var seq[seq[uint16]]): seq[seq[float32]] =
  ## first do within-sample normalization, then across sample normalization
  result = newSeqOfCap[seq[float32]](DPs.len)
  if DPs.len == 0: return
  var n_samples = DPs[0].len
  var rss = newSeq[RunningStat](n_samples)
  for row in DPs:
    result.add(newSeq[float32](n_samples))
    for j, r in row:
      result[result.high][j] = r.float32
      if r > 0'u16:
        rss[j].push(r.float32)

  var means = newSeq[float32](n_samples)
  for i, r in rss:
    means[i] = max(1, r.mean)

  for row in result.mitems:
    # intra-sample normalization with samples that have > 0 coverage
    var S = 0'f64
    var n: int
    for i, v in row.mpairs:
      # divide by column mean
      v = v / means[i]
      if v > 0:
        S += v.float64
        n += 1
    # inter-sample normalization, divide by row mean
    var rowMean = (S / n.float64).float32
    for v in row.mitems:
      v = v / rowMean

proc get(s:Sample, DPs: var seq[seq[float32]]): seq[float32] =
  result = newSeqOfCap[float32](DPs.len)
  for row in DPs:
    result.add(row[s.i])

proc median(vals: seq[float32]): float32 =
  if vals.len == 0: return 1
  var vals = vals
  sort(vals, system.cmp[float32])
  return vals[int(vals.high / 2)]

proc cull(duo: Duo, sites: seq[Site], DPs: var seq[seq[float32]], i_dist:int=25, genome_dist:int=200000, min_sites:int=2, min_size:int=20) =

  var dp_norm_kid = duo.kid.get(DPs)
  var dp_norm_parent = duo.parent.get(DPs)

  var seeds = newSeq[Site]()
  var nlq = 0
  for c in sites:
    if Transmit.LowQual in c.status:
      nlq.inc
      continue
    # only high-quality non-transmissions can seed.
    if Transmit.No in c.status: seeds.add(c)

  var used = initHashSet[int]()
  shallow(seeds)
  for j, site in seeds:
    if site.i in used: continue
    used.incl(site.i)
    var group = @[site]
    for k, other in seeds[j..seeds.high]:
      if other.i - seeds[seeds.high].i > i_dist: break
      if other.i in used: continue
      if group.near(other, i_dist, genome_dist):
        group.add(other)
    if group.len < min_sites: continue
    if group[group.high].start - group[0].start < min_size: continue
    for g in group:
      used.incl(g.i)
    var stats = group.stats(sites)
    if stats.hq_kid_hets / stats.supporting > 0.1: continue
    #for s in group:
    var kdps: seq[float32]
    var pdps: seq[float32]
    for g in group:
      kdps.add(dp_norm_kid[g.dps_i])
      pdps.add(dp_norm_parent[g.dps_i])
    var kid_med = median(kdps)
    var parent_med = median(pdps)

    #(hq_kid_hets: 0, kid_hets: 5, hq_parent_hets: 1, hq_parent_hom_alts: 2, supporting: 4, maybe: 16, supportingLowQ: 3, nearbyLowQ: 0, nearbyN: 0)
    var s = &"{site.chrom}\t{group[0].start}\t{group[group.high].start + 1}\t{group.len}\t{group[group.high].i - group[0].i + 1}\t{duo.kid.id}\t{duo.parent.id}\t{kid_med:.3f}\t{parent_med:.3f}\t"
    s &= &"{stats.hq_kid_hets}\t{stats.kid_hets}\t{stats.hq_parent_hets}\t{stats.hq_parent_hom_alts}"
    echo s

template min16*(i:int32):uint16 =
  max(0, min(uint16.high.int32, i)).uint16

proc toDP(AD: seq[int32]): seq[uint16] =
  result = newSeq[uint16](int(AD.len / 2))
  for i in 0..<result.len:
    result[i] = min16(AD[2*i]) + min16(AD[2*i+1])

include ./utils

proc main*(dropfirst:bool=false) =
  ## TODO: find homozygous deletions with ./.
  var p = newParser("slivar duodel"):
    help("""find denovo structural deletions in parent-child duos using non-transmission of SNPs
    see: https://github.com/brentp/slivar/wiki/finding-deletions-in-parent-child-duos
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
  if opts.help:
    echo p.help
    quit 0

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

  var sites = newSeq[seq[Site]](duos.len)

  var exclude: TableRef[string, Lapper[region]]
  if opts.exclude != "":
    exclude = read_bed(opts.exclude)

  var GQ: seq[int32]
  var AD: seq[int32]
  var x: seq[int32]
  var MQ: seq[float32]
  var DPs = newSeqOfCap[seq[uint16]](16384)
  var last_rid = -1
  var min_sites = parseInt(opts.min_sites)
  var min_size = parseInt(opts.min_size)
  echo "#chrom\tstart\tstop\tn_supporting_sites\tn_total_sites\tkid_id\tparent_id\tkid_median_dp\tparent_median_dp\thq_kid_hets\tkid_hets\thq_parent_hets\thq_parent_hom_alts"
  for v in ivcf:

    if v.rid != last_rid:
      if last_rid != -1:
        var fdps = DPs.normalize
        for duo in duos:
          duo.cull(sites[duo.i], fdps, min_sites=min_sites, min_size=min_size)
          sites[duo.i].setLen(0)
        DPs.setLen(0)
      last_rid = v.rid

    if v.CHROM == "X" or v.CHROM == "chrX" or v.CHROM == "MT" or v.CHROM == "chrM" or v.CHROM == "chrMT": continue
    if v.FILTER notin ["PASS", "", "."]: continue
    if len(v.ALT) > 1: continue
    if len(v.REF) != 1: continue
    if len(v.ALT[0]) != 1: continue
    if exclude != nil and stripChr(v.CHROM) in exclude and 0 != exclude[stripChr(v.CHROM)].count(v.start.int, v.start.int + 1):
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
    if AD.len != 2 * v.n_samples:
      # some incorrectly decomposed vcfs have this.
      continue

    DPs.add(AD.toDP)

    var alts = v.format.genotypes(x).alts

    for duo in duos:
      var m = transmit(duo, GQ, AD, alts)
      if Transmit.Uninformative in m: continue
      if LowQual:
        m.incl(Transmit.LowQual)
      var c = Site(chrom: $v.CHROM, start: v.start, status: m, dps_i: DPs.high, i: sites[duo.i].len)
      c.ads[0] = min16(AD[2*duo.kid.i])
      c.ads[1] = min16(AD[2*duo.kid.i+1])
      c.ads[2] = min16(AD[2*duo.parent.i])
      c.ads[3] = min16(AD[2*duo.parent.i+1])
      sites[duo.i].add(c)

  var fdps = DPs.normalize
  for duo in duos:
    duo.cull(sites[duo.i], fdps, min_sites=min_sites, min_size=min_size)

when isMainModule:
  main()
