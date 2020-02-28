import hts
import tables
import algorithm
import strutils
import lapper
export lapper

type
  region* = object
    chrom*: string
    start*: int
    stop*: int
    #name: string

template stripChr*[T:string|cstring](s:T): string =
  if s.len > 3 and ($s).startswith("chr"): ($s)[3..<s.len] else: $s

proc start*(r: region): int {.inline.} = r.start
proc stop*(r: region): int {.inline.} = r.stop

proc bed_line_to_region(line: string): region =
  var
    cse = line.strip().split('\t', 5)
  if len(cse) < 3:
    stderr.write_line("[slivar] skipping bad bed line:", line.strip())
    return
  result.chrom = stripChr(cse[0])
  result.start = parse_int(cse[1])
  result.stop = parse_int(cse[2])
   #if len(cse) > 3:
   #  reg.name = cse[3]

proc read_bed*(bed: string): TableRef[string, Lapper[region]] =
  var bed_regions = newTable[string, seq[region]]()
  var kstr = kstring_t(l:0, m: 0, s: nil)
  var hf = hts_open(cstring(bed), "r")
  while hts_getline(hf, cint(10), addr kstr) > 0:
    if kstr.s[0] == 't' and ($kstr.s).startswith("track "):
      continue
    if kstr.s[0] == '#':
      continue
    var v = bed_line_to_region($kstr.s)
    if v.chrom.len == 0: continue
    discard bed_regions.hasKeyOrPut(v.chrom, new_seq[region]())
    bed_regions[v.chrom].add(v)

  # since it is read into mem, can also well sort.
  result = newTable[string, Lapper[region]]()
  for chrom, ivs in bed_regions.mpairs:
    result[chrom] = lapify(ivs)

  hts.free(kstr.s)
  return result

