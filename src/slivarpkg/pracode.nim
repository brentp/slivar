import hashes

type pra = object
  pos {.bitsize: 28} : uint64 # 28
  flag {.bitsize: 4} : uint64 # 32
  rlen {.bitsize: 4} : uint64 # 36
  alen {.bitsize: 4} : uint64 # 40
  ra {.bitsize: 24} : uint64 # 64

type pfra* = ref object
  # position,flag,ref,alt
  position*: uint32
  flag*: uint8
  reference*: string
  alternate*: string

proc cmp_pfra*(a, b:pfra): int =
  if a.position != b.position:
    return cmp[uint32](a.position, b.position)
  if a.reference != b.reference:
    return cmp(a.reference, b.reference)
  if a.alternate != b.alternate:
    return cmp(a.alternate, b.alternate)
  return cmp(a.flag, b.flag)

proc `$`*(p:pfra): string {.inline.} =
  return $p[]

const lookup:array[256, uint8] = [0'u8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,1,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

const rlookup:array[4, char] = ['T', 'C', 'A', 'G']

proc hash(a: string, b:string): uint64 {.inline.} =
  ## Computes a Hash from `x`.
  var h: Hash = 0
  # Iterate over parts of `x`.
  for xAtom in a:
    # Mix the atom with the partial hash.
    h = h !& xAtom.int
  # Finish the hash.
  for xAtom in b:
    h = h !& xAtom.int
  result = (!$h).uint32.uint64

proc encode*(pos: uint32, ref_allele: string, alt_allele: string, flag:uint8=0): uint64 {.inline.} =
  var p = pra(pos:pos, flag: flag)
  if ref_allele.len + alt_allele.len > 12:
    p.ra = hash(ref_allele, alt_allele)
    return cast[uint64](p)

  p.rlen = ref_allele.len.uint64
  p.alen = alt_allele.len.uint64

  var ra = 0'u64
  for k in 0..<ref_allele.len:
    ra = ra * 4
    var c = lookup[ref_allele[k].int]
    ra += c

  for k in 0..<alt_allele.len:
    ra = ra * 4
    var c = lookup[alt_allele[k].int]
    ra += c
  p.ra = ra
  return cast[uint64](p)

proc decode*(e:uint64): pfra {.inline.} =
  var tmp = cast[pra](e)
  result = pfra(position:tmp.pos.uint32, flag: tmp.flag.uint8)
  if tmp.rlen == 0:
    return

  var e = tmp.ra

  result.alternate = newString(tmp.alen)
  for i in 0..<result.alternate.len:
    var index = e and 3
    e = e shr 2
    result.alternate[result.alternate.high - i] = rlookup[index]

  result.reference = newString(tmp.rlen)
  for i in 0..<result.reference.len:
    var index = e and 3
    e = e shr 2
    result.reference[result.reference.high - i] = rlookup[index]

when isMainModule:
  import unittest
  import times
  import random

  when defined(release):
      doAssert pra().sizeof == 8

  suite "pra test suite":

    test "that long ref alt works":
      check encode(1232434'u32, "AAAAAAAAAAAAAAAAA", "TTTTTTTTTTTTTTTTTTT") == 15766199596378934834'u64

    test "that flag with long ref alt gets set":
      var v = encode(1232434'u32, "AAAAAAAAAAAAAAAAA", "TTTTTTTTTTTTTTTTTTT", 3)
      check (cast[pra](v)).flag == 3
      v = encode(1232434'u32, "AAAAAAAAAAAAAAAAA", "TTTTTTTTTTTTTTTTTTT", 5)
      check (cast[pra](v)).flag == 5
      check (cast[pra](v)).pos == 1232434'u64
      check (cast[pra](v)).rlen == 0

      var d = v.decode
      check d.reference == ""
      check d.alternate == ""
      check d.position == 1232434

    test "that encode/decode roundtrip works":
      var v = encode(122434'u32, "AAAAAA", "TTTT", 3)
      echo v
      var d = v.decode
      echo d

    var t = cpuTime()
    var n = 5_000_000
    when not defined(release):
        n = 100_000

    var refs = ["AA", "TT", "CC", "GG", "GGG", "G", "CCC", "CGCGCCG", "TA",]
    var alts = ["T", "C", "GGGGG", "AAA", "GGC", "TTCAT", "AAAG", "GGG"]
    for i in 0..n:
      var p = rand(250_000_000).uint32
      var aref = rand(refs)
      var aalt = rand(alts)
      var e = encode(p, aref, aalt)
      var d = e.decode
      doAssert d.position == p
      doAssert d.reference == aref
      doAssert d.alternate == aalt

    echo int(n.float64 / (cpuTime() - t)), " encode/decodes per second"
