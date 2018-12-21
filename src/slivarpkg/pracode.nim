import hashes
import algorithm

type pra = object
  flag {.bitsize: 4} : uint64 # 56
  ra {.bitsize: 24} : uint64 # 52
  rlen {.bitsize: 4} : uint64 # 60
  alen {.bitsize: 4} : uint64 # 64

  pos {.bitsize: 28} : uint64 # 28

type pfra* = object
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

const lookup:array[256, uint8] = [3'u8,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,0,3,1,3,3,3,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3]

const rlookup:array[4, char] = ['A', 'C', 'G', 'T']

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

  var ra = 0'u32
  for k in 0..<ref_allele.len:
    ra = ra * 4
    var c = lookup[ref_allele[k].int]
    ra += c

  for k in 0..<alt_allele.len:
    ra = ra * 4
    var c = lookup[alt_allele[k].int]
    ra += c
  p.ra = ra.uint64
  return cast[uint64](p)

template encode*(p:pfra): uint64 =
  encode(p.position, p.reference, p.alternate, p.flag)

proc long_or_short(p:pfra): bool {.inline.} =
  return (p.reference.len == 0 and p.alternate.len == 0) or (p.reference.len + p.alternate.len > 12)

proc match(q:pfra, t:pfra): bool {.inline.} =
  if q.position != t.position: return false
  if long_or_short(q) and long_or_short(t): return true
  return q.reference == t.reference and q.alternate == t.alternate

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

proc find*(encoded:seq[uint64], q:pfra): int {.inline.} =
  var e = encode(q)
  var i = lowerBound[uint64](encoded, e)
  if i > encoded.high:
    i = encoded.high
  while i > 0:
    if encoded[i] == e:
      return i
    # have to go past so that we can get the lowest of this position.
    if (cast[pra](encoded[i])).pos.uint32 >= q.position:
      i -= 1
      continue
    else: break

  while i < encoded.len and cast[pra](encoded[i]).pos.uint32 < q.position:
    i += 1
  if cast[pra](encoded[i]).pos.uint32 > q.position:
    return -1

  var te = encoded[i].decode
  while i < encoded.len:
    if te.position > q.position: return -1
    if te.position < q.position:
      if i == encoded.high: return -1
      i += 1
      te = encoded[i].decode
      continue
    if q.match(te):
      return i

    if i == encoded.high: return - 1
    i += 1
    te = encoded[i].decode

  return - 1

when isMainModule:
  import unittest
  import times
  import random

  when defined(release):
    doAssert pra().sizeof == 8

  suite "pra test suite":

    test "that long ref alt works":
      check encode(1232434'u32, "AAAAAAAAAAAAAAAAA", "TTTTTTTTTTTTTTTTTTT") > 0'u64

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

    test "ordering":
      var aL = encode(1232434'u32, "AAAAAAAAAAAAAAAAA", "TTTTTTTTTTTTTTTTTTT", 3)
      var bL = encode(1232434'u32, "AAAAAAAAAAAAAAAAA", "TTTTTTTTTTTTTTTTTTT", 0)

      var a = encode(1232434'u32, "A", "T", 3)
      var b = encode(1232434'u32, "A", "T", 0)

      var cL = encode(1232434'u32, "A", "TTTTTTTTTTTTTTTTTTT", 0)

      check bL < aL
      check b < a

      check cmp_pfra(a.decode, b.decode) == 1
      check cmp_pfra(aL.decode, bL.decode) == 1

      check cL < a

    test "find with all exact":

      var haystack = @[
        encode(123'u32, "A", "T", 0'u8),
        encode(123'u32, "A", "T", 2'u8),
        encode(123'u32, "A", "CTTT", 0'u8),
        encode(123'u32, "A", "T", 3'u8),
        encode(123'u32, "G", "T", 3'u8),
        encode(123'u32, "A", "G", 0'u8),
        encode(124'u32, "A", "G", 0'u8),
        encode(124'u32, "A", "T", 0'u8)]
      sort(haystack, cmp[uint64])

      for i, h in haystack:
        check haystack[haystack.find(h.decode)].decode == h.decode
    test "find with changed flag ":

      var haystack = @[
        encode(123'u32, "A", "T", 0'u8),
        encode(123'u32, "A", "T", 2'u8),
        encode(123'u32, "A", "CTTT", 0'u8),
        encode(123'u32, "A", "T", 3'u8),
        encode(123'u32, "G", "T", 3'u8),
        encode(123'u32, "A", "G", 0'u8),
        encode(124'u32, "A", "G", 0'u8),
        encode(124'u32, "A", "T", 0'u8)]
      sort(haystack, cmp[uint64])
      for h in haystack:
        var v = h.decode
        v.flag = 6
        var found = haystack[haystack.find(v)].decode
        var exp = h.decode
        check found.position == exp.position
        check found.reference == exp.reference
        check found.alternate == exp.alternate

    test "find with large events":

      var haystack = @[
        encode(124'u32, "GGGCGGGGGGG", "TTTTTTTTTT", 0'u8),
        encode(124'u32, "GGGCGGGGGGG", "T", 0'u8)]
      sort(haystack, cmp[uint64])
      echo "["
      for e in haystack:
        echo e.decode, "->", e
      echo "]"

      var obsi = haystack.find(pfra(position: 124'u32, reference: "GGGCGGGGGGG", alternate:"TTTTTTTTTT", flag: 0'u8))
      var obs = haystack[obsi].decode
      check obs.reference == ""
      check obs.alternate == ""

      obsi = haystack.find(pfra(position: 124'u32, reference: "", alternate:"", flag: 0'u8))
      obs = haystack[obsi].decode
      check obs.position == 124
      check obs.reference == ""
      check obs.alternate == ""

      obsi = haystack.find(pfra(position: 124'u32, reference: "GGGCGGGGG", alternate:"TT", flag: 0'u8))
      check haystack[obsi].decode.position == 0


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


    var haystack = newSeq[uint64]()
    for i in 0..n:
      var p = rand(250_000_000).uint32
      var aref = rand(refs)
      var aalt = rand(alts)
      var e = encode(p, aref, aalt)
      haystack.add(e)

    sort(haystack)

    for h in haystack:
      check haystack[haystack.find(h.decode)].decode == h.decode
    t = cpuTime()
    for h in haystack:
      check haystack.find(h.decode) != -1
    echo int(n.float64 / (cpuTime() - t)), " finds per second"


    var needle = encode(1233333, "AAA", "T")
    t = cpuTime()
    for i in 0..haystack.high:
      check haystack.find(needle.decode) == -1
    echo int(n.float64 / (cpuTime() - t)), " non finds per second"
