import algorithm

type encoded = object
  filter {.bitsize: 2.} : uint64 ## whether a non-pass filter is set on the variant.
  enc {.bitsize: 26.} : uint64
  alen {.bitsize: 4.} : uint64
  rlen {.bitsize: 4.} : uint64
  pos {.bitsize: 28.} : uint64

const MaxCombinedLen* = 13

type pfra* = object
  # position,ref,alt
  position*: uint32
  filter*: bool
  reference*: string
  alternate*: string

proc cmp_pfra*(a, b:pfra): int =
  if a.position != b.position:
    return cmp[uint32](a.position, b.position)

  if a.reference != b.reference:
    if a.reference.len != b.reference.len:
      return cmp[int](a.reference.len, b.reference.len)
    return cmp(a.reference, b.reference)

  if a.alternate.len != b.alternate.len:
    return cmp[int](a.alternate.len, b.alternate.len)
  return cmp(a.alternate, b.alternate)

const lookup:array[128, uint8] = [3'u8,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,0,3,1,3,3,3,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3]

const rlookup:array[4, char] = ['A', 'C', 'G', 'T']

proc encode*(pos: uint32, ref_allele: string, alt_allele: string, filter:bool): uint64 {.inline.} =
  var p = encoded(pos:pos, rlen:ref_allele.len.uint64, alen:alt_allele.len.uint64, filter:filter.uint64)
  if ref_allele.len + alt_allele.len > MaxCombinedLen:
    p.rlen = 0; p.alen = 0
    return cast[uint64](p)

  var ra = 0'u32

  for k in 0..<ref_allele.len:
    ra = ra * 4
    ra += lookup[ref_allele[k].int]

  for k in 0..<alt_allele.len:
    ra = ra * 4
    ra += lookup[alt_allele[k].int]

  p.enc = ra.uint64
  return cast[uint64](p)

template encode*(p:pfra): uint64 =
  encode(p.position, p.reference, p.alternate, p.filter)

template is_long_allele(p:pfra): bool =
  (p.reference.len == 0 and p.alternate.len == 0) or (p.reference.len + p.alternate.len > MaxCombinedLen)

proc match(q:pfra, t:pfra): bool {.inline.} =
  if q.position != t.position: return false
  if q.is_long_allele and t.is_long_allele: return true
  return q.reference == t.reference and q.alternate == t.alternate

proc decode*(e:uint64): pfra {.inline.} =
  var tmp = cast[encoded](e)
  result = pfra(position:tmp.pos.uint32, filter: tmp.filter.bool)
  if tmp.rlen == 0:
    return

  var e = tmp.enc
  result.alternate = newString(tmp.alen)
  for i in 0..<result.alternate.len:
    let index = e and 3
    e = e shr 2
    result.alternate[result.alternate.high - i] = rlookup[index]

  result.reference = newString(tmp.rlen)
  for i in 0..<result.reference.len:
    let index = e and 3
    e = e shr 2
    result.reference[result.reference.high - i] = rlookup[index]


proc ilowerBound*(a: seq[uint64], key:uint64): int =
  # the lower bound in stdlib uses a closure for cmp. this uses < directly
  result = a.low
  var count = a.high - a.low + 1
  var step, pos: int
  while count != 0:
    step = count shr 1
    pos = result + step
    if a[pos] < key:
      result = pos + 1
      count -= step + 1
    else:
      count = step

proc find*(coded:seq[uint64], q:pfra): int {.inline.} =
  var e = encode(q)
  var i = ilowerBound(coded, e)
  if i == -1: return -1
  if i > coded.high:
    i = coded.high

  if cast[encoded](coded[i]).pos.uint32 != q.position:
    return -1
  if coded[i] == e:
    return i

  var te = coded[i].decode
  while i < coded.len:
    if q.match(te):
      return i
    if i == coded.high: return - 1
    i += 1
    if cast[encoded](coded[i]).pos.uint32 != q.position:
      return -1
    te = coded[i].decode

  return - 1

when isMainModule:
  import unittest
  import times
  import random

  when defined(release):
    doAssert encoded().sizeof == 8

  suite "pra test suite":

    test "that long ref alt works":
      check encode(1232434'u32, "AAAAAAAAAAAAAAAAA", "TTTTTTTTTTTTTTTTTTT", false) > 0'u64

    test "that flag with long ref alt gets set":
      var v = encode(1232434'u32, "AAAAAAAAAAAAAAAAA", "TTTTTTTTTTTTTTTTTTT", false)
      check (cast[encoded](v)).pos == 1232434'u64
      check (cast[encoded](v)).rlen == 0
      check (cast[encoded](v)).filter == 0
      v = encode(1232434'u32, "AAAAAAAAAAAAAAAAA", "TTTTTTTTTTTTTTTTTTT", true)
      check (cast[encoded](v)).filter == 1

      var d = v.decode
      check d.reference == ""
      check d.filter == true
      check d.alternate == ""
      check d.position == 1232434

    test "that encode/decode roundtrip works":
      var v = encode(122434'u32, "AAAAAA", "TTTT", true)
      echo v
      var d = v.decode
      check d.filter
      echo d

    test "ordering":
      var aL = encode(1232434'u32, "AAAAAAAAAAAAAAAAA", "TTTTTTTTTTTTTTTTTTT", false)
      var bL = encode(1232434'u32, "AAAAAAAAAAAAAAAAA", "TTTTTTTTTTTTTTTTTTT", false)

      var a = encode(1232434'u32, "A", "T", false)
      var b = encode(1232434'u32, "A", "T", false)

      var cL = encode(1232434'u32, "A", "TTTTTTTTTTTTTTTTTTT", false)

      check bL <= aL
      check b == a

      check cmp_pfra(a.decode, b.decode) == 0
      check cmp_pfra(aL.decode, bL.decode) == 0

      check cL < a

    test "find with all exact":

      var haystack = @[
        encode(123'u32, "A", "T", false),
        encode(123'u32, "A", "T", true),
        encode(123'u32, "A", "CTTT", false),
        encode(123'u32, "A", "T", false),
        encode(123'u32, "G", "T", false),
        encode(123'u32, "A", "G", false),
        encode(124'u32, "A", "G", false),
        encode(124'u32, "A", "T", false)]
      sort(haystack, cmp[uint64])

      for i, h in haystack:
        check haystack[haystack.find(h.decode)].decode == h.decode
    test "find with changed flag ":

      var haystack = @[
        encode(123'u32, "A", "T", false),
        encode(123'u32, "A", "T", true),
        encode(123'u32, "A", "CTTT", true),
        encode(123'u32, "A", "T", true),
        encode(123'u32, "G", "T", true),
        encode(123'u32, "A", "G", true),
        encode(124'u32, "A", "T", true)]
      sort(haystack, cmp[uint64])
      for h in haystack:
        var v = h.decode
        var found = haystack[haystack.find(v)].decode
        var exp = h.decode
        check found.position == exp.position
        check found.reference == exp.reference
        check found.alternate == exp.alternate

    test "find with large events":

      var haystack = @[
        encode(124'u32, "GGGCGGGGGGG", "TTTTTTTTTT", true),
        encode(124'u32, "GGGCGGGGGGG", "T", false)]
      sort(haystack, cmp[uint64])
      echo "["
      for e in haystack:
        echo e.decode, "->", e
      echo "]"

      var obsi = haystack.find(pfra(position: 124'u32, reference: "GGGCGGGGGGG", alternate:"TTTTTTTTTT"))
      var obs = haystack[obsi].decode
      check obs.reference == ""
      check obs.alternate == ""

      obsi = haystack.find(pfra(position: 124'u32, reference: "", alternate:""))
      obs = haystack[obsi].decode
      check obs.position == 124
      check obs.reference == ""
      check obs.alternate == ""

      obsi = haystack.find(pfra(position: 124'u32, reference: "GGGCGGGGG", alternate:"TT"))
      check obsi == -1

    test "with max size":

     for i in 0..100000:
       var alen = rand(1..<MaxCombinedLen)
       var rlen = MaxCombinedLen - alen
       check alen < MaxCombinedLen
       check alen + rlen == MaxCombinedLen

       var aref = ""
       for i in 0..<rlen:
         aref &= sample("ACTG")

       var aalt = ""
       for i in 0..<alen:
         aalt &= sample("ACTG")
       var p = rand(250_000_000).uint32
       var e = encode(p, aref, aalt, false)
       var d = e.decode
       check d.reference == aref
       check d.alternate == aalt
       check d.position == p

    test "with large and small same":
       var haystack = @[
         pfra(position: 874816, reference: "", alternate: "").encode,   # 60116897909992206
         pfra(position: 874816, reference: "", alternate: "").encode,   # 60116897957286235
         pfra(position: 874816, reference: "C", alternate: "T").encode # 60116902323683335
       ]
       check haystack.find(pfra(position: 874816, reference: "C", alternate: "T")) == 2

    var t = cpuTime()
    var n = 5_000_000
    when not defined(release):
        n = 100_000

    var refs = ["AA", "TT", "CC", "GG", "GGG", "G", "CCC", "CGCGCCG", "TA",]
    var alts = ["T", "C", "GGGGG", "AAA", "GGC", "TTCAT", "AAAG", "GGG"]
    for i in 0..n:
      var p = rand(250_000_000).uint32
      var aref = sample(refs)
      var aalt = sample(alts)
      var e = encode(p, aref, aalt, i mod 2 == 0)
      var d = e.decode
      doAssert d.position == p
      doAssert d.reference == aref
      doAssert d.alternate == aalt
      doAssert d.filter == (i mod 2 == 0)

    echo int(n.float64 / (cpuTime() - t)), " encode/decodes per second"


    var haystack = newSeq[uint64]()
    for i in 0..n:
      var p = rand(250_000_000).uint32
      var aref = sample(refs)
      var aalt = sample(alts)
      var e = encode(p, aref, aalt, false)
      haystack.add(e)

    sort(haystack)

    for h in haystack:
      check haystack[haystack.find(h.decode)].decode == h.decode
    t = cpuTime()
    for h in haystack:
      check haystack.find(h.decode) != -1
    echo int(n.float64 / (cpuTime() - t)), " finds per second"


    var needle = encode(1233333, "AAA", "T", false)
    t = cpuTime()
    for i in 0..haystack.high:
      check haystack.find(needle.decode) == -1
    echo int(n.float64 / (cpuTime() - t)), " non finds per second"
