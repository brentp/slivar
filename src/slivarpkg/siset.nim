## This module contains an `intSet` object that is very fast for set operations for small positive integers
##
## Nim's builtin sets are very fast for `incl` and `excl` operations, but quite slow for iteration and
## cardinality checking. This can be preferable to the builtin sets for those use-cases
##
## .. code-block::
##   import siset
##
##   var si = initIntSet(128)
##   for i in 0..10: si.incl(i)
##   assert 5 in si
##   si.excl(5)
##   assert 5 notin si
##   assert si.card == 10
##   so = initIntSet(128)
##   assert (si * so).card == 0
##
##   so.incl(22)
##   assert (si + so).card == 11
##   assert 22 in (si + so)
##   assert $so == "{22}"
##   so.incl(2)
##   assert $so == "{2, 22}"
##
##
## this includes code modified from https://github.com/willf/bitset
## the (BSD-3) License for that is included in this directory.

import strutils
import bitops

const wordSize = int(64)
const log2WordSize = int(6)

type intSet* = object
  values: seq[int64]

proc bitsizeof(x:typedesc):int {.inline.} =
  sizeof(x)*8

proc initIntSet*(maxsize:int): intSet =
  return intSet(values: newSeq[int64](int(0.5 + maxsize/bitsizeof(int64))))

proc initIntSetUn(maxsize:int): intSet =
  return intSet(values: newSeqUninitialized[int64](int(0.5 + maxsize/bitsizeof(int64))))

template incl*(b:intSet, s:SomeOrdinal) =
  assert s.int < b.values.len * bitsizeof(int64) and s.int >= 0
  b.values[s shr log2WordSize] = b.values[s shr log2WordSize] or (1.int shl (s.int and (wordSize - 1)))

template `+`*(b:intSet, s:SomeOrdinal) =
  b.incl(s)

template excl*(b:intSet, s:SomeOrdinal) =
  assert s.int < b.values.len * bitsizeof(int64) and s.int >= 0
  b.values[s shr log2WordSize] = b.values[s shr log2WordSize] and not (1.int shl (s.int and (wordSize - 1)))

template `-`*(b:intSet, s:SomeOrdinal) =
  b.excl(s)

template `==`*(a:intSet, b:intSet): bool =
  a.values == b.values

template contains*(b:intSet, s:SomeOrdinal): bool =
  (b.values[s shr log2WordSize] and (1 shl (s and (wordSize-1)))) != 0

template `[]`*(b:intSet, s:SomeInteger): bool =
  b.contains(s)

proc card*(b:intSet): int {.inline.} =
  for i in 0..b.values.high:
    result.inc(countSetbits(b.values[i]))

proc `*`*(a:intSet, b:intSet): intSet {.inline.} =
  ## set intersection. note that this is fairly slow due to heap allocation
  ## of a new Seq.
  result = initIntSetUn(max(a.values.len * bitsizeof(int64), b.values.len * bitsizeof(int64)))
  for i in 0..min(a.values.high, b.values.high):
    result.values[i] = a.values[i] and b.values[i]

proc `*=`*(a:var intSet, b:intSet): intSet {.inline.} =
  for i in 0..max(a.values.high, b.values.high):
    a.values[i] = a.values[i] and b.values[i]

iterator items*(b:intSet): int =
  var i = 0
  var x = 0
  while x < b.values.len:
    while x < b.values.len and b.values[x] == 0:
      i += 8
      x += 1
    var w = b.values[x] shr (i and (wordSize - 1))
    if w != 0:
      i += countTrailingZeroBits(w)
      yield i
      i += 1
    else:
      i += 8
    x = i shr log2WordSize

proc `$`*(b:intSet): string =
  var s = newSeqOfCap[string](16)
  for v in b:
    s.add($v)
  return '{' & join(s, ", ") & '}'

when isMainModule:
  import unittest
  suite "intSet":
    test "incl excl card contains":
      var b = initIntSet(100)
      b.incl(10)
      check 10 in b
      check 11 notin b
      check b.card == 1
      doAssert b.card == 1

      b.excl(10)
      check 10 notin b
      check 11 notin b
      check b.card == 0
      for i in 0..<60:
        b.incl(i)
      check b.card == 60

    test "card all":
      var b = initIntSet(100)
      for i in 0..99:
        b.incl(i)
      check b.card == 100


    test "+ - *":
      var a = initIntSet(100)
      var b = initIntSet(100)
      for i in 0..<100:
        if i mod 2 == 0:
          a.incl(i)
        else:
          b.incl(i)
      check (a * b).card == 0
      check a.card == 50
      check b.card == 50
      #check (a - b).card == 50
      #check (a - a).card == 0
      #check (a + a).card == 50
      #check (a + b).card == 100

    test "string method":
      var a = initIntSet(100)
      a.incl(22)
      a.incl(33)
      a.incl(44)
      check $a == "{22, 33, 44}"
      a.excl(21)
      check $a == "{22, 33, 44}"
      a.excl(22)
      check $a == "{33, 44}"

  import times

  template iteration_speed(): int {.dirty.} =
    var sum :int = 0
    for i in 0..N:
      sum = 0
      for k in iset:
        sum += k.int
      doAssert iset.card == 40
      doAssert sum == 890
    sum

  template incl_speed(): int {.dirty.} =
    ## incl then iter, then exclude.
    for i in 0..N:
      for v in [5'u8, 10, 15, 20, 25, 30]:
        iset.incl(v)
      if i == 0: doAssert iset.card == 6
      for v in iset:
        iset.excl(v)
    doAssert iset.card == 0
    iset.card

  template fill_set() {.dirty.} =
    for i in 0'u8..50:
      if i < 30'u8 or i > 40'u8:
        iset.incl(i)

  template intersect_speed(): int {.dirty.} =
    for i in 0..N:
      var t = (iset * jset)
      doAssert t.card == 9
      #var s = 0
      #doAssert t[41]
      #for v in t:
      #  s += v.int
      #  break
      #doAssert s == 41

    (iset * jset).card

  var N = 2000000
  when defined(release):
    N *= 10

  proc siset_iteration(): int =
    var iset = initIntSet(64)
    fill_set()
    iteration_speed()

  proc systemset_iteration(): int =
    var iset : set[uint8]
    fill_set()
    iteration_speed()

  proc siset_incl(): int =
    var iset = initIntSet(64)
    incl_speed()

  proc systemset_incl(): int =
    var iset : set[uint8]
    incl_speed()

  proc siset_intersect(): int =
    var iset = initIntSet(64)
    fill_set()
    var jset = initIntSet(64)
    for i in 41..49: jset.incl(i)
    intersect_speed()

  proc systemset_intersect(): int =
    var iset : set[uint8]
    fill_set()
    var jset : set[uint8]
    for i in 41..49: jset.incl(i.uint8)
    intersect_speed()

  var t = cpuTime()

  t = cpuTime()
  echo "siset iteration:", siset_iteration(), " ", (cpuTime() - t)
  t = cpuTime()
  echo "system set iteration:", systemset_iteration(), " ", (cpuTime() - t)

  t = cpuTime()
  echo "siset incl:", siset_incl(), " ", (cpuTime() - t)
  t = cpuTime()
  echo "system set incl:", systemset_incl(), " ", (cpuTime() - t)

  t = cpuTime()
  echo "siset intersect:", siset_intersect(), " ", (cpuTime() - t)
  t = cpuTime()
  echo "system set intersect:", systemset_intersect(), " ", (cpuTime() - t)
