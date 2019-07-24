
import bitops

iterator sub*[T: uint8|uint16](a: set[T], b: set[T]): T =

  when sizeof(T) == 2:
    let ac = cast[array[1024, uint64]](a)
    var bc = cast[array[1024, uint64]](b)
  elif sizeof(T) == 1:
    let ac = cast[array[4, uint64]](a)
    var bc = cast[array[4, uint64]](b)

  # modified from lemire: https://lemire.me/blog/2018/02/21/iterating-over-set-bits-quickly/
  for i, a64 in ac:
    if likely(a64 == bc[i]):
      continue
    var a64 = a64
    a64.clearMask(bc[i])

    while a64 != 0:
      let t = a64.bitAnd(0'u64 - a64)
      let r = a64.countTrailingZeroBits
      yield (i * 64 + r).T
      a64.flipMask(t)

when isMainModule:
  import unittest
  import sequtils

  suite "fast set subtraction":

    test "that subtraction works":

      var a: set[uint16]
      var b: set[uint16]

      for i in 0..20:
        a.incl(i.uint16)
        b.incl(i.uint16)

      #a.incl(21)
      for i in 22..50:
        a.incl(i.uint16)
        b.incl(i.uint16)



      a.incl(66)
      a.incl(61)
      a.incl(2198)

      #echo "66 added"
      check toSeq(a.sub(b)) == toSeq(a - b)




