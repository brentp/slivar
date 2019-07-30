import tables
import strutils

const order_x = staticRead("./default-order.txt")

proc adjustOrder*(order: string): TableRef[string, int] =
  result = newTable[string, int](16)
  for o in order.strip().split("\n"):
    if o[0] == '#': continue
    var n = o.toLowerAscii.strip()
    # keep IMPACT_CUTOFF as-is, all other impacts are lower-cased
    if n[0] == 'i' and n == "impact_cutoff": n = "IMPACT_CUTOFF"
    if o.len == 0: continue
    if n.endsWith("_variant"):
      n = n[0..n.high - 8]
    result[n] = result.len

  if "IMPACT_CUTOFF" notin result:
    stderr.write_line "[slivar] warning! no IMPACT_CUTOFF specified. INFO.impactful will not be useful"

var default_order* = adjustOrder(order_x)
