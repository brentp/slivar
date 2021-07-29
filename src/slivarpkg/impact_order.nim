import tables
import strutils
import os

const order_x = staticRead("./default-order.txt")

proc adjustOrder*(order: string): TableRef[string, int] =
  result = newTable[string, int](16)
  for o in order.strip().split("\n"):
    if o[0] == '#': continue
    var n = o.toLowerAscii.strip()
    # keep IMPACTFUL_CUTOFF as-is, all other impacts are lower-cased
    if n[0] == 'i' and n in ["impact_cutoff", "impactful_cutoff"]: n = "IMPACTFUL_CUTOFF"
    if o.len == 0: continue
    if n.endsWith("_variant"):
      n = n[0..n.high - 8]
    result[n] = result.len

  if "IMPACTFUL_CUTOFF" notin result:
    stderr.write_line "[slivar] warning! no IMPACTFUL_CUTOFF specified. INFO.impactful will not be useful"

var default_order* = adjustOrder(order_x)
if getEnv("SLIVAR_IMPACTFUL_ORDER") != "":
  if not fileExists(getEnv("SLIVAR_IMPACTFUL_ORDER")):
    raise newException(IOError, "[slivar] couldn't open file at:" & getEnv("SLIVAR_IMPACTFUL_ORDER") & " specified by env var 'SLIVAR_IMPACTFUL_ORDER'")
  default_order = adjustOrder(getEnv("SLIVAR_IMPACTFUL_ORDER").readFile)
