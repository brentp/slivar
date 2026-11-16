import hts/vcf
import hts/private/hts_concat
import hts/files
import math
import pedfile
import jsc
import json
import os
import ./impact_order
import ./tsv
import ./gnotate
import ./groups
import tables
import strformat
import strutils

proc clean(s:string): string =
  return s.replace("\"", "").replace("'", "")

type CompiledExpression* = object
  name*: string
  expr*: string

type NamedExpression* = object
  name*: string
  expr*: string

proc assertNewLabel(label: string, previous: seq[NamedExpression]) =
  for pe in previous:
    if pe.name == label:
      quit &"[slivar] ERROR duplicate label '{pe.name}'. Use unique labels"

proc getNamedExpressions*(ovcf:VCF, expressions:seq[string], invcf:string, asfloat:bool, previousExpressions: varargs[seq[NamedExpression]]): seq[NamedExpression] =
  ## parse (split on :) the expressions, return a sequence, and update the vcf header with a reasonable description.
  result = newSeq[NamedExpression]()
  for e in expressions:
    var t = e.split(seps={':'}, maxsplit=1)
    if t.len != 2:
      quit "must specify name:expression pairs. got:" & e

    assertNewLabel(t[0], result)
    for p in previousExpressions:
      assertNewLabel(t[0], p)

    result.add(NamedExpression(name:t[0], expr:t[1]))
    if asfloat:
      if ovcf.header.add_info(t[0], "1", "Float", &"added by slivar with expression: '{t[1].clean}'") != Status.OK:
        quit "error adding field to header"
    else:
      if ovcf.header.add_info(t[0], ".", "String", &"added by slivar with expression: '{t[1].clean}' from {invcf}") != Status.OK:
        quit "error adding field to header"

proc isDigit(s:string): bool {.inline.} =
  for c in s:
    if not c.isdigit: return false
  true

proc looks_like_region_file(f:string): bool =
  if not f.fileExists: return false
  if '-' in f and ':' in f: return false
  var fh:HTSFile
  if not open(fh, f):
    stderr.write_line &"[slivar] tried '{f}' as a region file but couldn't open. Trying as an actual region"
    return false
  defer:
    fh.close()
  for l in fh.lines:
    if l[0] == '#' or l.strip().len == 0: continue
    var toks = l.strip().split("\t")
    if toks.len >= 3 and toks[1].isdigit and toks[2].isdigit: return true
    stderr.write_line &"[slivar] tried '{f}' as a region file but it did not have proper format. Trying as an actual region"
    return false

iterator variants*(vcf:VCF, region:string): Variant =
  ## iterator over region or just the variants.
  if region == "" or region == "nil":
    for v in vcf: yield v
  elif fileExists(region) and looks_like_region_file(region):
    ## must be in bed format.
    for l in region.hts_lines(threads=0):
      if l[0] == '#' or l.strip().len == 0: continue
      var toks = l.strip().split(seps={'\t'})
      for v in vcf.query(&"{toks[0]}:{parseInt(toks[1]) + 1}-{toks[2]}"):
        yield v
  else:
    for v in vcf.query(region): yield v



type ISample* = ref object
  ped_sample*: pedfile.Sample

type exResult* = tuple[name:string, sampleList:seq[string], val:float32]
  ## a match yielded by evaluate: the expression name and the matching sample ids.
  ## val is float32.low for a family that matched without affected samples.

type idpair* = tuple[name:string, info:bool, format:bool]

type Unit = object
  ## a trio, a group row or a family
  ids: seq[string]   ## reported when an expression matches (kid, first sample, affected samples)
  family_id: string  ## families only

# ---------------------------------------------------------------------------
# shared memory between nim and javascript.
#
# each Buf is a C-allocated array that javascript sees as a typed array
# ($sl.<name>). When a buffer must grow, a new allocation is made and the typed
# array in javascript is replaced. The old memory is freed by JavaScriptCore
# once the old typed array is garbage collected (so a stale view is never
# left pointing at freed memory).

proc c_malloc(size: csize_t): pointer {.importc: "malloc", header: "<stdlib.h>".}
proc c_free(p: pointer) {.importc: "free", header: "<stdlib.h>".}

proc jsFree(bytes: pointer, ctx: pointer) {.cdecl.} = c_free(bytes)

type Buf[T] = object
  data: ptr UncheckedArray[T]
  cap: int

type Bridge = object
  ctx: JSContext
  sl: JSObjectRef  # the $sl object

proc bind_buf[T](b: Bridge, buf: var Buf[T], name: string, cap: int, index: int = -1) =
  ## (re-)allocate buf with capacity cap, preserving content, and expose it as $sl[name] (or $sl[name][index]).
  let ncap = max(cap, 1)
  let p = cast[ptr UncheckedArray[T]](c_malloc(csize_t(ncap * sizeof(T))))
  zeroMem(p, ncap * sizeof(T))
  if buf.data != nil and buf.cap > 0:
    copyMem(p, buf.data, min(buf.cap, ncap) * sizeof(T))
  var ex: JSValueRef
  let arr = JSObjectMakeTypedArrayWithBytesNoCopy(b.ctx.ctx, typedArrayKind(T), p, csize_t(ncap * sizeof(T)), jsFree, nil, ex.addr)
  b.ctx.ctx.check(ex, "creating typed array " & name)
  if index < 0:
    b.ctx.ctx.setProp(b.sl, name, JSValueRef(arr))
  else:
    let holder = b.ctx.ctx.getProp(b.sl, name).toObject(b.ctx.ctx)
    b.ctx.ctx[holder, index] = JSValueRef(arr)
  buf.data = p
  buf.cap = ncap

template ensure[T](b: Bridge, buf: var Buf[T], name: string, n: int, index: int = -1) =
  if unlikely(n > buf.cap):
    b.bind_buf(buf, name, max(n, 2 * buf.cap), index)

template `[]`[T](b: Buf[T], i: int): T = b.data[i]
template `[]=`[T](b: Buf[T], i: int, v: T) = b.data[i] = v
template `[]`[T](b: Buf[T], e: enum): T = b.data[ord(e)]
template `[]=`[T](b: Buf[T], e: enum, v: T) = b.data[ord(e)] = v

# ---------------------------------------------------------------------------
# layout of the shared buffers. This is the only definition: the matching
# javascript constants are generated from it (jsLayout) and inserted into
# runtime.js at compile time.
type
  # enums with string values are visible in javascript under those names
  # (with a prefix, e.g. META_GENERATION).
  MetaSlot = enum
    msGeneration = "GENERATION"  ## incremented for each variant
    msAlts = "ALTS"              ## ALTS (and so alts, het, ...) has been set
    msAb = "AB"                  ## AB (computed from AD) is set for this variant
    msCalculated = "CALCULATED"  ## the variant fields calculated from genotypes are set

  VariantNumber = enum
    ## slots in VNUM; the string values are the javascript attribute names.
    vnStart = "start", vnStop = "stop", vnPOS = "POS", vnQUAL = "QUAL",
    vnMultiallelic = "is_multiallelic",
    # calculated from the genotypes (msCalculated):
    vnAaf = "aaf", vnHwe = "hwe_score", vnCallRate = "call_rate", vnNumHomRef = "num_hom_ref",
    vnNumHet = "num_het", vnNumHomAlt = "num_hom_alt", vnNumUnknown = "num_unknown"

  VariantString = enum
    ## the first string slots of each variant
    vsCHROM = "CHROM", vsREF = "REF", vsALT = "ALT", vsID = "ID", vsFILTER = "FILTER"

  InfoKind = enum
    ## INFO_T value for each INFO field
    ikAbsent = "ABSENT", ikNumber = "NUMBER", ikString = "STRING", ikFlag = "FLAG", ikArray = "ARRAY"

  UnitKind = enum
    ## what expressions are evaluated on, in evaluation order
    ukTrio = "trio", ukGroup = "group", ukFamily = "family"

  InfoResult = enum
    ## returned by sl.info()
    irFail = "FAIL"               ## the --info expression was false
    irPass = "PASS"               ## it passed (or there was none)
    irError = "ERROR"             ## it threw; counted as a pass

const FmtString = -1'i32  ## FMT_N value for a string FORMAT field

proc jsLayout(): string {.compileTime.} =
  proc names(E: typedesc[enum]): string =
    var parts: seq[string]
    for e in E: parts.add "\"" & $e & "\""
    "[" & parts.join(", ") & "]"
  result = "  // generated from the enums in evaluator.nim\n"
  for e in MetaSlot: result.add &"  var META_{e} = {ord(e)};\n"
  for e in InfoKind: result.add &"  var INFO_{e} = {ord(e)};\n"
  for e in InfoResult: result.add &"  var RESULT_{e} = {ord(e)};\n"
  result.add &"  var VARIANT_NUMBERS = {names(VariantNumber)};\n"
  result.add &"  var VARIANT_STRINGS = {names(VariantString)};\n"
  result.add &"  var FIRST_CALCULATED = {ord(vnAaf)}, MULTIALLELIC = {ord(vnMultiallelic)};\n"
  result.add &"  var UNIT_KINDS = {names(UnitKind)};\n"
  result.add &"  var FMT_STRING = {FmtString};\n"

type FmtBuf = object
  ## the values of a numeric FORMAT field. Only one of these is ever used for a
  ## given id: htslib requires the type to match the header (Integer or Float).
  ints: Buf[int32]
  floats: Buf[float32]

type Evaluator* = ref object
  b: Bridge
  samples*: seq[Isample]
  gnos*: seq[Gnotater]

  # vcf header info
  field_names: seq[idpair]       # header id -> name and whether it's an INFO / FORMAT field
  gene_fields: seq[GeneIndexes]  # CSQ/ANN/BCSQ fields used to compute impactful, genic, ...
  impactful_id: int              # header ids of the fields slivar adds (-1 if absent)
  genic_id: int
  highest_impact_id: int

  # settings from the environment
  allow_fmt_strings: bool        # SLIVAR_FORMAT_STRINGS
  info_white_list: seq[string]   # SLIVAR_INFO_WHITELIST
  format_white_list: seq[string] # SLIVAR_FORMAT_WHITELIST

  # what expressions are evaluated on, and the expressions for each. javascript
  # reports matches as indexes into units (see sl.samples in runtime.js).
  # sample expressions are group expressions on a "sample" group.
  units: array[UnitKind, seq[Unit]]
  exprs: array[UnitKind, seq[CompiledExpression]]
  has_info_expression*: bool

  # javascript functions and objects
  js_info: JSObjectRef           # sl.info
  js_samples: JSObjectRef        # sl.samples
  js_errs: JSName                # "ERRS"

  # buffers shared with javascript (see runtime.js for what each holds)
  meta: Buf[int32]
  vnum: Buf[float64]
  alts: Buf[int8]
  ab: Buf[float32]
  fmt_n: Buf[int32]
  fmt_str: Buf[int32]
  fmt_bufs: seq[FmtBuf]
  info_n: Buf[int32]
  info_t: Buf[int8]
  info_off: Buf[int32]
  info_num: Buf[float64]
  str: Buf[uint8]
  str_off: Buf[int32]
  str_len: Buf[int32]
  out_buf: Buf[int32]

  # per-variant state
  n_info_num: int                # used length of info_num
  n_str: int                     # used length of str
  n_slots: int                   # number of strings (slots) in str
  last_info_ids: seq[int]        # INFO / FORMAT ids set for the current variant,
  last_fmt_ids: seq[int]         # reset before the next one

proc slivarDebug(ctx: JSContextRef, function: JSObjectRef, thisObject: JSObjectRef, argumentCount: csize_t,
                 arguments: ptr UncheckedArray[JSValueRef], exception: ptr JSValueRef): JSValueRef {.cdecl.} =
  var parts = newSeq[string](argumentCount.int)
  for i in 0..<argumentCount.int:
    parts[i] = arguments[i].toStringSafe(ctx)  # callbacks must not raise
  stderr.write_line parts.join(", ")
  return JSValueMakeUndefined(ctx)

proc trio_kids*(samples: seq[Sample], by_name: TableRef[string, ISample]): seq[Sample] =
  ## return all samples that have a mom and dad that are present in the vcf
  result = newSeqOfCap[Sample](16)
  for sample in samples:
    if sample.mom == nil or sample.mom.i == -1: continue
    if sample.dad == nil or sample.dad.i == -1: continue
    if by_name != nil:
      if sample.id notin by_name: continue
      if sample.mom.id notin by_name: continue
      if sample.dad.id notin by_name: continue
    result.add(sample)

const prelude = staticRead("prelude.js")
const runtime = block:
  const marker = "  /*@LAYOUT@*/\n"
  const src = staticRead("runtime.js")
  static: doAssert marker in src, "runtime.js is missing the layout marker"
  src.replace(marker, jsLayout())

proc getEnvNotEmpty(key: string): seq[string] =
  var s = getEnv(key)
  if s == "": return
  for f in s.split(","):
    result.add(f)

proc id2names*(h:Header): seq[idpair] =
  ## lookup of headerid -> name
  var hdr = h.hdr
  result = newSeq[idpair](hdr.n[0])
  let arr = cast[ptr UncheckedArray[bcf_idpair_t]](hdr.id[0])
  for i in 0..<hdr.n[0].int32:
    var idp = arr[i]
    if idp.val == nil or idp.key == nil: continue
    if idp.key.len == 0: continue
    var name = idp.key
    if idp.val.hrec[1] == nil and idp.val.hrec[2] == nil: continue
    if idp.val.id >= result.len:
      result.setLen(idp.val.id + 2)
    result[idp.val.id] = ($name, idp.val.hrec[1] != nil, idp.val.hrec[2] != nil)

proc js_add_expression(ev: Evaluator, kind, name, expr: string) =
  let ctx = ev.b.ctx
  let f = ctx.getProp(ev.b.sl, "addExpr").toObject(ctx)
  try:
    discard ctx.call(f, [ctx.toJS(kind), ctx.toJS(name), ctx.toJS(expr)], ev.b.sl)
  except JSError as e:
    raise newException(ValueError, e.msg & &"\nexpression was:'{expr}'")

proc newEvaluator*(ivcf:VCF, samples: seq[Sample], groups: seq[Group],
                   trio_expressions: seq[NamedExpression],
                   group_expressions: seq[NamedExpression],
                   family_expressions: seq[NamedExpression],
                   sample_expressions: seq[NamedExpression],
                   info_expr: string, gnos:seq[Gnotater], skip_non_variable_sites: bool): Evaluator =
  ## make a new evaluation context for the given string
  let ctx = newJSContext()
  result = Evaluator(b: Bridge(ctx: ctx))

  result.info_white_list = getEnvNotEmpty("SLIVAR_INFO_WHITELIST")
  result.format_white_list = getEnvNotEmpty("SLIVAR_FORMAT_WHITELIST")
  result.allow_fmt_strings = getEnv("SLIVAR_FORMAT_STRINGS") != ""

  var vcf_obj = newJObject()
  for f in ["ANN", "CSQ", "BCSQ"]:
    try:
      var gf:GeneIndexes
      var fields = ivcf.set_csq_fields(f, gf)
      if fields.len > 0:
        result.gene_fields.add(gf)
        vcf_obj[f] = %fields
    except KeyError:
      continue
  if result.gene_fields.len > 0:
    for k in ["impactful", "genic"]:
      if ivcf.header.add_info(k, "0", "Flag", &"variant is {k} according to slivar and impacts added by VEP/snpEff/bcftools (https://github.com/brentp/slivar/wiki/impactful)") != Status.OK:
        var ok = false
        try:
          var hi = ivcf.header.get(k, BCF_HL_INFO)
          if hi["Type"] == "Flag": ok = true
        except KeyError:
          discard

        if not ok:
          stderr.write_line &"[slivar] couldn't add '{k}' to header as it already existed"
    discard ivcf.header.add_info("highest_impact_order", "1", "Integer", "impact order (lower is higher) of this variant across all genes and transcripts it overlaps. this integer can be used as a look into the order list to get the actual impact")

  # the header may have changed above so the ids are calculated here.
  result.field_names = id2names(ivcf.header)
  for i, f in result.field_names:
    if not f.info: continue
    if f.name == "impactful": result.impactful_id = i
    elif f.name == "genic": result.genic_id = i
    elif f.name == "highest_impact_order": result.highest_impact_id = i

  let strict = getEnv("SLIVAR_NO_ATTRIBUTE_CHECKING") == ""
  if not strict:
    if getEnv("SLIVAR_QUIET") == "":
      stderr.write_line "[slivar] using non-strict mode. missing attributes will not show a warning or error"

  # need this because we can only have 1 object per sample id. this allows fast lookup by id.
  var by_name = newTable[string,ISample]()
  var sample_index = newTable[string, int]()
  for sample in samples:
    result.samples.add(ISample(ped_sample:sample))
    by_name[sample.id] = result.samples[result.samples.high]
    sample_index[sample.id] = result.samples.high

  discard ctx.eval(runtime, "runtime.js")
  ctx.setGlobal("debug", JSValueRef(ctx.newFunction("debug", slivarDebug)))
  discard ctx.eval(prelude, "prelude.js")
  result.b.sl = ctx.getGlobal("$sl").toObject(ctx)
  ctx.protect(result.b.sl)

  result.gnos = gnos

  # --- configuration for the javascript side. samples are referred to by
  # their index in ev.samples.
  var cfg = newJObject()
  cfg["strict"] = %strict
  cfg["skip_non_variable"] = %skip_non_variable_sites
  var info_names = newJArray()
  var fmt_names = newJArray()
  for f in result.field_names:
    info_names.add(%(if f.info: f.name else: ""))
    fmt_names.add(%(if f.format: f.name else: ""))
  cfg["info"] = info_names
  cfg["fmt"] = fmt_names
  proc idx(s: Sample): int = sample_index.getOrDefault(s.id, -1)
  var jsamples = newJArray()
  for sample in result.samples:
    let s = sample.ped_sample
    var kids = newJArray()
    for kid in s.kids:
      if kid.idx >= 0: kids.add(%kid.idx)
    jsamples.add(%*{"id": s.id, "i": s.i, "affected": s.affected, "phenotype": s.phenotype,
               "sex": (if s.sex == 2: "female" elif s.sex == 1: "male" else: "unknown"),
               "dad": (if s.dad != nil: s.dad.idx else: -1),
               "mom": (if s.mom != nil: s.mom.idx else: -1), "kids": kids})
  cfg["samples"] = jsamples

  # units: the globals to set (names, values) before evaluating an expression.
  # a value is a sample index, an array of them (plural group columns, families)
  # or null (an empty group column).
  var junits = newJObject()
  for kind in UnitKind: junits[$kind] = newJArray()
  template addUnit(kind: UnitKind, u: Unit, names: openArray[string], values: JsonNode, require_alt = newJNull()) =
    result.units[kind].add(u)
    junits[$kind].add(%*{"names": names, "values": values, "ids": u.ids, "require_alt": require_alt})

  for kid in samples.trio_kids(by_name):
    # with skip_non_variable_sites, a trio is skipped unless one of its samples (vcf indexes) has an alt.
    let require_alt = if skip_non_variable_sites: %*[kid.i, kid.dad.i, kid.mom.i] else: newJNull()
    addUnit(ukTrio, Unit(ids: @[kid.id]), ["kid", "dad", "mom"], %*[kid.idx, kid.dad.idx, kid.mom.idx], require_alt)
  if trio_expressions.len > 0:
    if result.units[ukTrio].len > 0:
      if getEnv("SLIVAR_QUIET") == "":
        stderr.write_line &"[slivar] evaluating on {result.units[ukTrio].len} trios"
    else:
      stderr.write_line &"[slivar] WARNING! specified --trio expressions without any trios"

  var all_groups = groups
  if sample_expressions.len > 0:
    # each sample is a single-column group named "sample"
    var g = Group(header: @["sample"], plural: @[false])
    for s in samples: g.rows.add(@[@[s]])
    all_groups.add(g)
  for g in all_groups:
    for row in g.rows:
      var values = newJArray()
      for k, col in row:
        if col.len == 0: values.add(newJNull())
        elif g.plural[k]:
          var arr = newJArray()
          for s in col: arr.add(%s.idx)
          values.add(arr)
        else: values.add(%col[0].idx)
      # a match reports the first sample of the row (the first column may be empty).
      let ids = if row[0].len > 0: @[row[0][0].id] else: @[]
      addUnit(ukGroup, Unit(ids: ids), g.header, values)

  if family_expressions.len > 0:
    var by_fam = newTable[string, seq[Sample]]()  # (iteration order determines output order)
    for s in samples:
      by_fam.mgetOrPut(s.family_id, @[]).add(s)
    for family_id, members in by_fam:
      var u = Unit(family_id: family_id)
      var arr = newJArray()
      for s in members:
        arr.add(%s.idx)
        if s.affected: u.ids.add(s.id)
      addUnit(ukFamily, u, ["fam"], %*[arr])
  cfg["units"] = junits
  var io = newJObject()
  for imp, idx in default_order: io[imp] = %idx
  cfg["impact_order"] = io
  cfg["vcf"] = vcf_obj

  let setup = ctx.getProp(result.b.sl, "setup").toObject(ctx)
  withJSString(js, $cfg):
    let jcfg = JSValueMakeFromJSONString(ctx.ctx, js)
    discard ctx.call(setup, [jcfg], result.b.sl)

  # --- shared buffers
  let ev = result
  let nvcf = max(1, ivcf.n_samples)
  let nids = max(1, ev.field_names.len)
  ev.b.bind_buf(ev.meta, "META", 8)
  ev.b.bind_buf(ev.vnum, "VNUM", 16)
  ev.b.bind_buf(ev.alts, "ALTS", nvcf)
  ev.b.bind_buf(ev.ab, "AB", nvcf)
  ev.b.bind_buf(ev.fmt_n, "FMT_N", nids)
  ev.b.bind_buf(ev.fmt_str, "FMT_STR", nids)
  ev.fmt_bufs = newSeq[FmtBuf](nids)
  ev.b.bind_buf(ev.info_n, "INFO_N", nids)
  ev.b.bind_buf(ev.info_t, "INFO_T", nids)
  ev.b.bind_buf(ev.info_off, "INFO_OFF", nids)
  ev.b.bind_buf(ev.info_num, "INFO_NUM", 256)
  ev.b.bind_buf(ev.str, "STR", 4096)
  ev.b.bind_buf(ev.str_off, "STR_OFF", 64)
  ev.b.bind_buf(ev.str_len, "STR_LEN", 64)
  let n_out = 1 + (trio_expressions.len * (ev.units[ukTrio].len + 1)) +
              ((group_expressions.len + sample_expressions.len) * (ev.units[ukGroup].len + 1)) +
              (family_expressions.len * (ev.units[ukFamily].len + 1))
  ev.b.bind_buf(ev.out_buf, "OUT", n_out)

  # --- expressions. sample expressions are group expressions (on the "sample" group).
  template add(kind: string, exprs: seq[NamedExpression], dest: var seq[CompiledExpression]) =
    for ex in exprs:
      ev.js_add_expression(kind, ex.name, ex.expr)
      dest.add(CompiledExpression(expr: ex.expr, name: ex.name))
  add($ukTrio, trio_expressions, ev.exprs[ukTrio])
  add($ukGroup, group_expressions, ev.exprs[ukGroup])
  add($ukGroup, sample_expressions, ev.exprs[ukGroup])
  add($ukFamily, family_expressions, ev.exprs[ukFamily])

  if info_expr != "" and info_expr != "nil":
    ev.js_add_expression("info", "", info_expr)
    ev.has_info_expression = true

  ev.js_info = ctx.getProp(ev.b.sl, "info").toObject(ctx)
  ev.js_samples = ctx.getProp(ev.b.sl, "samples").toObject(ctx)
  for o in [ev.js_info, ev.js_samples]: ctx.protect(o)
  ev.js_errs = jsName("ERRS")

proc load_js*(ev:Evaluator, code:string) =
  try:
    discard ev.b.ctx.eval(code, "--js")
  except JSError as e:
    quit "error evaluating code in:" & code & ":" & e.msg

# ---------------------------------------------------------------- strings

proc add_str(ev: Evaluator, s: openArray[char]): int {.inline.} =
  ## append a string to the arena and return its slot.
  result = ev.n_slots
  ev.b.ensure(ev.str_off, "STR_OFF", ev.n_slots + 1)
  ev.b.ensure(ev.str_len, "STR_LEN", ev.n_slots + 1)
  ev.b.ensure(ev.str, "STR", ev.n_str + s.len)
  if s.len > 0:
    copyMem(ev.str.data[ev.n_str].addr, s[0].unsafeAddr, s.len)
  ev.str_off[result] = ev.n_str.int32
  ev.str_len[result] = s.len.int32
  ev.n_str += s.len
  ev.n_slots.inc

proc add_str(ev: Evaluator, s: cstring): int {.inline.} =
  ev.add_str(toOpenArray(s, 0, s.len - 1))

proc add_str(ev: Evaluator, s: string): int {.inline.} =
  ev.add_str(toOpenArray(s, 0, s.len - 1))

# ---------------------------------------------------------------- per variant

proc set_variant_fields*(ev:Evaluator, variant:Variant) =
  ev.meta[msGeneration] = ev.meta[msGeneration] + 1
  ev.n_str = 0
  ev.n_slots = 0
  # the first string slots, in the order of VariantString
  let alts = variant.ALT
  for vs in VariantString:
    case vs
    of vsCHROM: discard ev.add_str(variant.CHROM)
    of vsREF: discard ev.add_str(variant.REF)
    of vsALT: discard ev.add_str(alts.join(","))
    of vsID: discard ev.add_str(variant.ID)
    of vsFILTER: discard ev.add_str(variant.FILTER)
  ev.vnum[vnStart] = variant.start.float64
  ev.vnum[vnStop] = variant.stop.float64
  ev.vnum[vnPOS] = variant.POS.float64
  ev.vnum[vnQUAL] = variant.QUAL
  ev.vnum[vnMultiallelic] = if alts.len > 1: 1.0 else: 0.0
  ev.meta[msCalculated] = 0
  if alts.len > 1:
    stderr.write_line &"[slivar] warning! found multiallelic site at {variant.CHROM}:{variant.start + 1}. slivar will do sane things with this, but it's better to decompose"

proc sum(counts: array[4, int]): int {.inline.} =
    return counts[0] + counts[1] + counts[2] + counts[3]

template aaf*(counts:array[4, int]): float64 =
  ## alternate allele frequency
  if counts[3] != counts.sum():
    float64(2 * counts[2] + counts[1]) / float64(2 * counts.sum - 2 * counts[3])
  else:
    0

proc hwe_score*(counts: array[4, int], aaf:float64): float64 {.inline.} =
  ## calculate the hardy-weinberg chi-sq deviation from expected. values > 6 are unlikely.
  ## counts is [num_hom_ref, hum_het, hum_hom_alt, num_unknown]
  var
    n_called = float64(counts[0] + counts[1] + counts[2])
    raf = 1 - aaf
    exp_hom_ref = (raf ^ 2) * n_called
    exp_het = (2.0 * (raf * aaf)) * n_called
    exp_hom_alt = (aaf ^ 2) * n_called

  result = ((counts[0].float64 - exp_hom_ref) ^ 2) / max(1, exp_hom_ref)
  result += ((counts[1].float64 - exp_het) ^ 2) / max(1, exp_het)
  result += ((counts[2].float64 - exp_hom_alt) ^ 2) / max(1, exp_hom_alt)

proc set_calculated_variant_fields*(ev:Evaluator, alts: var seq[int8]) =
  # homref, het, homalt, unknown (-1)
  if len(alts) == 0:
    return
  var counts = [0, 0, 0, 0]
  for a in alts:
    if unlikely(a < 0 or a > 2):
      counts[3] += 1
    else:
      counts[a] += 1

  var aaf = counts.aaf
  ev.meta[msCalculated] = 1
  ev.vnum[vnAaf] = aaf
  ev.vnum[vnHwe] = hwe_score(counts, aaf)
  ev.vnum[vnCallRate] = 1 - (counts[3].float64 / alts.len.float64)
  ev.vnum[vnNumHomRef] = counts[0].float64
  ev.vnum[vnNumHet] = counts[1].float64
  ev.vnum[vnNumHomAlt] = counts[2].float64
  ev.vnum[vnNumUnknown] = counts[3].float64

var info_warn = 0

proc write_warning(variant:Variant, msg: string, nerrors: var int) {.inline.} =
  nerrors.inc
  if nerrors < 10:
    stderr.write_line "[slivar] javascript error. this can some times happen when a field is missing."
    stderr.write_line msg
    stderr.write "[slivar] occured with variant:", variant.tostring()
    stderr.write_line "[slivar] continuing execution."
  if nerrors == 10:
    stderr.write_line "[slivar] not reporting further errors."

# set the value of INFO field `id` for the current variant.

proc set_info(ev: Evaluator, id: int, kind: InfoKind, off: int, n = 1) {.inline.} =
  ev.info_t[id] = kind.int8
  ev.info_off[id] = off.int32
  ev.info_n[id] = n.int32
  ev.last_info_ids.add(id)

proc add_info_num(ev: Evaluator, x: float64): int {.inline.} =
  ## append to info_num and return the offset.
  ev.b.ensure(ev.info_num, "INFO_NUM", ev.n_info_num + 1)
  result = ev.n_info_num
  ev.info_num[result] = x
  ev.n_info_num.inc

template set_info_nums[T](ev: Evaluator, id: int, vals: var seq[T], scalar: bool) =
  ## `scalar` comes from the header (Number=1): a Number=A field with a single
  ## value is still an array.
  let off = ev.n_info_num
  if scalar and vals.len > 0:
    discard ev.add_info_num(vals[0].float64)
    ev.set_info(id, ikNumber, off)
  else:
    for v in vals: discard ev.add_info_num(v.float64)
    ev.set_info(id, ikArray, off, vals.len)

proc set_info_flag(ev: Evaluator, id: int, value: bool) {.inline.} =
  ev.set_info(id, ikFlag, ev.add_info_num(if value: 1.0 else: 0.0))

proc set_info_string(ev: Evaluator, id: int, value: string) {.inline.} =
  ev.set_info(id, ikString, ev.add_str(value))

proc set_infos*(ev:var Evaluator, variant:Variant, ints: var seq[int32], floats: var seq[float32]) =
  var istr: string = ""
  var info = variant.info

  for id in ev.last_info_ids: ev.info_t[id] = ikAbsent.int8
  ev.last_info_ids.setLen(0)
  ev.n_info_num = 0

  # we found a csq field.
  var impact_found = false
  # we found an impactful variant, e.g. in CSQ, and dont' want to overwrite it
  # with non-impactful in ANN
  var impactful_found = false
  var genic_found = false
  var highest_impact = int32.high.int - 1
  let default_order = default_order

  for field in info.fields:
    if field.name == "impactful": continue
    if field.name == "genic": continue
    if ev.info_white_list.len > 0 and field.name notin ev.info_white_list:
      continue
    let id = field.i
    if field.vtype == BCF_TYPE.FLOAT:
      if info.get(field.name, floats) != Status.OK:
        quit "couldn't get field:" & field.name
      ev.set_info_nums(id, floats, field.n == 1)
    elif field.vtype == BCF_TYPE.CHAR:
      var ret = info.get(field.name, istr)
      if ret != Status.OK:
        quit "couldn't get field:" & field.name & " status:" & $ret
        # NOTE: all set as a single string for now.
      ev.set_info_string(id, istr)

      if field.name in ["CSQ", "BCSQ", "ANN"] and ev.gene_fields.len > 0:
        for g in ev.gene_fields:
          if g.csq_field != field.name: continue
          var imp = istr.get_highest_impact(g, default_order) #: tuple[impact: string, order: int, impactful: bool, genic: bool]
          impact_found = true
          highest_impact = min(highest_impact, imp.order)
          if imp.impactful:
            impactful_found = true
          if imp.genic:
            genic_found = true

    elif field.vtype in {BCF_TYPE.INT32, BCF_TYPE.INT16, BCF_TYPE.INT8}:
      if info.get(field.name, ints) != Status.OK:
        quit "couldn't get field:" & field.name
      if field.n == 1:
        if len(ints) == 0:
          if info_warn < 5:
            stderr.write_line &"[slivar] warning {field.name} had empty value for {variant.tostring()} setting to -127"
            info_warn.inc
            if info_warn == 5:
              stderr.write_line &"[slivar] not reporting further warnings of this type."
          ints.add(-127)
      ev.set_info_nums(id, ints, field.n == 1)
    elif field.vtype == BCF_TYPE.NULL:
      ev.set_info_flag(id, info.has_flag(field.name))

  if ev.gene_fields.len > 0:
    # newEvaluator added these to the header, so they have ids.
    ev.set_info_flag(ev.impactful_id, impactful_found)
    ev.set_info_flag(ev.genic_id, genic_found)
    ev.set_info(ev.highest_impact_id, ikNumber, ev.add_info_num(highest_impact.float64))
    discard variant.info.set("impactful", impactful_found)
    discard variant.info.set("genic", genic_found)
    discard variant.info.set("highest_impact_order", highest_impact)

proc js_errors(ev: Evaluator): seq[string] =
  let ctx = ev.b.ctx.ctx
  let errs = ctx[ev.b.sl, ev.js_errs].toObject(ctx)
  let n = ctx.getProp(errs, "length").toFloat(ctx).int
  for i in 0..<n:
    result.add ctx[errs, i].toString(ctx)

proc set_ab(ev: Evaluator, ints: var seq[int32], n_alleles:int=2) =
  var
    r: float32
    a: float32

  for s in ev.samples:
    let i = s.ped_sample.i
    r = ints[n_alleles * i].float32
    a = ints[n_alleles * i + 1].float32
    for j in 2..<n_alleles:
      a += ints[n_alleles * i + j].float32
    if r < 0 or a < 0:
      ev.ab[i] = -1.0
    else:
      ev.ab[i] = a / max(a + r, 1)
  ev.meta[msAb] = 1

proc set_format_field(ev: Evaluator, f:FormatField, fmt:FORMAT, gts:var Genotypes, ints: var seq[int32], floats: var seq[float32], strs: var seq[string]): bool =
  ## returns true if the field was set.
  let id = f.i
  template copy_numbers(vals: untyped, buf: untyped) =
    ## copy the values into this field's typed array.
    if f.n_per_sample == 0: return false
    ev.b.ensure(ev.fmt_bufs[id].buf, "FMT_ARR", vals.len, id)
    copyMem(ev.fmt_bufs[id].buf.data, vals[0].addr, vals.len * sizeof(vals[0]))
    ev.fmt_n[id] = f.n_per_sample.int32

  if f.name == "GT":
    # the allele values (0 for REF, -1 for unknown, ...); htslib stores
    # n_per_sample alleles for every sample.
    ints.setLen(gts.len * f.n_per_sample)
    for i in 0..<gts.len:
      for k, allele in gts[i]:
        ints[i * f.n_per_sample + k] = allele.value.int32
    copy_numbers(ints, ints)
    return true

  case f.vtype:
    of BCF_TYPE.INT32, BCF_TYPE.INT16, BCF_TYPE.INT8:
      if fmt.get(f.name, ints) != Status.OK:
        quit "couldn't get format field:" & f.name
      copy_numbers(ints, ints)
    of BCF_TYPE.FLOAT:
      if fmt.get(f.name, floats) != Status.OK:
        quit "couldn't get format field:" & f.name
      copy_numbers(floats, floats)
    of BCF_TYPE.CHAR:
      if not ev.allow_fmt_strings: return false
      if fmt.get(f.name, strs) != Status.OK:
        quit "couldn't get format field:" & f.name
      # one slot per vcf sample, in vcf order.
      ev.fmt_str[id] = ev.n_slots.int32
      for s in strs:
        discard ev.add_str(s)
      ev.fmt_n[id] = FmtString
    else:
      quit "Unknown field type:" & $f.vtype & " in field:" & f.name
  return true

proc set_format_fields*(ev:var Evaluator, v:Variant, gts: var Genotypes, alts: var seq[int8], ints: var seq[int32], floats: var seq[float32], strs: var seq[string]) =
  for id in ev.last_fmt_ids: ev.fmt_n[id] = 0
  ev.last_fmt_ids.setLen(0)
  ev.meta[msAb] = 0

  # without GT, alts (het, ...) are absent and no trio has an alt (skip_non_variable).
  ev.meta[msAlts] = int32(alts.len > 0)
  if alts.len > 0:
    copyMem(ev.alts.data, alts[0].addr, min(alts.len, ev.alts.cap))
  else:
    zeroMem(ev.alts.data, ev.alts.cap)

  var fmt = v.format
  # a real AB field is loaded like any other field (and takes precedence in
  # javascript); slivar's AB computed from AD is the fallback.
  var ab_computed = false
  for f in fmt.fields:
    if ev.format_white_list.len > 0 and f.name notin ev.format_white_list:
      continue
    if not ev.set_format_field(f, fmt, gts, ints, floats, strs): continue
    ev.last_fmt_ids.add(f.i)
    if f.name == "AD":
      let v_alt_len = max(1, v.ALT.len)
      if ev.samples.len > 0 and v.n_samples * (1 + v_alt_len) != ints.len:
        stderr.write_line &"""[slivar] error !!! please decompose and normalize after setting Number=A for the AD field in your VCF header"""
        stderr.write_line &"""         expected 2 values per sample for 'AD' field, but got {ints.len / v.n_samples:.1f} for variant: {v.CHROM}:{v.start + 1}:{v.REF}:{join(v.ALT, ",")}"""
        stderr.write_line """         see: https://github.com/brentp/slivar/wiki/decomposing-and-subsetting-vcfs"""
      if not ab_computed and v.n_samples * (1 + v_alt_len) == ints.len:
        ev.set_ab(ints, v.ALT.len + 1)
        ab_computed = true

template trio_expressions*(ev: Evaluator): seq[CompiledExpression] = ev.exprs[ukTrio]
template group_expressions*(ev: Evaluator): seq[CompiledExpression] = ev.exprs[ukGroup]
template family_expressions*(ev: Evaluator): seq[CompiledExpression] = ev.exprs[ukFamily]

proc has_sample_expressions*(ev: Evaluator): bool {.inline.} =
  for kind in UnitKind:
    if ev.exprs[kind].len > 0: return true

iterator evaluate*(ev:var Evaluator, variant:Variant, nerrors:var int): exResult =
  for gno in ev.gnos.mitems:
    discard gno.annotate(variant)

  var ints = newSeq[int32](2 * variant.n_samples)
  var floats = newSeq[float32](2 * variant.n_samples)
  var strs = newSeq[string]()

  ev.set_variant_fields(variant)
  var gts = variant.format.genotypes(ints)
  var alts = gts.alts
  ev.set_calculated_variant_fields(alts)
  ev.set_infos(variant, ints, floats)
  let ctx = ev.b.ctx.ctx

  if not ev.has_info_expression and not ev.has_sample_expressions:
      # if there's no info expression, they might just want gnomad annotations.
      yield ("", @[], 0.0'f32)
  else:
    ## the most expensive part is pulling out the format fields so we only do
    ## that if the info expression passes.
    var info_res = irPass
    if ev.has_info_expression:
      info_res = InfoResult(ctx.call(ev.js_info).toFloat(ctx).int)
      if info_res == irError:
        stderr.write_line "[slivar] error evaluating info expression (this can happen if a field is missing):"
        for e in ev.js_errors: stderr.write_line e

    if info_res != irFail:
      if not ev.has_sample_expressions:
        # if there are no sample expressions, they may just be using gnotate with an info filter.
        yield ("", @[], 0.0'f32)

      elif ev.units[ukTrio].len > 0 or ev.units[ukGroup].len > 0 or ev.units[ukFamily].len > 0:
        ev.set_format_fields(variant, gts, alts, ints, floats, strs)
        let nerr = ctx.call(ev.js_samples).toFloat(ctx).int
        if nerr > 0:
          let msgs = ev.js_errors
          for k in 0..<nerr:
            variant.write_warning(if k < msgs.len: msgs[k] else: "", nerrors)

        var o = 0
        for kind in UnitKind:
          for namedexpr in ev.exprs[kind]:
            let cnt = ev.out_buf[o].int
            o.inc
            var matching = newSeq[string]()
            for k in 0..<cnt:
              let u = ev.units[kind][ev.out_buf[o + k]]
              if kind == ukFamily and u.ids.len == 0:
                yield (namedexpr.name, @[u.family_id], float32.low)
                if nerrors < 10:
                  stderr.write_line &"[slivar] WARNING !!! using --family-expr without any affected samples. Adding family_id: '{u.family_id}' to the INFO list for '{namedexpr.name}'"
                  nerrors.inc
                  if nerrors == 10:
                    stderr.write_line &"[slivar] not reporting futher warnings"
              else:
                matching.add(u.ids)
            o += cnt
            if matching.len > 0:
              yield (namedexpr.name, matching, -1'f32)
