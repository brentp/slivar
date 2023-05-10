import hts/vcf
import hts/private/hts_concat
import hts/files
import math
import pedfile
import ./duko
import os
import ./fastsets
import ./impact_order
import ./tsv
import ./gnotate
import ./groups
import tables
import strformat
import strutils

proc clean(s:string): string =
  return s.replace("\"", "").replace("'", "")

type CompiledExpression = object
  name*: string
  expr: Dukexpr

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


type ISample = ref object
  ped_sample*: pedfile.Sample
  duk: Duko
  last_alts: int

type idpair* = tuple[name:string, info:bool, format:bool]

type Trio = array[3, ISample] ## kid, dad, mom

type Family = seq[ISample]

type IGroup* = object
  ## this is a copy of groups.Group, except we use ISample in place of Sample.
  header*: seq[string]
  plural*: seq[bool]
  # even a single colum is a @[Sample] so we need the triply nested level here.
  rows*: seq[seq[seq[ISample]]]

# at each iteration, we need to make sure there is no stale data left in any
# of the `duk` objects. this could happen if, e.g. variant 123 had a `popmax_AF`
# in the info field, but variant 124 did not.
# Due to the implementation, we need tmp to keep the difference without extra allocations.
# these sets keep the bcf_info_t.key and bcf_fmt_t.id
type FieldSets[T: uint8|uint16] = object
  last: set[T]
  curr: set[T]

type Evaluator* = ref object
  ctx: DTContext
  samples*: seq[Isample]

  field_names: seq[idpair]

  gene_fields: seq[GeneIndexes]

  info_field_sets: FieldSets[uint16]
  fmt_field_sets: FieldSets[uint8]
  allow_fmt_strings: bool

  empty: Duko

  # samples_ns is a name-space to store the samples.
  samples_ns: Duko

  VCF: Duko
  ImpactOrder: Duko
  INFO: Duko
  #info_impact_set : Duko # set of impacts seen in this variant
  variant: Duko
  gnos*:seq[Gnotater]

  trios: seq[Trio]
  families: seq[Family]
  trio_expressions*: seq[CompiledExpression]

  groups: seq[IGroup]
  group_expressions*: seq[CompiledExpression]

  sample_expressions: seq[CompiledExpression]
  family_expressions*: seq[CompiledExpression]

  float_expressions: seq[CompiledExpression]
  info_expression*: Dukexpr

  skip_non_variable_sites: bool

  info_white_list: seq[string]
  format_white_list: seq[string]

template fill[T: int8 | int32 | float32 | string](sample:ISample, name:string, values:var seq[T], nper:int) =
  when T is string:
      sample.duk[name] = values[sample.ped_sample.i]
  else:
    if nper == 1:
      sample.duk[name] = values[sample.ped_sample.i]
    elif nper >= 2:
      sample.duk[name] = values[(nper*sample.ped_sample.i)..<(nper*(sample.ped_sample.i+1))]


template set_with_obj*(ctx:DTContext, idx: int, key:string, value: bool) =
    ## set the property at key to a value to a heap-pointer at idx
    ctx.duk_push_boolean(value.duk_bool_t)
    doAssert ctx.duk_put_prop_lstring(idx, key, key.len.duk_size_t)

template set_with_obj*(ctx:DTContext, idx: int, key:string, value: int8) =
    ## set the property at key to a value to a heap-pointer at idx
    ctx.duk_push_int(value.duk_int_t)
    doAssert ctx.duk_put_prop_lstring(idx, key, key.len.duk_size_t)

proc set_hom_het_alts(ctx:Evaluator, alts: var seq[int8]) =
  # this is more complex than seems necessary because it tries to minimize
  # the number of api calls. so it unsets the last value and then sets the
  # current value. it also tries to avoid any calls if the current value
  # is the same as the last.
  for sample in ctx.samples:
    var alt = alts[sample.ped_sample.i]
    # don't need to reset if this alt is same as last one.
    if alt == sample.last_alts: continue
    var sctx = sample.duk.ctx
    var idx = sctx.duk_push_heapptr(sample.duk.vptr)
    case sample.last_alts:
      of 0:
        sctx.set_with_obj(idx, "hom_ref", false)
      of 1:
        sctx.set_with_obj(idx, "het", false)
      of 2:
        sctx.set_with_obj(idx, "hom_alt", false)
      of -1:
        sctx.set_with_obj(idx, "unknown", false)
      else:
        sctx.set_with_obj(idx, "hom_ref", false)
        sctx.set_with_obj(idx, "het", false)
        sctx.set_with_obj(idx, "hom_alt", false)
        sctx.set_with_obj(idx, "unknown", false)

    case alt:
      of 0:
        sctx.set_with_obj(idx, "hom_ref", true)
      of 1:
        sctx.set_with_obj(idx, "het", true)
      of 2:
        sctx.set_with_obj(idx, "hom_alt", true)
      else:
        sctx.set_with_obj(idx, "unknown", true)

    sctx.set_with_obj(idx, "alts", alt)
    sctx.pop()
    sample.last_alts = alt

var debug: DTCFunction = (proc (ctx: DTContext): cint {.stdcall.} =
  var nargs = ctx.duk_get_top()
  if nargs == 1:
    stderr.write_line $ctx.duk_safe_to_string(-1)
    return
  discard ctx.duk_push_string(", ")
  ctx.duk_insert(0)
  ctx.duk_join(nargs)
  stderr.write_line $ctx.duk_require_string(-1)
)

const samples_name = "$S"

proc set_sample_attributes(ev:Evaluator, by_name: TableRef[string, ISample]) =
  for sample in ev.samples:
    sample.duk["affected"] = sample.ped_sample.affected
    sample.duk["phenotype"] = sample.ped_sample.phenotype
    var sex = sample.ped_sample.sex
    sample.duk["sex"] = if sex == 2: "female" elif sex == 1: "male" else: "unknown"
    sample.duk["id"] = sample.ped_sample.id
    if sample.ped_sample.dad != nil:
      if sample.ped_sample.dad.id notin by_name: continue #
      sample.duk.alias(by_name[sample.ped_sample.dad.id].duk, "dad")
    if sample.ped_sample.mom != nil:
      if sample.ped_sample.mom.id notin by_name: continue #
      sample.duk.alias(by_name[sample.ped_sample.mom.id].duk, "mom")
    if ev.ctx.duk_peval_string_noresult(cstring(&"{samples_name}[\"{sample.ped_sample.id}\"].kids = []")) != 0:
        quit "error setting sample kids"
    for kid in sample.ped_sample.kids:
      #stderr.write_line &"{samples_name}[\"{sample.ped_sample.id}\"].kids.push({samples_name}[\"{kid.id}\"])"
      if kid.id notin by_name: continue #
      if ev.ctx.duk_peval_string(cstring(&"{samples_name}[\"{sample.ped_sample.id}\"].kids.push({samples_name}[\"{kid.id}\"])")) != 0:
        var err = $ev.ctx.duk_safe_to_string(-1)
        quit &"error setting sample kid: {kid.id} for parent: {sample.ped_sample.id}" & "\n" & err
      else:
        ev.ctx.pop()

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

proc make_one_row(ev:Evaluator, grps: seq[seq[Sample]], by_name: TableRef[string, ISample]): seq[seq[ISample]] =
  for col in grps:
    var v = newSeq[Isample](col.len)
    for i, sample in col:
      v[i] = by_name[sample.id]
    result.add(v)

proc make_igroups(ev:Evaluator, groups: seq[Group], by_name:TableRef[string, ISample], sample_expressions: seq[NamedExpression]): seq[IGroup] =
  ## just copy the groups, but turn each sample into an ISample
  for g in groups:
    var ig = IGroup(header:g.header, plural:g.plural)

    for row in g.rows:
      ig.rows.add(ev.make_one_row(row, by_name))

    result.add(ig)
  if sample_expressions.len == 0: return

  # turn samples into a single-column group for each sample
  var ig = IGroup(header: @["sample"], plural: @[false], rows: @[])
  for sample in ev.samples:
    ig.rows.add(@[@[by_name[sample.ped_sample.id]]])
  result.add(ig)

const prelude = staticRead("prelude.js")

proc getEnvNotEmpty(key: string): seq[string] =
  var s = getEnv(key)
  if s == "": return
  for f in s.split(","):
    result.add(f)

proc newEvaluator*(ivcf:VCF, samples: seq[Sample], groups: seq[Group], float_expressions: seq[NamedExpression],
                   trio_expressions: seq[NamedExpression],
                   group_expressions: seq[NamedExpression],
                   family_expressions: seq[NamedExpression],
                   sample_expressions: seq[NamedExpression],
                   info_expr: string, gnos:seq[Gnotater], field_names:seq[idpair], skip_non_variable_sites: bool): Evaluator =
  ## make a new evaluation context for the given string
  var my_fatal: duk_fatal_function = (proc (udata: pointer, msg:cstring) {.stdcall.} =
    stderr.write_line "slivar fatal error:"
    quit $msg
  )

  result = Evaluator(ctx:duk_create_heap(nil, nil, nil, nil, my_fatal), skip_non_variable_sites:skip_non_variable_sites)
  result.ctx.duk_require_stack_top(50)
  result.field_names = field_names

  result.info_white_list = getEnvNotEmpty("SLIVAR_INFO_WHITELIST")
  result.format_white_list = getEnvNotEmpty("SLIVAR_FORMAT_WHITELIST")
  result.allow_fmt_strings = getEnv("SLIVAR_FORMAT_STRINGS") != ""
  result.VCF = result.ctx.newObject("VCF")

  for f in ["ANN", "CSQ", "BCSQ"]:
    try:
      var gf:GeneIndexes
      var fields = ivcf.set_csq_fields(f, gf)
      if fields.len > 0:
        result.gene_fields.add(gf)
        result.VCF[f] = fields
      # add this to the field names so we can clear it as needed
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


  # when a column is empty, we use this empty object as a place-holder.
  result.empty = result.ctx.newObject("empty")

  let strict = getEnv("SLIVAR_NO_ATTRIBUTE_CHECKING") == ""
  if not strict:
    if getEnv("SLIVAR_QUIET") == "":
      stderr.write_line "[slivar] using non-strict mode. missing attributes will not show a warning or error"

  # need this because we can only have 1 object per sample id. this allows fast lookup by id.
  var by_name = newTable[string,ISample]()

  if result.ctx.duk_peval_string(strictO):
    var err = $result.ctx.duk_safe_to_string(-1)
    raise newException(ValueError, err)

  if result.ctx.duk_peval_string(prelude):
    var err = $result.ctx.duk_safe_to_string(-1)
    raise newException(ValueError, err)

  result.ImpactOrder = result.ctx.newStrictObject("ImpactOrder")
  for imp, idx in defaultOrder: result.ImpactOrder[imp] = idx

  if strict:
    result.samples_ns = result.ctx.newStrictObject(samples_name)

    for sample in samples:
      result.samples.add(ISample(ped_sample:sample, last_alts: -2, duk:result.samples_ns.newStrictObject(sample.id)))
      by_name[sample.id] = result.samples[result.samples.high]
  else:
    result.samples_ns = result.ctx.newObject(samples_name)

    for sample in samples:
      result.samples.add(ISample(ped_sample:sample, last_alts: -2, duk:result.samples_ns.newObject(sample.id)))
      by_name[sample.id] = result.samples[result.samples.high]

  result.gnos = gnos
  discard result.ctx.duk_push_c_function(debug, -1.cint)
  discard result.ctx.duk_put_global_string("debug")

  for ex in trio_expressions:
    result.trio_expressions.add(CompiledExpression(expr: result.ctx.compile(ex.expr), name: ex.name))
  for ex in group_expressions:
    result.group_expressions.add(CompiledExpression(expr: result.ctx.compile(ex.expr), name: ex.name))

  for ex in family_expressions:
    result.family_expressions.add(CompiledExpression(expr: result.ctx.compile(ex.expr), name: ex.name))
  # sample expressions get added to group_expressions
  for ex in sample_expressions:
    result.group_expressions.add(CompiledExpression(expr: result.ctx.compile(ex.expr), name: ex.name))

  for ex in float_expressions:
    result.float_expressions.add(CompiledExpression(expr: result.ctx.compile(ex.expr), name: ex.name))

  # get family groups.
  if result.family_expressions.len > 0:
    var by_fam = newTable[string, seq[Sample]]()
    for s in samples:
      by_fam.mgetOrPut(s.family_id, @[]).add(s)
    for family_id, samples in by_fam:
      var fam = newSeq[Isample]()
      for s in samples:
        fam.add(by_name[s.id])
      result.families.add(fam)


  if info_expr != "" and info_expr != "nil":
    result.info_expression = result.ctx.compile(info_expr)

  result.groups = result.make_igroups(groups, by_name, sample_expressions)

  for kid in samples.trio_kids(by_name):
    result.trios.add([by_name[kid.id], by_name[kid.dad.id], by_name[kid.mom.id]])

  if trio_expressions.len > 0:
    if result.trios.len > 0:
      if getEnv("SLIVAR_QUIET") == "":
        stderr.write_line &"[slivar] evaluating on {result.trios.len} trios"
    else:
      stderr.write_line &"[slivar] WARNING! specified --trio expressions without any trios"

  if strict:
    result.INFO = result.ctx.newStrictObject("INFO")
    result.variant = result.ctx.newStrictObject("variant")
    #if getEnv("SLIVAR_IMPACT_SET") != "":
    #  result.info_impact_set = result.INFO.newStrictObject("info_impact_set")
  else:
    result.INFO = result.ctx.newObject("INFO")
    result.variant = result.ctx.newObject("variant")
    #if getEnv("SLIVAR_IMPACT_SET") != "":
    #  result.info_impact_set = result.INFO.newStrictObject("info_impact_set")
  result.set_sample_attributes(by_name)

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


proc c_memset*(p: pointer, value: cint, size: csize_t): pointer {.
  importc: "memset", header: "<string.h>", discardable.}


proc set_ab(ctx: Evaluator, fmt:FORMAT, ints: var seq[int32], n_alleles:int=2) =
  var
    r: float32
    a: float32

  for s in ctx.samples:
    r = ints[n_alleles * s.ped_sample.i].float32
    a = ints[n_alleles * s.ped_sample.i + 1].float32
    for j in 2..<n_alleles:
      a += ints[n_alleles * s.ped_sample.i + j].float32
    if r < 0 or a < 0:
      s.duk["AB"] = -1.0
    else:
      s.duk["AB"] = a / max(a + r, 1)

proc load_js*(ev:Evaluator, code:string) =
    discard ev.ctx.duk_push_string(code)
    if ev.ctx.duk_peval() != 0:
      var err = ev.ctx.duk_safe_to_string(-1)
      quit "error evaluating code in:" & code & ":" & $err
    ev.ctx.pop()



proc set_format_field(ctx: Evaluator, f:FormatField, fmt:FORMAT, gts:var Genotypes, ints: var seq[int32], floats: var seq[float32], strs: var seq[string]) =

  if f.name == "GT":
      var ialleles = newSeq[int32]()
      for sample in ctx.samples.mitems:
          let gti = gts[sample.ped_sample.i]
          ialleles.setLen(gti.len)
          for i, allele in gti:
              ialleles[i] = allele.value.int32
          sample.duk["GT"] = ialleles
      return

  case f.vtype:
    of BCF_TYPE.INT32, BCF_TYPE.INT16, BCF_TYPE.INT8:
      if fmt.get(f.name, ints) != Status.OK:
        quit "couldn't get format field:" & f.name
      for sample in ctx.samples.mitems:
        sample.fill(f.name, ints, f.n_per_sample)
    of BCF_TYPE.FLOAT:
      if fmt.get(f.name, floats) != Status.OK:
        quit "couldn't get format field:" & f.name
      for sample in ctx.samples.mitems:
        sample.fill(f.name, floats, f.n_per_sample)
    of BCF_TYPE.CHAR:
      if ctx.allow_fmt_strings:
        if fmt.get(f.name, strs) != Status.OK:
          quit "couldn't get format field:" & f.name
        for sample in ctx.samples.mitems:
          sample.fill(f.name, strs, f.n_per_sample)
    else:
      quit "Unknown field type:" & $f.vtype & " in field:" & f.name

proc set_variant_fields*(ctx:Evaluator, variant:Variant) =
  ctx.variant["FILTER"] = variant.FILTER
  ctx.variant["CHROM"] = variant.CHROM
  ctx.variant["start"] = variant.start
  ctx.variant["stop"] = variant.stop
  ctx.variant["POS"] = variant.POS
  ctx.variant["QUAL"] = variant.QUAL
  ctx.variant["REF"] = variant.REF
  ctx.variant["ALT"] = variant.ALT
  ctx.variant["is_multiallelic"] = (len(variant.ALT) > 1)
  ctx.variant["ID"] = variant.ID
  if variant.ALT.len > 1:
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

proc set_calculated_variant_fields*(ctx:Evaluator, alts: var seq[int8]) =
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
  ctx.variant["aaf"] = aaf
  ctx.variant["hwe_score"] = hwe_score(counts, aaf)
  ctx.variant["call_rate"] = 1 - (counts[3].float64 / alts.len.float64)
  ctx.variant["num_hom_ref"] = counts[0]
  ctx.variant["num_het"] = counts[1]
  ctx.variant["num_hom_alt"] = counts[2]
  ctx.variant["num_unknown"] = counts[3]

proc clear_unused_infos(ev: Evaluator, f:FieldSets) {.inline.} =

  for idx in f.last.sub(f.curr):
    ev.INFO.del(ev.field_names[idx].name)

var info_warn = 0

proc write_warning(variant:Variant, nerrors: var int) {.inline.} =
  nerrors.inc
  if nerrors < 10:
    stderr.write_line "[slivar] javascript error. this can some times happen when a field is missing."
    stderr.write_line  getCurrentExceptionMsg()
    stderr.write "[slivar] occured with variant:", variant.tostring()
    stderr.write_line "[slivar] continuing execution."
  if nerrors == 10:
    stderr.write_line "[slivar] not reporting further errors."


proc set_infos*(ev:var Evaluator, variant:Variant, ints: var seq[int32], floats: var seq[float32]) =
  var istr: string = ""
  var info = variant.info

  swap(ev.info_field_sets.last, ev.info_field_sets.curr)
  #zeroMem(ev.info_field_sets.curr.addr, sizeof(ev.info_field_sets.curr))
  ev.info_field_sets.curr = {}

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
    if field.vtype == BCF_TYPE.FLOAT:
      if info.get(field.name, floats) != Status.OK:
        quit "couldn't get field:" & field.name
      if field.n == 1:
          ev.INFO[field.name] = floats[0]
      else:
          ev.INFO[field.name] = floats
    elif field.vtype == BCF_TYPE.CHAR:
      var ret = info.get(field.name, istr)
      if ret != Status.OK:
        quit "couldn't get field:" & field.name & " status:" & $ret
        # NOTE: all set as a single string for now.
      ev.INFO[field.name] = istr

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
            info_warn.inc
          ints.add(-127)
        ev.INFO[field.name] = ints[0]
      else:
        ev.INFO[field.name] = ints
    elif field.vtype == BCF_TYPE.NULL:
      ev.INFO[field.name] = info.has_flag(field.name)
    ev.info_field_sets.curr.incl(field.i.uint16)

  # clear any field in last variant but not in this one.
  ev.clear_unused_infos(ev.info_field_sets)
  ev.variant.alias(ev.INFO, "INFO")

  if ev.gene_fields.len > 0:
    ev.INFO["impactful"] = impactful_found
    ev.INFO["highest_impact_order"] = highest_impact
    ev.INFO["genic"] = genic_found
    discard variant.info.set("impactful", impactful_found)
    discard variant.info.set("genic", genic_found)
    discard variant.info.set("highest_impact_order", highest_impact)

type exResult* = tuple[name:string, sampleList:seq[string], val:float32]

iterator evaluate_floats(ev:Evaluator, nerrors: var int, variant:Variant): exResult =
  for i, namedexpr in ev.float_expressions:
    try:
      var val = namedexpr.expr.asfloat()
      ev.INFO[namedexpr.name] = val
      yield (namedexpr.name, @[], val)
    except:
      stderr.write_line "[slivar] problematic variant:" & variant.tostring()
      raise

iterator evaluate_trios(ctx:Evaluator, nerrors: var int, variant:Variant): exResult =
    var alts: seq[int8]
    if ctx.skip_non_variable_sites:
      var ints:seq[int32]
      alts = variant.format.genotypes(ints).alts
    for i, namedexpr in ctx.trio_expressions:
      var matching_samples = newSeq[string]()
      for trio in ctx.trios:

        if ctx.skip_non_variable_sites and alts[trio[0].ped_sample.i] <= 0 and alts[trio[1].ped_sample.i] <= 0 and alts[trio[2].ped_sample.i] <= 0: continue

        trio[0].duk.alias("kid")
        trio[1].duk.alias("dad")
        trio[2].duk.alias("mom")
        try:
            if namedexpr.expr.check():
              matching_samples.add(trio[0].ped_sample.id)
        except:
          variant.write_warning(nerrors)
      if len(matching_samples) > 0:
        # set INFO of this result so subsequent expressions can use it.
        ctx.INFO[namedexpr.name] = join(matching_samples, ",")
        yield (namedexpr.name, matching_samples, -1'f32)
      elif ctx.INFO.hasKey(namedexpr.name):
        ctx.INFO.del(namedexpr.name)

proc alias_objects*(ctx:DTContext, os: seq[ISample], copyname:string) {.inline.} =
  ## add an array of objects and alias them to a new name.
  var idx = ctx.duk_push_array()
  for i, o in os:
    doAssert ctx.duk_push_heapptr(o.duk.vptr) >= 0
    discard ctx.duk_put_prop_index(idx, i.duk_uarridx_t)

  doAssert ctx.duk_put_global_lstring(copyname, copyname.len.duk_size_t)

iterator evaluate_families(ev:Evaluator, nerrors: var int, variant:Variant): exResult =
  for i, namedexpr in ev.family_expressions:
    var matching = newSeq[string]()
    var fam_matching = newSeq[string]()
    var val = -1'f32
    for fam in ev.families:
      fam_matching.setLen(0)
      # TODO: write alias family.
      ev.ctx.alias_objects(fam, "fam")
      try:
        if namedexpr.expr.check:
          for s in fam:
            if s.ped_sample.affected:
              fam_matching.add(s.ped_sample.id)
          if fam_matching.len == 0:
            yield ($namedexpr.name, @[fam[0].ped_sample.family_id], float32.low)
            if nerrors < 10:
              stderr.write_line &"[slivar] WARNING !!! using --family-expr without any affected samples. Adding family_id: '{fam[0].ped_sample.family_id}' to the INFO list for '{namedexpr.name}'"
              nerrors.inc
              if nerrors == 10:
                stderr.write_line &"[slivar] not reporting futher warnings"
          else:
            matching.add(fam_matching)

      except ValueError:
        variant.write_warning(nerrors)

    if len(matching) > 0:
      ev.INFO[namedexpr.name] = join(matching, ",")
      yield ($namedexpr.name, matching, val)
    elif ev.INFO.hasKey(namedexpr.name):
      ev.INFO.del(namedexpr.name)

iterator evaluate_groups(ev:Evaluator, nerrors: var int, variant:Variant): exResult =
    ## note that every group expression is currently applied to every group.
    ## we may want certain expressions applied to certain groups, but how to let the user
    ## specify the connection?
    for i, namedexpr in ev.group_expressions:
      var matching_groups = newSeq[string]()
      for group in ev.groups:
        for row in group.rows:
          for k, col in row:
            if len(col) == 0:
              # if the column is empty, we alias it to our empty object
              # this prevents us from accessing the object of the previous
              # group.
              ev.empty.alias(group.header[k])
            elif not group.plural[k]:
              # TODO: allow skipping if all are hom-ref
              col[0].duk.alias(group.header[k])
            else:
              ev.ctx.alias_objects(col, group.header[k])
          try:
            if namedexpr.expr.check():
              matching_groups.add(row[0][0].ped_sample.id)
          except:
            variant.write_warning(nerrors)
      if len(matching_groups) > 0:
        # set INFO of this result so subsequent expressions can use it.
        ev.INFO[namedexpr.name] = join(matching_groups, ",")
        yield ($namedexpr.name, matching_groups, -1'f32)
      elif ev.INFO.hasKey(namedexpr.name):
        ev.INFO.del(namedexpr.name)

template clear_unused_formats(ev:Evaluator) =
  #for idx in ev.fmt_field_sets.last - ev.fmt_field_sets.curr:
  for idx in ev.fmt_field_sets.last.sub(ev.fmt_field_sets.curr):
    #if idx in ev.fmt_field_sets.curr: continue
    for sample in ev.samples:
      sample.duk.del(ev.field_names[idx].name)

proc set_format_fields*(ev:var Evaluator, v:Variant, gts: var Genotypes, alts: var seq[int8], ints: var seq[int32], floats: var seq[float32], strs: var seq[string]) =
  # fill the format fields
  swap(ev.fmt_field_sets.last, ev.fmt_field_sets.curr)
  ev.fmt_field_sets.curr = {}

  if alts.len > 0:
    ev.set_hom_het_alts(alts)

  var fmt = v.format
  var has_ad = false
  var has_ab = false
  for f in fmt.fields:
    if ev.format_white_list.len > 0 and f.name notin ev.format_white_list:
      continue
    ev.fmt_field_sets.curr.incl(f.i.uint8)
    if has_ab and f.name == "AB": continue
    ev.set_format_field(f, fmt, gts, ints, floats, strs)
    if f.name == "AD":
      has_ad = true
      let v_alt_len = max(1, v.ALT.len)
      if ev.samples.len * (1 + v_alt_len) != ints.len:
        stderr.write_line &"""[slivar] error !!! please decompose and normalize after setting Number=A for the AD field in your VCF header"""
        stderr.write_line &"""         expected 2 values per sample for 'AD' field, but got {ints.len / ev.samples.len:.1f} for variant: {v.CHROM}:{v.start + 1}:{v.REF}:{join(v.ALT, ",")}"""
        stderr.write_line """         see: https://github.com/brentp/slivar/wiki/decomposing-and-subsetting-vcfs"""
      if not has_ab:
        if ev.samples.len * (1 + valt_len) == ints.len:
          ev.set_ab(fmt, ints, v.ALT.len + 1)
        else:
          ev.fmt_field_sets.curr.excl(f.i.uint8)
        has_ab = true

    elif f.name == "AB":
      has_ab = true

  ev.clear_unused_formats()

proc has_sample_expressions*(ev: Evaluator): bool {.inline.} =
  return ev.trio_expressions.len > 0 or ev.group_expressions.len > 0 or ev.family_expressions.len > 0 or ev.sample_expressions.len > 0

proc try_info_check(ev: var Evaluator): bool {.inline.} =
  try:
    return ev.info_expression.check
  except:
    stderr.write_line "[slivar] error evaluating info expression (this can happen if a field is missing):"
    stderr.write_line getCurrentExceptionMsg()

    # NOTE if e.g. a field is missing, we still call it pass
    return true

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

  if ev.info_expression.ctx == nil and not ev.has_sample_expressions:
      # if there's no info expression, they might just want gnomad annotations.
      yield ("", @[], 0.0'f32)
  else:

    ## the most expensive part is pulling out the format fields so we pull all fields
    ## and set values for all samples.
    ## once all that is done, we evaluate the expressions.
    ## the field_sets make it so that we only clear fields from the duk objects that were
    ## set last variant, but not this variant.

    if ev.info_expression.ctx == nil or ev.try_info_check:
      for r in ev.evaluate_floats(nerrors, variant): yield r

      if not ev.has_sample_expressions:
        # if there are no sample expressions, they may just be using gnotate with an info filter.
        yield ("", @[], 0.0'f32)

      elif ev.trios.len > 0 or ev.groups.len > 0 or ev.families.len > 0:
        ev.set_format_fields(variant, gts, alts, ints, floats, strs)
        for r in ev.evaluate_trios(nerrors, variant): yield r
        for r in ev.evaluate_groups(nerrors, variant): yield r
        for r in ev.evaluate_families(nerrors, variant): yield r

