// slivar's javascript runtime.
//
// Per-variant data is never copied into javascript objects. nim decodes it into
// typed arrays that view nim memory (the fields of `sl` below) and the objects
// that expressions use (INFO, variant, kid, mom, sample, ...) are Proxies that
// read from those arrays when an attribute is accessed.
//
// Per variant, nim calls sl.info() for the --info expression and, if that
// passes, sl.samples() for all trio, group and family expressions.
//
// Everything internal is under the global `$sl`. The constants in CAPITALS that
// are not defined here (META_*, INFO_*, VARIANT_*, RESULT_*, ...) are generated
// from evaluator.nim and inserted at the marker below.
var $sl = (function() {
  "use strict";
  /*@LAYOUT@*/

  // the buffers shared with nim. nim replaces a typed array (`$sl.X = ...`)
  // whenever its buffer has to grow.
  var sl = {
    META: null,       // Int32Array, indexed by META_*
    VNUM: null,       // Float64Array: numeric variant fields, in the order of VARIANT_NUMBERS
    ALTS: null,       // Int8Array: number of alt alleles, by vcf sample index
    AB: null,         // Float32Array: allele balance computed from AD, by vcf sample index
    FMT_N: null,      // Int32Array: values per sample for each format id (0: absent, FMT_STRING: a string)
    FMT_ARR: [],      // the values (a typed array) of each numeric format field (and GT)
    FMT_STR: null,    // Int32Array: first string slot of each string format field
    INFO_T: null,     // Int8Array: INFO_* kind of each info id
    INFO_N: null,     // Int32Array: length of number arrays
    INFO_OFF: null,   // Int32Array: offset into INFO_NUM, or the string slot
    INFO_NUM: null,   // Float64Array
    STR: null,        // Uint8Array: utf-8 bytes of all strings of the current variant
    STR_OFF: null,    // Int32Array: offset into STR of each string slot
    STR_LEN: null,    // Int32Array: length of each string slot. The first slots are VARIANT_STRINGS
    OUT: null,        // Int32Array: matches, written by sl.samples
    ERRS: [],         // error messages from the last call
  };

  var STRICT = true;     // throw on access to missing attributes (unless SLIVAR_NO_ATTRIBUTE_CHECKING)
  var MISSING = {};      // returned internally for "no value"
  var SAMPLE_INDEX = Symbol("slivar-sample-index");  // index in the vcf
  var SAMPLE_SLOT = Symbol("slivar-sample-slot");    // index in cfg.samples
  function generation() { return sl.META[META_GENERATION]; }

  function unknown(p) {
    // toJSON is looked up by JSON.stringify; allow that for debugging.
    if (STRICT && typeof p === "string" && p !== "toJSON") { throw "unknown attribute:" + p; }
    return undefined;
  }
  // a plain object that throws on access to missing attributes.
  var strictHandler = { get: function(t, p) { return p in t ? t[p] : unknown(p); } };
  function strictObject(o) { return STRICT ? new Proxy(o, strictHandler) : o; }

  // ------------------------------------------------------------- reading the buffers
  var strCacheGeneration = -1, strCache = [];
  function str(slot) {  // the string in a slot (decoded once per variant)
    if (generation() !== strCacheGeneration) { strCacheGeneration = generation(); strCache = []; }
    var s = strCache[slot];
    if (s === undefined) { s = strCache[slot] = decode(sl.STR_OFF[slot], sl.STR_LEN[slot]); }
    return s;
  }
  function decode(off, len) {  // utf-8 bytes from STR
    var s = "";
    for (var i = off; i < off + len; i += 8192) {  // in chunks: apply() has an argument limit
      s += String.fromCharCode.apply(null, sl.STR.subarray(i, Math.min(off + len, i + 8192)));
    }
    if (/[\x80-\xff]/.test(s)) {  // not ascii: decode as utf-8
      try { s = decodeURIComponent(escape(s)); } catch (e) { /* invalid utf-8: keep the bytes as latin-1 */ }
    }
    return s;
  }
  function numArray(arr, off, n) {  // a plain Array copy of arr[off..off+n)
    var r = new Array(n);
    for (var k = 0; k < n; k++) { r[k] = arr[off + k]; }
    return r;
  }

  // INFO: ids are vcf header ids.
  function infoValue(t, id) {
    var off = sl.INFO_OFF[id];
    switch (sl.INFO_T[id]) {
      case INFO_NUMBER: return sl.INFO_NUM[off];
      case INFO_STRING: return str(off);
      case INFO_FLAG: return sl.INFO_NUM[off] !== 0;
      case INFO_ARRAY: return numArray(sl.INFO_NUM, off, sl.INFO_N[id]);
    }
    return MISSING;
  }

  // FORMAT: ids are vcf header ids, followed by ids for slivar's own sample
  // attributes (DERIVED). When the vcf has a real FORMAT field with the same
  // name as a derived one, the real field is used, and the derived attribute
  // (FMT_FALLBACK) only for variants that don't have the field.
  var DERIVED = ["alts", "hom_ref", "het", "hom_alt", "unknown", "AB"];
  var KIND_FIELD = 0, KIND_GT = 1;  // other kinds are 2 + the index in DERIVED
  var KIND_ALTS = 2, KIND_HOM_REF = 3, KIND_HET = 4, KIND_HOM_ALT = 5, KIND_UNKNOWN = 6, KIND_AB = 7;
  var FMT_KIND = [], FMT_FALLBACK = [];
  function fmtValue(f, i) {  // value of format id f for vcf sample i
    var kind = FMT_KIND[f];
    if (kind === KIND_FIELD) {
      var n = sl.FMT_N[f];
      if (n === 0) { return MISSING; }
      if (n === 1) { return sl.FMT_ARR[f][i]; }
      if (n === FMT_STRING) { return str(sl.FMT_STR[f] + i); }
      return numArray(sl.FMT_ARR[f], i * n, n);
    }
    if (kind === KIND_GT) {  // always an array, even for haploid calls
      var p = sl.FMT_N[f];
      return p === 0 ? MISSING : numArray(sl.FMT_ARR[f], i * p, p);
    }
    if (kind === KIND_AB) { return sl.META[META_AB] ? sl.AB[i] : MISSING; }
    if (!sl.META[META_ALTS]) { return MISSING; }
    var alts = sl.ALTS[i];
    switch (kind) {
      case KIND_ALTS: return alts;
      case KIND_HOM_REF: return alts === 0;
      case KIND_HET: return alts === 1;
      case KIND_HOM_ALT: return alts === 2;
      case KIND_UNKNOWN: return alts < 0 || alts > 2;
    }
    return MISSING;
  }
  function sampleValue(t, f) {
    var i = t[SAMPLE_INDEX], v = fmtValue(f, i);
    if (v === MISSING && FMT_FALLBACK[f] !== undefined) { v = fmtValue(FMT_FALLBACK[f], i); }
    return v;
  }

  // variant: ids index VARIANT_KEYS.
  var VARIANT_KEYS = VARIANT_NUMBERS.concat(VARIANT_STRINGS, ["INFO"]);
  var ALT_SLOT = VARIANT_STRINGS.indexOf("ALT");
  function variantValue(t, k) {
    if (k < VARIANT_NUMBERS.length) {
      if (k === MULTIALLELIC) { return sl.VNUM[k] !== 0; }
      if (k >= FIRST_CALCULATED && !sl.META[META_CALCULATED]) { return MISSING; }  // only with genotypes
      return sl.VNUM[k];
    }
    var slot = k - VARIANT_NUMBERS.length;
    if (slot === ALT_SLOT) { return sl.STR_LEN[slot] === 0 ? [] : str(slot).split(","); }
    if (slot < VARIANT_STRINGS.length) { return str(slot); }
    return INFO;
  }

  // ------------------------------------------------------------- INFO, variant and sample objects
  // These are Proxies over a plain target object. A "field" is a name whose
  // value comes from the current variant (an INFO/FORMAT field or a variant
  // attribute); any other name is an ordinary, persistent property of the target.
  //
  // Setting or deleting a field from javascript overrides it until the next
  // variant (like the old behavior, where slivar re-set all fields per variant).
  // Arrays are cached per variant, so changes to e.g. kid.AD are seen by later
  // reads. Both are kept in `cached`, valid while `cachedAt` (same index) is
  // the current generation.
  //
  //   env:   {MISSING, unknown, generation} from this file
  //   ids:   field name -> field id
  //   names: all field names (for Object.keys)
  //   value(target, id):  the current value of a field, or MISSING
  //   slot(target, id):   index into `cached` (distinct for each target and field)
  //
  // Each kind of object gets its own copy of this function (see makeHandler),
  // so it may only use its arguments.
  function fieldHandler(env, ids, names, value, slot, nslots) {
    "use strict";  // copies made by makeHandler don't inherit the outer directive
    var MISSING = env.MISSING, unknown = env.unknown, generation = env.generation;
    var cachedAt = new Int32Array(nslots).fill(-1), cached = new Array(nslots).fill(null);

    function current(t, id) {
      var k = slot(t, id);
      if (cachedAt[k] === generation()) { return cached[k]; }
      var v = value(t, id);
      if (v !== MISSING && typeof v === "object") { cachedAt[k] = generation(); cached[k] = v; }
      return v;
    }
    function override(t, id, v) {
      var k = slot(t, id);
      cachedAt[k] = generation();
      cached[k] = v;
    }
    function field(t, p) {  // the value of field p, or MISSING (also when p is not a field)
      var id = ids[p];
      return id === undefined ? MISSING : current(t, id);
    }
    return {
      get: function(t, p) {
        var v = field(t, p);
        if (v !== MISSING) { return v; }
        return p in t ? t[p] : unknown(p);
      },
      set: function(t, p, v) {
        var id = ids[p];
        if (id === undefined) { t[p] = v; } else { override(t, id, v); }
        return true;
      },
      deleteProperty: function(t, p) {
        var id = ids[p];
        if (id !== undefined) { override(t, id, MISSING); }
        return delete t[p];
      },
      defineProperty: function(t, p, d) {
        var id = ids[p];
        if (id === undefined || !("value" in d)) { return Reflect.defineProperty(t, p, d); }
        override(t, id, d.value);
        return true;
      },
      has: function(t, p) { return field(t, p) !== MISSING || p in t; },
      ownKeys: function(t) {
        var keys = names.filter(function(n) { return field(t, n) !== MISSING; });
        return keys.concat(Object.keys(t).filter(function(k) { return keys.indexOf(k) === -1; }));
      },
      getOwnPropertyDescriptor: function(t, p) {
        var v = field(t, p);
        if (v === MISSING) { return Reflect.getOwnPropertyDescriptor(t, p); }
        return {value: v, writable: true, enumerable: true, configurable: true};
      },
    };
  }

  // a new copy of fieldHandler (from its source), so that the JIT optimizes the
  // handler of each kind of object separately: sharing one copy made the
  // 300-sample cohort benchmark ~8% slower. The source text of each copy must
  // differ, or JSC's code cache shares the compiled code anyway.
  var handlerCopies = 0;
  function makeHandler(ids, names, value, slot, nslots) {
    var copy = (0, eval)("(" + fieldHandler.toString() + ")/*copy " + (handlerCopies++) + "*/");
    return copy({MISSING: MISSING, unknown: unknown, generation: generation}, ids, names, value, slot, nslots);
  }

  var INFO = null, variant = null;

  // ------------------------------------------------------------- setup (called once by nim)
  // units[kind] are the trios, group rows or families (built by nim) that the
  // expressions in exprs[kind] are evaluated on.
  var units = {}, exprs = {info: []};
  UNIT_KINDS.forEach(function(kind) { units[kind] = []; exprs[kind] = []; });

  sl.setup = function(cfg) {
    STRICT = cfg.strict;

    var infoIds = Object.create(null), infoNames = [];
    cfg.info.forEach(function(name, id) { if (name) { infoIds[name] = id; infoNames.push(name); } });

    var fmtIds = Object.create(null), fmtNames = [];
    cfg.fmt.forEach(function(name, f) {
      if (name) { fmtIds[name] = f; fmtNames.push(name); FMT_KIND[f] = name === "GT" ? KIND_GT : KIND_FIELD; }
    });
    DERIVED.forEach(function(name, k) {
      var id = cfg.fmt.length + k;
      FMT_KIND[id] = KIND_ALTS + k;
      if (fmtIds[name] !== undefined) { FMT_FALLBACK[fmtIds[name]] = id; }
      else { fmtIds[name] = id; fmtNames.push(name); }
    });
    var nfmt = cfg.fmt.length + DERIVED.length;
    sl.FMT_ARR = new Array(cfg.fmt.length);

    var variantIds = Object.create(null);
    VARIANT_KEYS.forEach(function(k, i) { variantIds[k] = i; });

    // the INFO target holds expression results set by slivar and values set by users.
    INFO = new Proxy({}, makeHandler(infoIds, infoNames, infoValue,
      function(t, id) { return id; }, cfg.info.length + 1));
    variant = new Proxy({}, makeHandler(variantIds, VARIANT_KEYS, variantValue,
      function(t, k) { return k; }, VARIANT_KEYS.length));
    var sampleHandler = makeHandler(fmtIds, fmtNames, sampleValue,
      function(t, f) { return t[SAMPLE_SLOT] * nfmt + f; }, cfg.samples.length * nfmt);

    // samples: cfg.samples[j].i is the vcf index; dad, mom and kids refer to indexes in cfg.samples.
    var S = [], byId = {};
    cfg.samples.forEach(function(c, j) {
      var t = {affected: c.affected, phenotype: c.phenotype, sex: c.sex, id: c.id};
      Object.defineProperty(t, SAMPLE_INDEX, {value: c.i, configurable: true});
      Object.defineProperty(t, SAMPLE_SLOT, {value: j, configurable: true});
      S.push(byId[c.id] = new Proxy(t, sampleHandler));
    });
    cfg.samples.forEach(function(c, j) {
      if (c.dad >= 0) { S[j].dad = S[c.dad]; }
      if (c.mom >= 0) { S[j].mom = S[c.mom]; }
      S[j].kids = c.kids.map(function(k) { return S[k]; });
    });

    globalThis.INFO = INFO;
    globalThis.variant = variant;
    globalThis.$S = strictObject(byId);
    globalThis.empty = {};
    globalThis.ImpactOrder = strictObject(cfg.impact_order);
    globalThis.VCF = cfg.vcf;

    function toSample(v) {  // a sample index, an array of them, or null for an empty group column
      if (v === null) { return globalThis.empty; }
      return Array.isArray(v) ? v.map(toSample) : S[v];
    }
    UNIT_KINDS.forEach(function(kind) {
      units[kind] = cfg.units[kind].map(function(u) {
        var values = u.values.map(toSample);
        return {set: setter(u.names, values), values: values, ids: u.ids, requireAlt: u.require_alt};
      });
    });
  };

  // a unit sets the globals named in u.names (kid, dad, mom; the group header;
  // or fam) before an expression is evaluated, using a generated setter such as
  // `globalThis["kid"] = v[0]; ...`, which made the 180-sample cohort benchmark
  // ~10% faster than a loop over the names. Arrays (fam, plural group columns)
  // are copied so that an expression that changes one doesn't affect the next.
  var setters = Object.create(null);  // (no inherited keys: an alias may be named e.g. "constructor")
  function setter(names, values) {
    var lines = names.map(function(n, k) {
      return "globalThis[" + JSON.stringify(n) + "] = v[" + k + "]" + (Array.isArray(values[k]) ? ".slice();" : ";");
    });
    var key = lines.join("\n");
    return setters[key] || (setters[key] = new Function("v", key));
  }

  // ------------------------------------------------------------- expressions
  // returns a function, or null when expr is not a single expression. Then it is
  // evaluated as a script and its completion value is used; nim checks the
  // script syntax (javascript has no way to do that without running it).
  function compile(expr) {
    var src = expr.trim().replace(/;+\s*$/, "");
    try { return new Function("return (" + src + "\n)"); } catch (e) { return null; }
  }

  // returns true if src is evaluated as a script (and its syntax must be checked).
  sl.addExpr = function(kind, name, src) {
    if (!(kind in exprs)) { throw "unknown expression kind:" + kind; }
    var f = compile(src), script = f === null;
    if (script) { f = function() { return (0, eval)(src); }; }
    exprs[kind].push({name: name, src: src, f: f});
    return script;
  };

  function errorString(e, ex) {
    var where = e && e.line ? " (line " + e.line + (e.sourceURL ? " of " + e.sourceURL : "") + ")" : "";
    return "javascript error: " + String(e) + where + " for expression:" + ex.src + "\n";
  }

  // evaluates the --info expression. Returns a RESULT_* code.
  sl.info = function() {
    sl.ERRS.length = 0;
    var ex = exprs.info[0];
    try {
      // like slivar with duktape, only `true` passes (not other truthy values).
      return ex.f() === true ? RESULT_PASS : RESULT_FAIL;
    } catch (e) {
      sl.ERRS.push(errorString(e, ex));
      return RESULT_ERROR;  // treated as a pass
    }
  };

  // evaluates all trio, then group, then family expressions. For each expression
  // OUT gets [number of matches, index of matching unit, ...], and INFO[name] is
  // set to the ids of the matching units (through the INFO proxy, so a result
  // shadows an input field of the same name). Returns the number of errors.
  sl.samples = function() {
    var OUT = sl.OUT, ALTS = sl.ALTS, o = 0, nerrors = 0;
    sl.ERRS.length = 0;
    for (var kk = 0; kk < UNIT_KINDS.length; kk++) {
      var us = units[UNIT_KINDS[kk]], es = exprs[UNIT_KINDS[kk]];
      for (var e = 0; e < es.length; e++) {
        var ex = es[e], countAt = o++, ids = [];
        for (var u = 0; u < us.length; u++) {
          var unit = us[u], r = unit.requireAlt;
          // with skip_non_variable, skip trios where no sample has an alt allele.
          if (r !== null && ALTS[r[0]] <= 0 && ALTS[r[1]] <= 0 && ALTS[r[2]] <= 0) { continue; }
          unit.set(unit.values);
          try {
            if (ex.f() === true) {
              OUT[o++] = u;
              for (var j = 0; j < unit.ids.length; j++) { ids.push(unit.ids[j]); }
            }
          } catch (err) {
            if (nerrors++ < 10) { sl.ERRS.push(errorString(err, ex)); }
          }
        }
        OUT[countAt] = o - countAt - 1;
        if (ids.length > 0) { INFO[ex.name] = ids.join(","); }
        else if (ex.name in INFO) { delete INFO[ex.name]; }
      }
    }
    return nerrors;
  };

  globalThis.strictObject = function() { return strictObject({}); };
  return sl;
})();
