
#[
# groups are user-specified groups. a group file might look like:
#
# #kid	mom	dad
# smp1	smp2	smp3
# smp4	smp5	smp6
# smp7  smp8    smp9
#
# This would set each row of `smp`s to have the column header. For this example, groups
# are not needed because slivar provides a --trio mode. But we can support more general cases like:
#
# #kid	mom	dad sib
# smp1	smp2	smp3 smp1_sib
# smp4	smp5	smp6 smp4_sib
# smp7  smp8    smp9 smp7_sib
#
# such that now, the 4th column sill be available as "sib" in the javascript expression.
#
# column headers ending with 's' are assumed to be list of samples, so, e.g:
#
# #normal tumors
# N1      T1a,T1b,T1c
# N2      T2a,T2b,T2c
# N3      T3a,T3b,T3c
# N4      T4a   # note only a single here, but will still be a list.
#
# then tumors will be a list that can be accessed in javascript as e.g. tumors[0], tumors[1],
# tumors[2]. Note that even for rows with a single value in the column, it will be presented as a
# list as this is a critical way to facilitate the handling of "ragged" arrays.
#
# multiple groups can be specified per file. Any line starting with '#' is assumed to be a header
# and any line without '#' is assumed to fall under the preceding header line. empty lines are
# ignored.
#
]#

import pedfile
import strutils
import strformat
import tables

type Group* = object
  header*: seq[string]
  plural*: seq[bool]
  # even a single colum is a @[Sample] so we need the triply nested level here.
  rows*: seq[seq[seq[Sample]]]

proc parse_group_header_line(groups:var seq[Group], line:string) =
  var group = Group(header: line.split(seps={'\t'}))
  for h in group.header.mitems:
      h = h.strip()
      group.plural.add(h[h.high].toLowerAscii == 's')
  groups.add(group)

proc parse_group_line(groups:var seq[Group], line:string, sample_lookup:TableRef[string,Sample]) =
  if groups.len == 0:
    raise newException(ValueError, "unexpected line:" & line & " occurred before any header")
  var toks = line.split(seps={'\t'})
  if toks.len != groups[groups.high].header.len:
    var sheader = join(groups[groups.high].header, "\t")
    raise newException(ValueError, "unexpected number of samples for line:'" & line & "' with header:'" &
       sheader & "'\n" & "add empty columns where needed")
  # even a single column goes in as a seq for uniformity.
  var g = newSeqOfCap[seq[Sample]](len(toks))
  for i, snames in toks:
    var col_samples = newSeq[Sample]()
    if snames == "":
      g.add(col_samples)
      continue
    for sname in snames.split(seps={','}):
      if sname notin sample_lookup:
        raise newException(ValueError, &"unknown sample {sname} in line '{line}' for groups")
      col_samples.add(sample_lookup[sname])

    if col_samples.len > 1 and not groups[groups.high].plural[i]:
      let h = groups[groups.high].header[i]
      var msg = &"slivar/groups:got > 1 sample in line {line}, column {i}. {$col_samples}"
      msg &= &"either use a single sample or change header '{h}' to '{h}s' to indicate many samples in a single column"
      raise newException(ValueError, msg)
    g.add(col_samples)
  groups[groups.high].rows.add(g)

proc to_lookup(samples:seq[Sample]): TableRef[string, Sample] =
  ## allow to lookup by sample name.
  result = newTable[string, Sample]()
  for s in samples:
    result[s.id] = s

proc parse_groups*(path: string, samples:seq[Sample]): seq[Group] =
  result = newSeq[Group]()
  var sample_lookup = samples.to_lookup
  for l in path.lines:
    var line = l.strip(chars={'\n', '\r'})
    if line.len == 0: continue
    if line[0] == '#':
      result.parse_group_header_line(line.strip(chars={'#'}, trailing=false))
      continue
    result.parse_group_line(line, sample_lookup)

when isMainModule:

  const tmpFile = "__tmp.groups"

  proc toFile(s:string) =
    var f:File
    if not open(f, tmpFile, fmWrite): quit "couldn't open test file"
    f.write(s)
    f.close()

  import unittest

  suite "groups suite":
    test "that file without header gives error":
      "asdf\t123".toFile
      expect ValueError:
        var samples = @[Sample(id:"asdf")]
        discard parse_groups(tmpFile, samples)

    test "that file with non-matching column numbers gives error":
      "#mom\tdad\tkid\na\tb".toFile
      expect ValueError:
        var samples = @[Sample(id:"asdf")]
        discard parse_groups(tmpFile, samples)

    test "that file with matching columns but unknown sample errors":
      "#mom\tdad\tkid\na\tb\tc\n".toFile
      expect ValueError:
        var samples = @[Sample(id:"a"), Sample(id:"b")]
        var groups = parse_groups(tmpFile, samples)
        echo groups

    test "that file with matching columns and samples passes":
      "#mom\tdad\tkid\na\tb\tc\n".toFile
      var samples = @[Sample(id:"a"), Sample(id:"b"), Sample(id:"c")]
      var groups = parse_groups(tmpFile, samples)
      check groups[0].header == @["mom", "dad", "kid"]
      check groups[0].plural == @[false, false, false]
      echo groups[0]
      check groups[0].rows.len == 1
      check groups[0].rows[0].len == 3
      check groups[0].rows[0][0].len == 1


    test "that multiple groups errors on bad number of samples":
      "#mom\tdad\tkid\na\tb\tc\n#tumor\tnormal\nt1\tt2\tn1".toFile
      var samples = @[Sample(id:"a"), Sample(id:"b"), Sample(id:"c"), Sample(id:"t1"), Sample(id:"t2"), Sample(id:"n1")]
      expect ValueError:
        var groups = parse_groups(tmpFile, samples)
        discard groups

    test "that multiple groups works correctly":
      "#mom\tdad\tkid\na\tb\tc\n#tumor\tnormal\nt1a\tt1b\nt2a\tt2b".toFile
      var samples = @[Sample(id:"a"), Sample(id:"b"), Sample(id:"c"), Sample(id:"t1a"), Sample(id:"t1b"), Sample(id:"t2a"), Sample(id:"t2b")]
      var groups = parse_groups(tmpFile, samples)
      check groups.len == 2
      check groups[0].header == @["mom", "dad", "kid"]
      check groups[0].rows.len == 1

      check groups[1].plural == @[false, false]
      check groups[1].header == @["tumor", "normal"]
      check groups[1].rows.len == 2

    test "that plurality works correctly":
      "#moms\tdad\tkids\na,b\tb\tc\n".toFile
      var samples = @[Sample(id:"a"), Sample(id:"b"), Sample(id:"c")]
      var groups = parse_groups(tmpFile, samples)
      check groups[0].plural == @[true, false, true]
      echo groups

    test "that plural without proper header gives error":
      "#mom\tdad\tkid\na,b\tb\tc\n".toFile
      var samples = @[Sample(id:"a"), Sample(id:"b"), Sample(id:"c")]
      expect ValueError:
        discard parse_groups(tmpFile, samples)

    test "that empty groups are ok":
      var samples = @[Sample(id:"a"), Sample(id:"b"), Sample(id:"c")]
      "#aff\tuns\n\ta,b,c\n".toFile
      var groups = parse_groups(tmpFile, samples)
      check groups.len == 1
      var p = groups[0]
      check p.rows.len == 1
      check p.rows[0][0].len == 0
      check p.rows[0][1].len == 3

    test "that truncated rows are ok":
      let txt = """#affected	carrier	healthy
s1	s2	s3
s4	s5
s6	s7
s8
s9
"""
      txt.toFile
      var samples = @[Sample(id:"s1"),
                    Sample(id:"s2"),
                    Sample(id:"s3"),
                    Sample(id:"s4"),
                    Sample(id:"s5"),
                    Sample(id:"s6"),
                    Sample(id:"s7"),
                    Sample(id:"s8"),
                    Sample(id:"s9"),
                    ]
      expect ValueError:
        var groups = parse_groups(tmpFile, samples)
