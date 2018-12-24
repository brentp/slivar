
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

import bpbiopkg/pedfile
import strutils
import strformat
import tables

type Group* = object
  header: seq[string]
  plural: seq[bool]
  groups: seq[seq[Sample]]

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
    raise newException(ValueError, "unexpected number of samples for line:" & line & " with header:" & sheader)
  # even a single column goes in as a seq for uniformity.
  for i, snames in toks:
    var col_samples = newSeq[Sample]()
    for sname in snames.split(seps={','}):
      if sname notin sample_lookup:
        raise newException(ValueError, &"unknown sample {sname} in line '{line}' for groups")
      col_samples.add(sample_lookup[sname])

    if col_samples.len > 1 and not groups[groups.high].plural[i]:
      quit &"slivar/groups:got > 1 sample in line {line}, column {i}. {col_samples}"
    groups[groups.high].groups.add(col_samples)

proc to_lookup(samples:seq[Sample]): TableRef[string, Sample] =
  ## allow to lookup by sample name.
  result = newTable[string, Sample]()
  for s in samples:
    result[s.id] = s

proc parse_groups*(path: string, samples:seq[Sample]): seq[Group] =
  result = newSeq[Group]()
  var sample_lookup = samples.to_lookup
  for l in path.lines:
    var line = l.strip()
    if line.len == 0: continue
    if line[0] == '#':
      result.parse_group_header_line(line)
      continue
    result.parse_group_line(line, sample_lookup)
