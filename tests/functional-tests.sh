#!/bin/bash
test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest

. ssshtest

set -o nounset

set -e
nim c -d:debug --boundChecks:on -x:on src/slivar
set +e
exe=./src/slivar

run check_help_works $exe expr --help
assert_exit_code 0
assert_in_stdout "Options:"
rm -f xx.bcf


run check_ragged_groups $exe expr --alias tests/ragged.group -v tests/ashk-trio.vcf.gz \
	--group-expr "missing_dad:dad.DP > 100" \
	--group-expr "ok_mom:mom.DP > 100" \
	-o xx.bcf
assert_exit_code 0
# in groups, dad is empty in 2nd row, so we skip in the first column, but mom can have variants
# in the 2nd case
assert_in_stderr "sample	missing_dad	ok_mom
HG002	4295	4552
HG003	0	0
HG004	0	4552"

run check_error_repeated_name $exe expr -v tests/ashk-trio.vcf.gz --sample-expr "high_depth:sample.DP > 800" --trio "high_depth:kid.DP > 800" --ped tests/ashk-trio.ped -o xx.bcf
assert_exit_code 1
assert_in_stderr "duplicate label"

## sample expressions
run check_sample_expressions $exe expr -v tests/ashk-trio.vcf.gz --pass-only --sample-expr "high_depth:sample.DP > 800" --ped tests/ashk-trio.ped -o xx.bcf
assert_in_stderr "sample	high_depth
HG002	45
HG003	31
HG004	40"
assert_exit_code 0

run check_denovo $exe expr -v tests/ashk-trio.vcf.gz --js js/slivar-functions.js --pass-only --trio "denovo:kid.het && mom.hom_ref && dad.hom_ref && (mom.AD[1] + dad.AD[1]) < 2 && kid.GQ > 10 && mom.GQ > 10 && dad.GQ > 10 && kid.DP > 10 && mom.DP > 10 && dad.DP > 10" --ped tests/ashk-trio.ped -o xx.bcf
assert_exit_code 0
assert_equal 0 $(bcftools view -H -i 'FMT/GT[0] == "RR"' xx.bcf | wc -l)
assert_equal "0" "$(bcftools view -H -i 'FMT/DP[0] <= 10 || FMT/DP[1] <= 10 || FMT/DP[2] <= 10' xx.bcf | wc -l)"
assert_equal 11 $(bcftools view -H xx.bcf | wc -l)


run check_denovo_with_filter $exe expr --js js/my.js -v tests/ashk-trio.vcf.gz --pass-only --trio "denovo:variant.FILTER == 'PASS' && kid.alts == 1 && mom.alts == 0 && dad.alts == 0 && (mom.AD[1] + dad.AD[1]) < 2 && kid.GQ > 10 && mom.GQ > 10 && dad.GQ > 10 && kid.DP > 10 && mom.DP > 10 && dad.DP > 10" --ped tests/ashk-trio.ped -o xx.bcf

assert_exit_code 0
assert_equal 4 $(bcftools view -H xx.bcf | wc -l)
#rm xx.bcf

run check_denovo_with_info_filter $exe expr --pass-only -v tests/ashk-trio.vcf.gz --info "variant.FILTER == 'PASS'" --trio "denovo:kid.alts == 1 && mom.alts == 0 && dad.alts == 0 && (mom.AD[1] + dad.AD[1]) < 2 && kid.GQ > 10 && mom.GQ > 10 && dad.GQ > 10 && kid.DP > 10 && mom.DP > 10 && dad.DP > 10" --ped tests/ashk-trio.ped -o xx.bcf
assert_exit_code 0
assert_equal 4 $(bcftools view -H xx.bcf | wc -l)

run check_denovo_with_info_filter_include_all $exe expr -v tests/ashk-trio.vcf.gz --info "variant.FILTER == 'PASS'" --trio "denovo:kid.alts == 1 && mom.alts == 0 && dad.alts == 0 && (mom.AD[1] + dad.AD[1]) < 2 && kid.GQ > 10 && mom.GQ > 10 && dad.GQ > 10 && kid.DP > 10 && mom.DP > 10 && dad.DP > 10" --ped tests/ashk-trio.ped -o xx.bcf
assert_exit_code 0
assert_equal 9870 $(bcftools view -H xx.bcf | wc -l)


run check_aliased_info $exe expr -v tests/ashk-trio.vcf.gz --info "variant.FILTER == 'PASS' && variant.INFO.DP > 10" --pass-only --trio "denovo:kid.alts == 1 && mom.alts == 0 && dad.alts == 0" --ped tests/ashk-trio.ped -o xx.bcf
assert_in_stderr "HG002	12"
assert_equal 12 $(bcftools view -H xx.bcf | wc -l)
assert_exit_code 0

run check_missing_zip $exe expr -g XXXXXXXXXX.zip -v tests/ashk-trio.vcf.gz
assert_in_stderr "error opening XXXXXXXXXX.zip"
assert_exit_code 1


### groups

# this should have same result as `check_denovo_with_info_filter`
run check_simple_groups $exe expr --pass-only -v tests/ashk-trio.vcf.gz -a tests/trio.group --group-expr 'denovo:variant.FILTER == '"'"'PASS'"'"' && kid.alts == 1 && mom.alts == 0 && dad.alts == 0 && (  mom.AD[1] + dad.AD[1]) < 2 && kid.GQ > 10 && mom.GQ > 10 && dad.GQ > 10 && kid.DP > 10 && mom.DP > 10 && dad.DP > 10' -o xx.bcf
assert_exit_code 0
assert_equal 4 $(bcftools view -H xx.bcf | wc -l)


#### filter

rm -f xx.bcf
run check_slivar_gnotate $exe expr --info "variant.call_rate < 0" -o xx.bcf -v tests/ashk-trio.vcf.gz
assert_exit_code 0
assert_equal 0 $(bcftools view -H xx.bcf | wc -l)

rm -f xx.bcf
run check_slivar_gnotate_all $exe expr --info "variant.call_rate > -1" -o xx.bcf -v tests/ashk-trio.vcf.gz
assert_exit_code 0
assert_equal $(bcftools view -H tests/ashk-trio.vcf.gz | wc -l) $(bcftools view -H xx.bcf | wc -l)

rm -f xx.bcf
run check_slivar_gnotate_load $exe expr --js tests/test-functions.js --info "call_rate(variant)" -o xx.bcf -v tests/ashk-trio.vcf.gz
assert_equal 9834 $(bcftools view -H xx.bcf | wc -l)
rm -f xx.bcf

run check_compound_hets $exe compound-hets -v tests/comphet.vcf --ped tests/ashk-trio.ped -o ixx.vcf
assert_exit_code 0
assert_equal $(grep -c expect=yes tests/comphet.vcf) $(grep -c expect=yes ixx.vcf)
assert_equal $(grep -c expect=no ixx.vcf) 0
rm -f ixx.vcf

run check_compound_hets_fields $exe compound-hets -s BCSQ -v tests/comphet.vcf --ped tests/ashk-trio.ped -o ixx.vcf
assert_exit_code 0
assert_equal 0 $(grep -cv ^# ixx.vcf)
rm -f ixx.vcf

run check_compound_hets_fields_matching $exe compound-hets -s ch_samples -v tests/comphet.vcf --ped tests/ashk-trio.ped -o ixx.vcf
assert_exit_code 0
assert_equal 9 $(grep -cv ^# ixx.vcf)
rm -f ixx.vcf


run check_bug $exe expr --vcf tests/bug.vcf --ped tests/bug.ped -o xx.bcf
assert_exit_code 0
rm -f xx.bcf

# TSV
run check_tsv $exe tsv -i DP -c BCSQ --csq-column dna_change -s denovo --ped tests/ashk-trio.ped  tests/with-bcsq.vcf -o xx.tsv
assert_exit_code 0
assert_equal $(grep -c $'DP\t' xx.tsv) 1
assert_equal 4 $(grep -c ^denovo xx.tsv)


run check_family_expr $exe expr -p tests/ashk-trio.ped -v tests/ashk-trio.vcf.gz --trio "dn:mom.hom_ref && kid.het && dad.hom_ref" --family-expr 'dnf:fam.every(function(s) {return s.het == s.affected && s.hom_ref == !s.affected})' -o /dev/null
assert_exit_code 0
assert_in_stderr "sample	dn	dnf
HG002	26	26"

rm xx.tsv
