#!/bin/bash
test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest

. ssshtest

set -o nounset

set -e
nim c --boundChecks:on -x:on src/slivar
set +e
exe=./src/slivar

run check_help_works $exe expr --help
assert_exit_code 0
assert_in_stdout "Options:"
rm -f xx.bcf

run check_denovo $exe expr -v tests/ashk-trio.vcf.gz --trio "denovo:kid.alts == 1 && mom.alts == 0 && dad.alts == 0 && (mom.AD[1] + dad.AD[1]) < 2 && kid.GQ > 10 && mom.GQ > 10 && dad.GQ > 10 && kid.DP > 10 && mom.DP > 10 && dad.DP > 10" --ped tests/ashk-trio.ped -o xx.bcf
assert_exit_code 0
assert_equal 0 $(bcftools view -H -i 'FMT/GT[0] == "RR"' xx.bcf | wc -l)
assert_equal "0" "$(bcftools view -H -i 'FMT/DP[0] <= 10 || FMT/DP[1] <= 10 || FMT/DP[2] <= 10' xx.bcf | wc -l)"
assert_equal 11 $(bcftools view -H xx.bcf | wc -l)


run check_denovo_with_filter $exe expr -v tests/ashk-trio.vcf.gz --trio "denovo:variant.FILTER == 'PASS' && kid.alts == 1 && mom.alts == 0 && dad.alts == 0 && (mom.AD[1] + dad.AD[1]) < 2 && kid.GQ > 10 && mom.GQ > 10 && dad.GQ > 10 && kid.DP > 10 && mom.DP > 10 && dad.DP > 10" --ped tests/ashk-trio.ped -o xx.bcf

assert_exit_code 0
assert_equal 4 $(bcftools view -H xx.bcf | wc -l)
#rm xx.bcf

run check_denovo_with_info_filter $exe expr -v tests/ashk-trio.vcf.gz --info "variant.FILTER == 'PASS'" --trio "denovo:kid.alts == 1 && mom.alts == 0 && dad.alts == 0 && (mom.AD[1] + dad.AD[1]) < 2 && kid.GQ > 10 && mom.GQ > 10 && dad.GQ > 10 && kid.DP > 10 && mom.DP > 10 && dad.DP > 10" --ped tests/ashk-trio.ped -o xx.bcf

assert_exit_code 0
assert_equal 4 $(bcftools view -H xx.bcf | wc -l)


### groups

# this should have same result as `check_denovo_with_info_filter`
run check_simple_groups $exe expr -v tests/ashk-trio.vcf.gz -a tests/trio.group --group-expr 'denovo:variant.FILTER == '"'"'PASS'"'"' && kid.alts == 1 && mom.alts == 0 && dad.alts == 0 && (  mom.AD[1] + dad.AD[1]) < 2 && kid.GQ > 10 && mom.GQ > 10 && dad.GQ > 10 && kid.DP > 10 && mom.DP > 10 && dad.DP > 10' -o xx.bcf
assert_exit_code 0
assert_equal 4 $(bcftools view -H xx.bcf | wc -l)


#### filter

rm -f xx.bcf
run check_slivar_gnotate $exe gnotate --expr "variant.call_rate < 0" -o xx.bcf tests/ashk-trio.vcf.gz
assert_exit_code 0
assert_equal 0 $(bcftools view -H xx.bcf | wc -l)

rm -f xx.bcf
run check_slivar_gnotate_all $exe gnotate --expr "variant.call_rate > -1" -o xx.bcf tests/ashk-trio.vcf.gz
assert_exit_code 0
assert_equal $(bcftools view -H tests/ashk-trio.vcf.gz | wc -l) $(bcftools view -H xx.bcf | wc -l)

rm -f xx.bcf
run check_slivar_gnotate_load $exe gnotate --js tests/test-functions.js --expr "call_rate(variant)" -o xx.bcf tests/ashk-trio.vcf.gz
assert_equal 9834 $(bcftools view -H xx.bcf | wc -l)
rm -f xx.bcf

run check_compound_hets $exe compound-hets -v tests/comphet.vcf --ped tests/ashk-trio.ped -o _xx.vcf
assert_exit_code 0
assert_equal $(grep -c expect=yes tests/comphet.vcf) $(grep -c expect=yes _xx.vcf)
assert_equal $(grep -c expect=no _xx.vcf) 0
rm -f _xx.vcf
