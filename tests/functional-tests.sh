#!/bin/bash
test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest

. ssshtest

set -o nounset

set -e
nim c --boundChecks:on -x:on src/slivar
set +e
exe=./src/slivar

run check_help_works $exe trio --help
assert_exit_code 0
assert_in_stdout "Options:"

run check_denovo $exe trio -v tests/ashk-trio.vcf.gz --trio "denovo:kid.alts == 1 && mom.alts == 0 && dad.alts == 0 && (mom.AD[1] + dad.AD[1]) < 2 && kid.GQ > 10 && mom.GQ > 10 && dad.GQ > 10 && kid.DP > 10 && mom.DP > 10 && dad.DP > 10" --ped tests/ashk-trio.ped -o xx.bcf
assert_exit_code 0
assert_equal 0 $(bcftools view -H -i 'FMT/GT[0] == "RR"' xx.bcf | wc -l)
assert_equal "0" "$(bcftools view -H -i 'FMT/DP[0] <= 10 || FMT/DP[1] <= 10 || FMT/DP[2] <= 10' xx.bcf | wc -l)"
assert_equal 11 $(bcftools view -H xx.bcf | wc -l)


run check_denovo_with_filter $exe trio -v tests/ashk-trio.vcf.gz --trio "denovo:variant.FILTER == 'PASS' && kid.alts == 1 && mom.alts == 0 && dad.alts == 0 && (mom.AD[1] + dad.AD[1]) < 2 && kid.GQ > 10 && mom.GQ > 10 && dad.GQ > 10 && kid.DP > 10 && mom.DP > 10 && dad.DP > 10" --ped tests/ashk-trio.ped -o xx.bcf

assert_exit_code 0
assert_equal 4 $(bcftools view -H xx.bcf | wc -l)
#rm xx.bcf

run check_denovo_with_info_filter $exe trio -g gnomad-2.1.zip -v tests/ashk-trio.vcf.gz --info "variant.FILTER == 'PASS'" --trio "denovo:kid.alts == 1 && mom.alts == 0 && dad.alts == 0 && (mom.AD[1] + dad.AD[1]) < 2 && kid.GQ > 10 && mom.GQ > 10 && dad.GQ > 10 && kid.DP > 10 && mom.DP > 10 && dad.DP > 10" --ped tests/ashk-trio.ped -o xx.bcf

assert_exit_code 0
assert_equal 4 $(bcftools view -H xx.bcf | wc -l)
