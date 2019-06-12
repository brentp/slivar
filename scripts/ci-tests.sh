#!/bin/sh
set -e
apk add -q bash
mkdir test-tmp && cd test-tmp
git clone --depth 1 https://github.com/samtools/htslib
git clone --depth 1 https://github.com/samtools/bcftools
cd bcftools && make -j4 install && cd ../.. && rm -rf test-tmp
nimble install -y && nimble test
bash tests/functional-tests.sh
nim c -d:release -o:/test/slivar_static -d:nsb_static src/slivar
