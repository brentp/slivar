#!/bin/sh
set -e
apk add -q bash
mkdir test-tmp && cd test-tmp
git clone https://github.com/samtools/htslib
cd htslib && git checkout 1.9  && cd ..
git clone https://github.com/samtools/bcftools
cd bcftools && git checkout 1.9 && make -j4 install && cd ../.. && rm -rf test-tmp
nimble install -y && nimble test
nim c -d:danger -o:/test/slivar_static -d:nsb_static src/slivar
