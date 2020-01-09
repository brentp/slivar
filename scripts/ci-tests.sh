#!/bin/sh
set -e
apk add -q bash curl-dev
mkdir test-tmp && cd test-tmp
git clone https://github.com/samtools/htslib
cd htslib && git checkout 1.10  && cd ..
git clone https://github.com/samtools/bcftools
cd bcftools && git checkout 1.10 && autoheader && autoconf && ./configure --disable-libcurl && make -j4 install && cd ../.. && rm -rf test-tmp
nimble install -y && nimble test
nim c -d:danger -o:/test/slivar_static -d:nsb_static src/slivar
