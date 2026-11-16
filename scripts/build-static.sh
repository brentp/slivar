#!/usr/bin/env bash
# Build a fully static slivar binary (linux, glibc, x86_64 or aarch64).
#
# Requirements: a C/C++ toolchain (gcc/g++), make, curl, bzip2, nim + nimble,
# and static zlib, libbz2, liblzma and libdeflate (on debian/ubuntu:
#   apt-get install build-essential curl bzip2 zlib1g-dev libbz2-dev liblzma-dev libdeflate-dev).
#
# htslib is built from source without libcurl so it can be linked statically.
# JavaScriptCore is fetched automatically by jsc-nim on first compile (~260MB download,
# cached in ~/.cache/jsc-nim; set JSC_PREFIX to use another location).
#
# usage: scripts/build-static.sh [output-path]   (default: ./slivar_static)
set -euo pipefail

out=${1:-slivar_static}
htslib_version=${HTSLIB_VERSION:-1.21}
work=${BUILD_DIR:-$(pwd)/.static-build}
mkdir -p "$work"

if [[ ! -f "$work/htslib-$htslib_version/libhts.a" ]]; then
  curl -fsSL -o "$work/htslib.tar.bz2" \
    "https://github.com/samtools/htslib/releases/download/$htslib_version/htslib-$htslib_version.tar.bz2"
  tar -xjf "$work/htslib.tar.bz2" -C "$work"
  (cd "$work/htslib-$htslib_version" && \
    ./configure --disable-libcurl --disable-gcs --disable-s3 --with-libdeflate >/dev/null && \
    make -j"$(nproc)" lib-static >/dev/null)
fi

nim c -d:danger --hints:off \
  --dynlibOverride:libhts --dynlibOverride:hts \
  --passL:"$work/htslib-$htslib_version/libhts.a" \
  --passL:"-static -ldeflate -lz -lbz2 -llzma -lpthread" \
  -o:"$out" src/slivar

strip "$out"
case "$out" in /*) exe="$out";; *) exe="./$out";; esac
"$exe" expr --help >/dev/null
echo "built $out:"
file "$out" 2>/dev/null || ls -la "$out"
