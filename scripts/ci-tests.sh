#!/bin/sh
set -e
apk add -q bash curl-dev
nimble install -y https://github.com/brentp/duktape-nim@#dev
nimble install -y && nimble test
nim c -d:danger -o:/test/slivar_static -d:nsb_static src/slivar
