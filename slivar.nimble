import ospaths
template thisModuleFile: string = instantiationInfo(fullPaths = true).filename

when fileExists(thisModuleFile.parentDir / "src/slivar.nim"):
  # In the git repository the Nimble sources are in a ``src`` directory.
  import src/slivarpkg/version as _
else:
  # When the package is installed, the ``src`` directory disappears.
  import slivarpkg/version as _

# Package

version       = slivarVersion
author        = "Brent Pedersen"
description   = "expression on variants for great good"
license       = "MIT"


# Dependencies

requires "hts >= 0.2.20", "lapper", "nimgen", "binaryheap", "zip", "https://github.com/brentp/duktape-nim#dev", "docopt", "https://github.com/brentp/nim-minizip >= 0.0.2 ", "argparse >= 0.7.0"
srcDir = "src"
installExt = @["nim"]

bin = @["slivar"]

skipDirs = @["tests"]

import ospaths,strutils

task test, "run the tests":
  exec "bash tests/functional-tests.sh"
  exec "nim c --lineDir:on --debuginfo -r --threads:on src/slivarpkg/duko"
  exec "nim c --lineDir:on --debuginfo -r --threads:on src/slivarpkg/pracode"
  exec "nim c --lineDir:on --debuginfo -r --threads:on src/slivarpkg/groups"
  exec "nim c --lineDir:on --debuginfo -r --threads:on src/slivarpkg/comphet"
  exec "nim c --lineDir:on --debuginfo -r --threads:on src/slivarpkg/fastsets"
  exec "nim c --lineDir:on --debuginfo -r --threads:on src/slivarpkg/pedfile"

task docs, "Builds documentation":
  mkDir("docs"/"slivar")
  for file in listfiles("src/slivar"):
    if file.endswith("value.nim"): continue
    if splitfile(file).ext == ".nim":
      exec "nim doc2 --verbosity:0 --hints:off -o:" & "docs" /../ file.changefileext("html").split("/", 1)[1] & " " & file

