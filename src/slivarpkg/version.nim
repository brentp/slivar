import strutils

const slivarVersion* = "0.3.4"
const slivarGitCommit* = staticExec("git rev-parse --verify HEAD 2>/dev/null").strip()
