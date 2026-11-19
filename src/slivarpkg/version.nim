import strutils

const slivarVersion* = "0.5.0"
const slivarGitCommit* = staticExec("git rev-parse --verify HEAD 2>/dev/null").strip()
