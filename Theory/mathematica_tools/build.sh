#!/bin/zsh
# Build an evaluated Mathematica notebook from a marked-up source.
# usage: build.sh source.src out.nb [out.pdf|none]
# Requires python3 and a Wolfram kernel; set WolframKernel if wolframscript cannot find one.
here=${0:A:h}
src=${1:A}; nb=${2:A}; pdf=${3:-none}; [[ $pdf != none ]] && pdf=${pdf:A}
work=$(mktemp -d)
python3 "$here/mbuild.py" "$src" "$work/cells.wl" || exit 1
: ${WolframKernel:=/Applications/Wolfram.app/Contents/MacOS/WolframKernel}
export WolframKernel
wolframscript -file "$here/buildnb.wl" "$work/cells.wl" "$nb" none > "$work/log" 2>&1 &
pid=$!
# the front end sometimes lingers after Exit[]; stop once the build has reported
for i in $(seq 1 1800); do sleep 1; grep -q "cells with messages" "$work/log" 2>/dev/null && break; kill -0 $pid 2>/dev/null || break; done
sleep 3; kill $pid 2>/dev/null; pkill -f "buildnb.wl" 2>/dev/null; true
cat "$work/log"
# export the PDF from a fresh front end: exporting inside the build session can deadlock on dynamic content
if [[ $pdf != none ]]; then
  wolframscript -file "$here/exportpdf.wl" "$nb" "$pdf" > "$work/plog" 2>&1 &
  pid=$!
  for i in $(seq 1 300); do sleep 1; grep -q "exported" "$work/plog" 2>/dev/null && break; kill -0 $pid 2>/dev/null || break; done
  sleep 2; kill $pid 2>/dev/null; true
  cat "$work/plog"
fi
