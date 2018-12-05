#!/bin/bash

BASE=()
while [ "$1" != "--" ]
do
  BASE+=("$1")
  shift
done
shift

powerset() {
  [ $# -eq 0 ] && { echo; return; }
  ( shift; powerset "$@" ) | while read a
  do
    echo $1 $a
    echo $a
  done
}

powerset "$@" | xargs -L 1 bash .ci-build.sh "${BASE[@]}"
