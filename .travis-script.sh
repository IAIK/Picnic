#!/bin/bash

powerset() {
  [ $# -eq 0 ] && { echo; return; }
  ( shift; powerset "$@" ) | while read a
  do
    echo $1 $a
    echo $a
  done
}

powerset "$@" | xargs -L 1 -d '\n' bash .travis-build.sh
