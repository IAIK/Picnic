#!/bin/bash
set -e

directory=build-$(sha1sum <<< "$@" | awk '{print $1}')
set -x

mkdir -p "$directory"
cd "$directory"
if [[ $# -eq 1 ]] && [[ -z $1 ]]
then
  cmake ..
else
  cmake .. "$@"
fi
cd ..

cmake --build "$directory"
cmake --build "$directory" --target test
