#!/bin/bash
set -e

# OS X has shasum, but no sha1sum.
if [[ "$TRAVIS_OS_NAME" == "osx" ]]
then
  function sha1sum
  {
    shasum "$@"
  }
fi

directory=build-$(sha1sum <<< "$@" | awk '{print $1}')
set -x

function run_cmake
{
  if [[ -n "$CMAKE_GENERATOR" ]]
  then
    cmake -G "$CMAKE_GENERATOR" "$@"
  else
    cmake "$@"
  fi
}

mkdir -p "$directory"
cd "$directory"
if [[ $# -eq 1 ]] && [[ -z $1 ]]
then
  run_cmake ..
else
  run_cmake .. "$@"
fi
cd ..

cmake --build "$directory"
case "$CMAKE_GENERATOR" in
  Visual*)
    cmake --build "$directory" --target RUN_TESTS
    ;;
  *)
    cmake --build "$directory" --target test
    ;;
esac
