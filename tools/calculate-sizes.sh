#!/bin/bash
set -e

insts=(lowmc_128_128_182.c.o lowmc_128_128_20.c.o lowmc_192_192_284.c.o lowmc_192_192_30.c.o
lowmc_256_256_363.c.o lowmc_256_256_38.c.o)

mkdir -p size/plain
pushd size/plain
cmake ../../.. -DWITH_LOWMC_OPT=OFF -DWITH_LTO=OFF
popd
make -C size/plain -j picnic

for f in ${insts[@]}
do
  nm --print-size --size-sort --radix=d size/plain/CMakeFiles/picnic.dir/$f > plain.$f.size
done

mkdir -p size/rrkc
pushd size/rrkc
cmake ../../.. -DWITH_LOWMC_OPT=ORKC -DWITH_LTO=OFF
popd
make -C size/rrkc -j picnic

for f in ${insts[@]}
do
  nm --print-size --size-sort --radix=d size/rrkc/CMakeFiles/picnic.dir/$f > rrkc.$f.size
done

mkdir -p size/ollc
pushd size/ollc
cmake ../../.. -DWITH_LOWMC_OPT=OLLE -DWITH_LTO=OFF
popd
make -C size/ollc -j picnic

for f in ${insts[@]}
do
  nm --print-size --size-sort --radix=d size/ollc/CMakeFiles/picnic.dir/$f > ollc.$f.size
done

for f in ${insts[@]}
do
  echo $f
  echo -n "plain: "
  egrep -v "(lowmc|rounds)" plain.$f.size | awk '{s+=$2} END {print s}'
  echo -n "orkc: "
  egrep -v "(lowmc|rounds)" rrkc.$f.size | awk '{s+=$2} END {print s}'
  echo -n "olle: "
  egrep -v "(lowmc|rounds)" ollc.$f.size | awk '{s+=$2} END {print s}'
  echo -n "olle constants: "
  egrep "precomputed_constant" ollc.$f.size | awk '{s+=$2} END {print s}'
  echo -n "olle keys: "
  egrep "precomputed_round_key" ollc.$f.size | awk '{s+=$2} END {print s}'
  echo -n "olle linear layer: "
  egrep "(Zi_|Z_|Ri_)" ollc.$f.size | awk '{s+=$2} END {print s}'
done
