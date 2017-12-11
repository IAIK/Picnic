Picnic: Post-Quantum Signatures
===============================

The Picnic signature scheme is a family of digital signature schemes secure
against attacks by quantum computers. This repository contains an optimized
implementation of these schemes. The scheme and parameter sets are specified in
the [Picnic Specification Document](https://github.com/Microsoft/Picnic/blob/master/spec.pdf).
The public API of the library and the serialization format is compatible with
the [reference implementation](https://github.com/Microsoft/Picnic).

A research paper describing the signature scheme is also available:
**Post-Quantum Zero-Knowledge and Signatures from Symmetric-Key Primitives**
Melissa Chase and David Derler and Steven Goldfeder and Claudio Orlandi and
Sebastian Ramacher and Christian Rechberger and Daniel Slamanig and Greg
Zaverucha.
*In Proceedings of ACM CCS 2017*.
*[Cryptology ePrint Archive: Report 2017/279](http://eprint.iacr.org/2017/279)*

Building
--------

First configure the build with `cmake` and then run `make`:
```sh
mkdir build; cd build
cmake ..
make
```

The cmake based build system supports the following flags:
 * ``WITH_SIMD_OPT``:  Enable SIMD optimizations.
 * ``WITH_AVX2``: Use AVX2 if available.
 * ``WITH_SSE2``: Use SSE2 if available.
 * ``WITH_NEON``: Use NEON if available.
 * ``WITH_MARCH_NATIVE``: Build with -march=native -mtune=native (if supported).
 * ``WITH_LTO``: Enable link-time optimization (if supported).
 * ``WITH_MUL_M4RI``: Use methods of four russians for matrix multiplication.
 * ``WITH_REDUCED_LINEAR_LAYER``: Enable partial pre-computation of round key.

On Windows the code can be built using Visual Studio and cmake's Visual Sutdio
solution generator as follows: Open the "Developer Command Prompt for VS 2017"
and from the source folder, run:
```sh
mkdir build; cd build
cmake -G "Visual Studio 15 2017 Win64" ..
msbuild  /t:Rebuild /p:Configuration=Release picnic.sln
```

After running `cmake`, one can also open and build the solution directly with
the Visual Studio IDE. The code was tested using `cmake' for Windows version
3.10 and Visual Studio 15 2017.

License
-------

The code is licensed under the MIT license.
