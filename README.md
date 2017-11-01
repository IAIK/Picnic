Picnic: Post-Quantum Signatures
===============================

The Picnic signature scheme is a family of digital signature schemes secure
against attacks by quantum computers. This repository contains an optimized
implementation of these schemes. The scheme and parameter sets are specified in
the [Picnic Specification Document](https://github.com/Microsoft/Picnic/blob/master/spec.pdf).

A research paper describing the signature scheme is also available:
**Post-Quantum Zero-Knowledge and Signatures from Symmetric-Key Primitives**
Melissa Chase and David Derler and Steven Goldfeder and Claudio Orlandi and
Sebastian Ramacher and Christian Rechberger and Daniel Slamanig and Greg
Zaverucha.
*In Proceedings of ACM CCS 2017*.
*Cryptology ePrint Archive: Report 2017/279*
<http://eprint.iacr.org/2017/279>

Building
--------

First configure the build cmake and then run make:

```sh
mkdir build; cd build
cmake ..
make
```

License
-------

The code is licensed under the MIT license.
