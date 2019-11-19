/*
 *  This file is part of the optimized implementation of the Picnic signature scheme.
 *  See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */

#if defined(LOWMC_INSTANCE)
#define N_LOWMC CONCAT(LOWMC, LOWMC_M)
#define SBOX_FUNC CONCAT(sbox_layer, LOWMC_M)
#define SBOX(x) SBOX_FUNC(BLOCK(x, 0))
#include "lowmc_impl.c.i"

#if defined(WITH_ZKBPP)
#undef N_LOWMC
#define N_LOWMC CONCAT(CONCAT(LOWMC, LOWMC_M), store)
#define RECORD_STATE
#include "lowmc_impl.c.i"
#endif

#if defined(WITH_KKW)
#undef N_LOWMC
#undef RECORD_STATE
#undef SBOX
#undef SBOX_FUNC
#define SBOX_FUNC CONCAT(CONCAT(sbox_layer, LOWMC_M), aux)
#define SBOX(x, tapes) SBOX_FUNC(BLOCK(x, 0), tapes)
#define N_LOWMC CONCAT(CONCAT(LOWMC, LOWMC_M), compute_aux)
#define PICNIC2_AUX_COMPUTATION
#include "lowmc_impl.c.i"
#endif

#undef LOWMC_M
#undef N_LOWMC
#undef RECORD_STATE
#undef PICNIC2_AUX_COMPUTATION
#undef SBOX
#undef SBOX_FUNC
#endif

// vim: ft=c
