/*
 *  This file is part of the optimized implementation of the Picnic signature scheme.
 *  See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */

#if defined(LOWMC_INSTANCE_10)
#define M_FIXED_10
#define SBOX(x) sbox_layer_10_uint64(&FIRST_ROW(x)[(LOWMC_N / (sizeof(word) * 8)) - 1])
#define N_LOWMC CONCAT(LOWMC, 10)
#define MUL_R MUL_R_10
#define MUL_Z MUL_Z_10
#define LOWMC_R LOWMC_R_10
#define LOWMC_M 10
#define MZD_SHUFFLE CONCAT(SHUFFLE, 30)
#define LOWMC_INSTANCE LOWMC_INSTANCE_10
#include "lowmc_impl.c.i"

#undef N_LOWMC
#define N_LOWMC CONCAT(LOWMC, store_10)
#define RECORD_STATE
#include "lowmc_impl.c.i"
#undef M_FIXED_10
#undef SBOX
#undef RECORD_STATE
#undef MUL_R
#undef MUL_Z
#undef LOWMC_R
#undef LOWMC_INSTANCE
#undef N_LOWMC
#undef LOWMC_M
#undef MZD_SHUFFLE
#endif

#if defined(WITH_LOWMC_M1) && defined(LOWMC_INSTANCE_1)
#define M_FIXED_1
#define SBOX(x) sbox_layer_1_uint64(&FIRST_ROW(x)[(LOWMC_N / (sizeof(word) * 8)) - 1])
#define N_LOWMC CONCAT(LOWMC, 1)
#define MUL_R MUL_R_1
#define MUL_Z MUL_Z_1
#define LOWMC_R LOWMC_R_1
#define LOWMC_M 1
#define MZD_SHUFFLE CONCAT(SHUFFLE, 3)
#define LOWMC_INSTANCE LOWMC_INSTANCE_1
#include "lowmc_impl.c.i"

#undef N_LOWMC
#define N_LOWMC CONCAT(LOWMC, store_1)
#define RECORD_STATE
#include "lowmc_impl.c.i"
#undef RECORD_STATE
#undef M_FIXED_1
#undef SBOX
#undef MUL_R
#undef MUL_Z
#undef LOWMC_R
#undef LOWMC_INSTANCE
#undef N_LOWMC
#undef LOWMC_M
#undef MZD_SHUFFLE
#endif

// vim: ft=c
