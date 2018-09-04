/*
 *  This file is part of the optimized implementation of the Picnic signature scheme.
 *  See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */

#define M_FIXED_10
#define SBOX(x, mask) sbox_layer_10_uint64(x)
#define N_LOWMC CONCAT(LOWMC, 10)
#define MUL_R MUL_R_10
#define MUL_Z MUL_Z_10
#define LOWMC_R LOWMC_R_10
#define LOWMC_M 10
#if defined(LOWMC_INSTANCE_10)
#define LOWMC_INSTANCE LOWMC_INSTANCE_10
#endif
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

#define M_FIXED_1
#define SBOX(x, mask) sbox_layer_1_uint64(x)
#define N_LOWMC CONCAT(LOWMC, 1)
#define MUL_R MUL_R_1
#define MUL_Z MUL_Z_1
#define LOWMC_R LOWMC_R_1
#define LOWMC_M 1
#if defined(LOWMC_INSTANCE_1)
#define LOWMC_INSTANCE LOWMC_INSTANCE_1
#endif
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

#if defined(WITH_CUSTOM_INSTANCES)
#define SBOX SBOX_IMPL
#define N_LOWMC LOWMC
#define MUL_R MUL_R_1
#define MUL_Z MUL_Z_1
#define LOWMC_M (lowmc->m)
#include "lowmc_impl.c.i"

#undef N_LOWMC
#define N_LOWMC CONCAT(LOWMC, store)
#define RECORD_STATE
#include "lowmc_impl.c.i"
#undef SBOX
#undef RECORD_STATE
#undef MUL_R
#undef MUL_Z
#undef N_LOWMC
#undef LOWMC_M
#endif

// vim: ft=c
