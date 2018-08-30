/*
 *  This file is part of the optimized implementation of the Picnic signature scheme.
 *  See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */

#define M_FIXED_10
#define N_SIGN CONCAT(SIGN, 10)
#define N_VERIFY CONCAT(VERIFY, 10)
#define MUL_Z MUL_Z_10
#define MUL_A MUL_A_10
#define LOWMC_R LOWMC_R_10
#if defined(LOWMC_INSTANCE_10)
#define LOWMC_INSTANCE LOWMC_INSTANCE_10
#endif
#include "mpc_lowmc_impl.c.i"
#undef MUL_Z
#undef MUL_A
#undef LOWMC_R
#undef LOWMC_INSTANCE
#undef M_FIXED_10

#define M_FIXED_1
#define N_SIGN CONCAT(SIGN, 1)
#define N_VERIFY CONCAT(VERIFY, 1)
#define MUL_Z MUL_Z_1
#define MUL_A MUL_A_1
#define LOWMC_R LOWMC_R_1
#if defined(LOWMC_INSTANCE_1)
#define LOWMC_INSTANCE LOWMC_INSTANCE_1
#endif
#include "mpc_lowmc_impl.c.i"
#undef MUL_Z
#undef MUL_A
#undef LOWMC_R
#undef LOWMC_INSTANCE
#undef M_FIXED_1

#if defined(WITH_CUSTOM_INSTANCES)
#define SBOX_SIGN SIGN_SBOX
#define SBOX_VERIFY VERIFY_SBOX
#define SBOX_ARGS SBOX_NUM_ARGS
#define N_SIGN SIGN
#define N_VERIFY VERIFY
#include "mpc_lowmc_impl.c.i"
#undef SBOX_SIGN
#undef SBOX_VERIFY
#undef SBOX_ARGS
#endif

#undef SIGN
#undef VERIFY

// vim: ft=c
