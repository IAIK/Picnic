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
#include "mpc_lowmc_impl.c.i"
#undef M_FIXED_10

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
