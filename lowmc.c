/*
 *  This file is part of the optimized implementation of the Picnic signature scheme.
 *  See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "io.h"
#include "lowmc.h"
#include "lowmc_pars.h"
#include "mzd_additional.h"

#if defined(WITH_OPT)
#include "simd.h"
#endif

#if !defined(_MSC_VER)
#include <stdalign.h>
#endif
#include <string.h>

static uint64_t sbox_layer_10_bitsliced_uint64(uint64_t in) {
  // a, b, c
  const uint64_t x0s = (in & MASK_X0I) << 2;
  const uint64_t x1s = (in & MASK_X1I) << 1;
  const uint64_t x2m = in & MASK_X2I;

  // (b & c) ^ a
  const uint64_t t0 = (x1s & x2m) ^ x0s;
  // (c & a) ^ a ^ b
  const uint64_t t1 = (x0s & x2m) ^ x0s ^ x1s;
  // (a & b) ^ a ^ b ^c
  const uint64_t t2 = (x0s & x1s) ^ x0s ^ x1s ^ x2m;

  return (in & MASK_MASK) ^ (t0 >> 2) ^ (t1 >> 1) ^ t2;
}

/**
 * S-box for m = 10
 */
static void sbox_layer_10_uint64(mzd_local_t* x) {
  uint64_t* d = &FIRST_ROW(x)[x->width - 1];
  *d          = sbox_layer_10_bitsliced_uint64(*d);
}

#if defined(WITH_LOWMC_M1)
static uint64_t sbox_layer_1_bitsliced_uint64(uint64_t in) {
  // a, b, c
  const uint64_t x0s = (in & MASK_X0I_1) << 2;
  const uint64_t x1s = (in & MASK_X1I_1) << 1;
  const uint64_t x2m = in & MASK_X2I_1;

  // (b & c) ^ a
  const uint64_t t0 = (x1s & x2m) ^ x0s;
  // (c & a) ^ a ^ b
  const uint64_t t1 = (x0s & x2m) ^ x0s ^ x1s;
  // (a & b) ^ a ^ b ^c
  const uint64_t t2 = (x0s & x1s) ^ x0s ^ x1s ^ x2m;

  return (in & MASK_MASK_1) ^ (t0 >> 2) ^ (t1 >> 1) ^ t2;
}

/**
 * S-box for m = 1
 */
static void sbox_layer_1_uint64(mzd_local_t* x) {
  uint64_t* d = &FIRST_ROW(x)[x->width - 1];
  *d          = sbox_layer_1_bitsliced_uint64(*d);
}
#endif

#if defined(WITH_LOWMC_128_128_20)
#include "lowmc_128_128_20.h"
#endif
#if defined(WITH_LOWMC_192_192_30)
#include "lowmc_192_192_30.h"
#endif
#if defined(WITH_LOWMC_256_256_38)
#include "lowmc_256_256_38.h"
#endif
#if defined(WITH_LOWMC_128_128_182)
#include "lowmc_128_128_182.h"
#endif
#if defined(WITH_LOWMC_192_192_284)
#include "lowmc_192_192_284.h"
#endif
#if defined(WITH_LOWMC_256_256_363)
#include "lowmc_256_256_363.h"
#endif

// uint64 based implementation
#define XOR mzd_xor_uint64
#define MUL SELECT_V_VL(mzd_mul_v_uint64, mzd_mul_vl_uint64)
#define ADDMUL SELECT_V_VL(mzd_addmul_v_uint64, mzd_addmul_vl_uint64)
#define XOR_MC mzd_xor_uint64
#define MUL_MC SELECT_V_VL(mzd_mul_v_uint64, mzd_mul_vl_uint64)
#define SHUFFLE mzd_shuffle

#define MUL_R_1  mzd_addmul_v_uint64_3
#define MUL_R_10 mzd_addmul_v_uint64_30
#define MUL_Z_1  mzd_mul_v_parity_uint64_128_3
#define MUL_Z_10 mzd_mul_v_parity_uint64_128_30

#define LOWMC_N LOWMC_L1_N
#define LOWMC_R_10 LOWMC_L1_R
#define LOWMC_R_1 LOWMC_L1_1_R
#if defined(WITH_LOWMC_128_128_20)
#define LOWMC_INSTANCE_10 (&lowmc_128_128_20)
#endif
#if defined(WITH_LOWMC_128_128_182)
#define LOWMC_INSTANCE_1 (&lowmc_128_128_182)
#endif

#define LOWMC lowmc_uint64_128
#include "lowmc.c.i"

#undef MUL_Z_1
#undef MUL_Z_10
#define MUL_Z_1  mzd_mul_v_parity_uint64_192_3
#define MUL_Z_10 mzd_mul_v_parity_uint64_192_30

#undef LOWMC
#undef LOWMC_N
#undef LOWMC_R_10
#undef LOWMC_R_1
#undef LOWMC_INSTANCE_10
#undef LOWMC_INSTANCE_1

#define LOWMC_N LOWMC_L3_N
#define LOWMC_R_10 LOWMC_L3_R
#define LOWMC_R_1 LOWMC_L3_1_R
#if defined(WITH_LOWMC_192_192_30)
#define LOWMC_INSTANCE_10 (&lowmc_192_192_30)
#endif
#if defined(WITH_LOWMC_192_192_284)
#define LOWMC_INSTANCE_1 (&lowmc_192_192_284)
#endif
#define LOWMC lowmc_uint64_192
#include "lowmc.c.i"

#undef MUL_Z_1
#undef MUL_Z_10
#define MUL_Z_1  mzd_mul_v_parity_uint64_256_3
#define MUL_Z_10 mzd_mul_v_parity_uint64_256_30

#undef LOWMC
#undef LOWMC_N
#undef LOWMC_R_10
#undef LOWMC_R_1
#undef LOWMC_INSTANCE_10
#undef LOWMC_INSTANCE_1

#define LOWMC_N LOWMC_L5_N
#define LOWMC_R_10 LOWMC_L5_R
#define LOWMC_R_1 LOWMC_L5_1_R
#if defined(WITH_LOWMC_256_256_38)
#define LOWMC_INSTANCE_10 (&lowmc_256_256_38)
#endif
#if defined(WITH_LOWMC_256_256_363)
#define LOWMC_INSTANCE_1 (&lowmc_256_256_363)
#endif
#define LOWMC lowmc_uint64_256
#include "lowmc.c.i"

#undef LOWMC
#undef LOWMC_N
#undef LOWMC_R_10
#undef LOWMC_R_1
#undef LOWMC_INSTANCE_10
#undef LOWMC_INSTANCE_1

#if defined(WITH_OPT)
#if defined(WITH_SSE2)
#undef XOR_MC
#undef MUL_MC
#define XOR_MC mzd_xor_sse
#define MUL_MC SELECT_V_VL(mzd_mul_v_sse, mzd_mul_vl_sse)
#define FN_ATTR ATTR_TARGET("sse2")

// L1 using SSE2
#undef XOR
#undef MUL
#undef ADDMUL
#define XOR mzd_xor_sse_128
#define MUL SELECT_V_VL(mzd_mul_v_sse_128, mzd_mul_vl_sse_128)
#define ADDMUL SELECT_V_VL(mzd_addmul_v_sse_128, mzd_addmul_vl_sse_128)

#undef LOWMC_INSTANCE_1
#undef LOWMC_INSTANCE_10
#undef LOWMC_N
#undef LOWMC_R_1
#undef LOWMC_R_10
#if defined(WITH_LOWMC_128_128_20)
#define LOWMC_INSTANCE_10 (&lowmc_128_128_20)
#endif
#if defined(WITH_LOWMC_128_128_182)
#define LOWMC_INSTANCE_1 (&lowmc_128_128_182)
#endif
#define LOWMC_N LOWMC_L1_N
#define LOWMC_R_10 LOWMC_L1_R
#define LOWMC_R_1 LOWMC_L1_1_R

#undef MUL_R_1
#undef MUL_R_10
#undef MUL_Z_1
#undef MUL_Z_10
#define MUL_R_1  mzd_addmul_v_sse_3_128
#define MUL_R_10 mzd_addmul_v_sse_30_128
#define MUL_Z_1  mzd_mul_v_parity_uint64_128_3
#define MUL_Z_10 mzd_mul_v_parity_uint64_128_30

#undef LOWMC
#define LOWMC lowmc_sse_128
#include "lowmc.c.i"

// L3 using SSE2
#undef XOR
#undef MUL
#undef ADDMUL
#define XOR mzd_xor_sse_256
#define MUL SELECT_V_VL(mzd_mul_v_sse_192, mzd_mul_vl_sse_192)
#define ADDMUL SELECT_V_VL(mzd_addmul_v_sse_192, mzd_addmul_vl_sse_192)

#undef LOWMC_INSTANCE_1
#undef LOWMC_INSTANCE_10
#undef LOWMC_N
#undef LOWMC_R_1
#undef LOWMC_R_10
#if defined(WITH_LOWMC_192_192_30)
#define LOWMC_INSTANCE_10 (&lowmc_192_192_30)
#endif
#if defined(WITH_LOWMC_192_192_284)
#define LOWMC_INSTANCE_1 (&lowmc_192_192_284)
#endif
#define LOWMC_N LOWMC_L3_N
#define LOWMC_R_10 LOWMC_L3_R
#define LOWMC_R_1 LOWMC_L3_1_R

#undef MUL_R_1
#undef MUL_R_10
#undef MUL_Z_1
#undef MUL_Z_10
#define MUL_R_1  mzd_addmul_v_sse_3_192
#define MUL_R_10 mzd_addmul_v_sse_30_192
#define MUL_Z_1  mzd_mul_v_parity_uint64_192_3
#define MUL_Z_10 mzd_mul_v_parity_uint64_192_30

#undef LOWMC
#define LOWMC lowmc_sse_192
#include "lowmc.c.i"

// L5 using SSE2
#undef MUL
#undef ADDMUL
#define MUL SELECT_V_VL(mzd_mul_v_sse_256, mzd_mul_vl_sse_256)
#define ADDMUL SELECT_V_VL(mzd_addmul_v_sse_256, mzd_addmul_vl_sse_256)

#undef LOWMC_INSTANCE_1
#undef LOWMC_INSTANCE_10
#undef LOWMC_N
#undef LOWMC_R_1
#undef LOWMC_R_10
#if defined(WITH_LOWMC_256_256_38)
#define LOWMC_INSTANCE_10 (&lowmc_256_256_38)
#endif
#if defined(WITH_LOWMC_256_256_363)
#define LOWMC_INSTANCE_1 (&lowmc_256_256_363)
#endif
#define LOWMC_N LOWMC_L5_N
#define LOWMC_R_10 LOWMC_L5_R
#define LOWMC_R_1 LOWMC_L5_1_R

#undef MUL_R_1
#undef MUL_R_10
#undef MUL_Z_1
#undef MUL_Z_10
#define MUL_R_1  mzd_addmul_v_sse_3_256
#define MUL_R_10 mzd_addmul_v_sse_30_256
#define MUL_Z_1  mzd_mul_v_parity_uint64_256_3
#define MUL_Z_10 mzd_mul_v_parity_uint64_256_30

#undef LOWMC
#define LOWMC lowmc_sse_256
#include "lowmc.c.i"

#undef FN_ATTR
#endif

#if defined(WITH_SSE2) && defined(WITH_POPCNT)
#undef XOR_MC
#undef MUL_MC
#define XOR_MC mzd_xor_sse
#define MUL_MC SELECT_V_VL(mzd_mul_v_sse, mzd_mul_vl_sse)
#define FN_ATTR ATTR_TARGET("sse2,popcnt")

// L1 using SSE2
#undef XOR
#undef MUL
#undef ADDMUL
#define XOR mzd_xor_sse_128
#define MUL SELECT_V_VL(mzd_mul_v_sse_128, mzd_mul_vl_sse_128)
#define ADDMUL SELECT_V_VL(mzd_addmul_v_sse_128, mzd_addmul_vl_sse_128)

#undef LOWMC_INSTANCE_1
#undef LOWMC_INSTANCE_10
#undef LOWMC_N
#undef LOWMC_R_1
#undef LOWMC_R_10
#if defined(WITH_LOWMC_128_128_20)
#define LOWMC_INSTANCE_10 (&lowmc_128_128_20)
#endif
#if defined(WITH_LOWMC_128_128_182)
#define LOWMC_INSTANCE_1 (&lowmc_128_128_182)
#endif
#define LOWMC_N LOWMC_L1_N
#define LOWMC_R_10 LOWMC_L1_R
#define LOWMC_R_1 LOWMC_L1_1_R

#undef MUL_R_1
#undef MUL_R_10
#undef MUL_Z_1
#undef MUL_Z_10
#define MUL_R_1  mzd_addmul_v_sse_3_128
#define MUL_R_10 mzd_addmul_v_sse_30_128
#define MUL_Z_1  mzd_mul_v_parity_popcnt_128_3
#define MUL_Z_10 mzd_mul_v_parity_popcnt_128_30

#undef LOWMC
#define LOWMC lowmc_sse_popcnt_128
#include "lowmc.c.i"

// L3 using SSE2
#undef XOR
#undef MUL
#undef ADDMUL
#define XOR mzd_xor_sse_256
#define MUL SELECT_V_VL(mzd_mul_v_sse_192, mzd_mul_vl_sse_192)
#define ADDMUL SELECT_V_VL(mzd_addmul_v_sse_192, mzd_addmul_vl_sse_192)

#undef LOWMC_INSTANCE_1
#undef LOWMC_INSTANCE_10
#undef LOWMC_N
#undef LOWMC_R_1
#undef LOWMC_R_10
#if defined(WITH_LOWMC_192_192_30)
#define LOWMC_INSTANCE_10 (&lowmc_192_192_30)
#endif
#if defined(WITH_LOWMC_192_192_284)
#define LOWMC_INSTANCE_1 (&lowmc_192_192_284)
#endif
#define LOWMC_N LOWMC_L3_N
#define LOWMC_R_10 LOWMC_L3_R
#define LOWMC_R_1 LOWMC_L3_1_R

#undef MUL_R_1
#undef MUL_R_10
#undef MUL_Z_1
#undef MUL_Z_10
#define MUL_R_1  mzd_addmul_v_sse_3_192
#define MUL_R_10 mzd_addmul_v_sse_30_192
#define MUL_Z_1  mzd_mul_v_parity_popcnt_192_3
#define MUL_Z_10 mzd_mul_v_parity_popcnt_192_30

#undef LOWMC
#define LOWMC lowmc_sse_popcnt_192
#include "lowmc.c.i"

// L5 using SSE2
#undef MUL
#undef ADDMUL
#define MUL SELECT_V_VL(mzd_mul_v_sse_256, mzd_mul_vl_sse_256)
#define ADDMUL SELECT_V_VL(mzd_addmul_v_sse_256, mzd_addmul_vl_sse_256)

#undef LOWMC_INSTANCE_1
#undef LOWMC_INSTANCE_10
#undef LOWMC_N
#undef LOWMC_R_1
#undef LOWMC_R_10
#if defined(WITH_LOWMC_256_256_38)
#define LOWMC_INSTANCE_10 (&lowmc_256_256_38)
#endif
#if defined(WITH_LOWMC_256_256_363)
#define LOWMC_INSTANCE_1 (&lowmc_256_256_363)
#endif
#define LOWMC_N LOWMC_L5_N
#define LOWMC_R_10 LOWMC_L5_R
#define LOWMC_R_1 LOWMC_L5_1_R

#undef MUL_R_1
#undef MUL_R_10
#undef MUL_Z_1
#undef MUL_Z_10
#define MUL_R_1  mzd_addmul_v_sse_3_256
#define MUL_R_10 mzd_addmul_v_sse_30_256
#define MUL_Z_1  mzd_mul_v_parity_popcnt_256_3
#define MUL_Z_10 mzd_mul_v_parity_popcnt_256_30

#undef LOWMC
#define LOWMC lowmc_sse_popcnt_256
#include "lowmc.c.i"

#undef FN_ATTR
#endif

#if defined(WITH_AVX2)
#undef XOR_MC
#undef MUL_MC
#undef SHUFFLE
#define XOR_MC mzd_xor_avx
#define MUL_MC SELECT_V_VL(mzd_mul_v_avx, mzd_mul_vl_avx)
#define SHUFFLE mzd_shuffle_pext
#define FN_ATTR ATTR_TARGET("avx2,bmi2")

// L1 using AVX2
#undef XOR
#undef MUL
#undef ADDMUL
#define XOR mzd_xor_sse_128
#define MUL SELECT_V_VL(mzd_mul_v_avx_128, mzd_mul_vl_avx_128)
#define ADDMUL SELECT_V_VL(mzd_addmul_v_avx_128, mzd_addmul_vl_avx_128)

#undef LOWMC_INSTANCE_1
#undef LOWMC_INSTANCE_10
#undef LOWMC_N
#undef LOWMC_R_1
#undef LOWMC_R_10
#if defined(WITH_LOWMC_128_128_20)
#define LOWMC_INSTANCE_10 (&lowmc_128_128_20)
#endif
#if defined(WITH_LOWMC_128_128_182)
#define LOWMC_INSTANCE_1 (&lowmc_128_128_182)
#endif
#define LOWMC_N LOWMC_L1_N
#define LOWMC_R_10 LOWMC_L1_R
#define LOWMC_R_1 LOWMC_L1_1_R

#undef MUL_R_1
#undef MUL_R_10
#undef MUL_Z_1
#undef MUL_Z_10
#define MUL_R_1  mzd_addmul_v_sse_3_128
#define MUL_R_10 mzd_addmul_v_avx_30_128
#define MUL_Z_1  mzd_mul_v_parity_uint64_128_3
#define MUL_Z_10 mzd_mul_v_parity_uint64_128_30

#undef LOWMC
#define LOWMC lowmc_avx_128
#include "lowmc.c.i"

// L3 using AVX2
#undef XOR
#undef MUL
#undef ADDMUL
#define XOR mzd_xor_avx_256
#define MUL SELECT_V_VL(mzd_mul_v_avx_192, mzd_mul_vl_avx_192)
#define ADDMUL SELECT_V_VL(mzd_addmul_v_avx_192, mzd_addmul_vl_avx_192)

#undef LOWMC_INSTANCE_1
#undef LOWMC_INSTANCE_10
#undef LOWMC_N
#undef LOWMC_R_1
#undef LOWMC_R_10
#if defined(WITH_LOWMC_192_192_30)
#define LOWMC_INSTANCE_10 (&lowmc_192_192_30)
#endif
#if defined(WITH_LOWMC_192_192_284)
#define LOWMC_INSTANCE_1 (&lowmc_192_192_284)
#endif
#define LOWMC_N LOWMC_L3_N
#define LOWMC_R_10 LOWMC_L3_R
#define LOWMC_R_1 LOWMC_L3_1_R

#undef MUL_R_1
#undef MUL_R_10
#undef MUL_Z_1
#undef MUL_Z_10
#define MUL_R_1  mzd_addmul_v_avx_3_192
#define MUL_R_10 mzd_addmul_v_avx_30_192
#define MUL_Z_1  mzd_mul_v_parity_uint64_192_3
#define MUL_Z_10 mzd_mul_v_parity_uint64_192_30

#undef LOWMC
#define LOWMC lowmc_avx_192
#include "lowmc.c.i"

// L5 using AVX2
#undef MUL
#undef ADDMUL
#define MUL SELECT_V_VL(mzd_mul_v_avx_256, mzd_mul_vl_avx_256)
#define ADDMUL SELECT_V_VL(mzd_addmul_v_avx_256, mzd_addmul_vl_avx_256)

#undef LOWMC_INSTANCE_1
#undef LOWMC_INSTANCE_10
#undef LOWMC_N
#undef LOWMC_R_1
#undef LOWMC_R_10
#if defined(WITH_LOWMC_256_256_38)
#define LOWMC_INSTANCE_10 (&lowmc_256_256_38)
#endif
#if defined(WITH_LOWMC_256_256_363)
#define LOWMC_INSTANCE_1 (&lowmc_256_256_363)
#endif
#define LOWMC_N LOWMC_L5_N
#define LOWMC_R_10 LOWMC_L5_R
#define LOWMC_R_1 LOWMC_L5_1_R

#undef MUL_R_1
#undef MUL_R_10
#undef MUL_Z_1
#undef MUL_Z_10
#define MUL_R_1  mzd_addmul_v_avx_3_256
#define MUL_R_10 mzd_addmul_v_avx_30_256
#define MUL_Z_1  mzd_mul_v_parity_uint64_256_3
#define MUL_Z_10 mzd_mul_v_parity_uint64_256_30

#undef LOWMC
#define LOWMC lowmc_avx_256
#include "lowmc.c.i"

#undef FN_ATTR

#if defined(WITH_POPCNT)
#define FN_ATTR ATTR_TARGET("avx2,bmi2,popcnt")

// L1 using AVX2
#undef XOR
#undef MUL
#undef ADDMUL
#define XOR mzd_xor_sse_128
#define MUL SELECT_V_VL(mzd_mul_v_avx_128, mzd_mul_vl_avx_128)
#define ADDMUL SELECT_V_VL(mzd_addmul_v_avx_128, mzd_addmul_vl_avx_128)

#undef LOWMC_INSTANCE_1
#undef LOWMC_INSTANCE_10
#undef LOWMC_N
#undef LOWMC_R_1
#undef LOWMC_R_10
#if defined(WITH_LOWMC_128_128_20)
#define LOWMC_INSTANCE_10 (&lowmc_128_128_20)
#endif
#if defined(WITH_LOWMC_128_128_182)
#define LOWMC_INSTANCE_1 (&lowmc_128_128_182)
#endif
#define LOWMC_N LOWMC_L1_N
#define LOWMC_R_10 LOWMC_L1_R
#define LOWMC_R_1 LOWMC_L1_1_R

#undef MUL_R_1
#undef MUL_R_10
#undef MUL_Z_1
#undef MUL_Z_10
#define MUL_R_1  mzd_addmul_v_sse_3_128
#define MUL_R_10 mzd_addmul_v_avx_30_128
#define MUL_Z_1  mzd_mul_v_parity_popcnt_128_3
#define MUL_Z_10 mzd_mul_v_parity_popcnt_128_30

#undef LOWMC
#define LOWMC lowmc_avx_popcnt_128
#include "lowmc.c.i"

// L3 using AVX2
#undef XOR
#undef MUL
#undef ADDMUL
#define XOR mzd_xor_avx_256
#define MUL SELECT_V_VL(mzd_mul_v_avx_192, mzd_mul_vl_avx_192)
#define ADDMUL SELECT_V_VL(mzd_addmul_v_avx_192, mzd_addmul_vl_avx_192)

#undef LOWMC_INSTANCE_1
#undef LOWMC_INSTANCE_10
#undef LOWMC_N
#undef LOWMC_R_1
#undef LOWMC_R_10
#if defined(WITH_LOWMC_192_192_30)
#define LOWMC_INSTANCE_10 (&lowmc_192_192_30)
#endif
#if defined(WITH_LOWMC_192_192_284)
#define LOWMC_INSTANCE_1 (&lowmc_192_192_284)
#endif
#define LOWMC_N LOWMC_L3_N
#define LOWMC_R_10 LOWMC_L3_R
#define LOWMC_R_1 LOWMC_L3_1_R

#undef MUL_R_1
#undef MUL_R_10
#undef MUL_Z_1
#undef MUL_Z_10
#define MUL_R_1  mzd_addmul_v_avx_3_192
#define MUL_R_10 mzd_addmul_v_avx_30_192
#define MUL_Z_1  mzd_mul_v_parity_popcnt_192_3
#define MUL_Z_10 mzd_mul_v_parity_popcnt_192_30

#undef LOWMC
#define LOWMC lowmc_avx_popcnt_192
#include "lowmc.c.i"

// L5 using AVX2
#undef MUL
#undef ADDMUL
#define MUL SELECT_V_VL(mzd_mul_v_avx_256, mzd_mul_vl_avx_256)
#define ADDMUL SELECT_V_VL(mzd_addmul_v_avx_256, mzd_addmul_vl_avx_256)

#undef LOWMC_INSTANCE_1
#undef LOWMC_INSTANCE_10
#undef LOWMC_N
#undef LOWMC_R_1
#undef LOWMC_R_10
#if defined(WITH_LOWMC_256_256_38)
#define LOWMC_INSTANCE_10 (&lowmc_256_256_38)
#endif
#if defined(WITH_LOWMC_256_256_363)
#define LOWMC_INSTANCE_1 (&lowmc_256_256_363)
#endif
#define LOWMC_N LOWMC_L5_N
#define LOWMC_R_10 LOWMC_L5_R
#define LOWMC_R_1 LOWMC_L5_1_R

#undef MUL_R_1
#undef MUL_R_10
#undef MUL_Z_1
#undef MUL_Z_10
#define MUL_R_1  mzd_addmul_v_avx_3_256
#define MUL_R_10 mzd_addmul_v_avx_30_256
#define MUL_Z_1  mzd_mul_v_parity_popcnt_256_3
#define MUL_Z_10 mzd_mul_v_parity_popcnt_256_30

#undef LOWMC
#define LOWMC lowmc_avx_popcnt_256
#include "lowmc.c.i"

#undef FN_ATTR
#endif

#undef SHUFFLE
#define SHUFFLE mzd_shuffle
#endif

#if defined(WITH_NEON)
#undef XOR_MC
#undef MUL_MC
#define XOR_MC mzd_xor_neon
#define MUL_MC SELECT_V_VL(mzd_mul_v_neon, mzd_mul_vl_neon)

// L1 using NEON
#undef XOR
#undef MUL
#undef ADDMUL
#define XOR mzd_xor_neon_128
#define MUL SELECT_V_VL(mzd_mul_v_neon_128, mzd_mul_vl_neon_128)
#define ADDMUL SELECT_V_VL(mzd_addmul_v_neon_128, mzd_addmul_vl_neon_128)

#undef LOWMC_INSTANCE_1
#undef LOWMC_INSTANCE_10
#undef LOWMC_N
#undef LOWMC_R_1
#undef LOWMC_R_10
#if defined(WITH_LOWMC_128_128_20)
#define LOWMC_INSTANCE_10 (&lowmc_128_128_20)
#endif
#if defined(WITH_LOWMC_128_128_182)
#define LOWMC_INSTANCE_1 (&lowmc_128_128_182)
#endif
#define LOWMC_N LOWMC_L1_N
#define LOWMC_R_10 LOWMC_L1_R
#define LOWMC_R_1 LOWMC_L1_1_R

#undef MUL_R_1
#undef MUL_R_10
#undef MUL_Z_1
#undef MUL_Z_10
#define MUL_R_1  mzd_addmul_v_neon_3_128
#define MUL_R_10 mzd_addmul_v_neon_30_128
#define MUL_Z_1  mzd_mul_v_parity_uint64_128_3
#define MUL_Z_10 mzd_mul_v_parity_uint64_128_30

#undef LOWMC
#define LOWMC lowmc_neon_128
#include "lowmc.c.i"

// L3 using NEON
#undef XOR
#undef MUL
#undef ADDMUL
#define XOR mzd_xor_neon_256
#define MUL SELECT_V_VL(mzd_mul_v_neon_192, mzd_mul_vl_neon_192)
#define ADDMUL SELECT_V_VL(mzd_addmul_v_neon_192, mzd_addmul_vl_neon_192)

#undef LOWMC_INSTANCE_1
#undef LOWMC_INSTANCE_10
#undef LOWMC_N
#undef LOWMC_R_1
#undef LOWMC_R_10
#if defined(WITH_LOWMC_192_192_30)
#define LOWMC_INSTANCE_10 (&lowmc_192_192_30)
#endif
#if defined(WITH_LOWMC_192_192_284)
#define LOWMC_INSTANCE_1 (&lowmc_192_192_284)
#endif
#define LOWMC_N LOWMC_L3_N
#define LOWMC_R_10 LOWMC_L3_R
#define LOWMC_R_1 LOWMC_L3_1_R

#undef MUL_R_1
#undef MUL_R_10
#undef MUL_Z_1
#undef MUL_Z_10
#define MUL_R_1  mzd_addmul_v_neon_3_192
#define MUL_R_10 mzd_addmul_v_neon_30_192
#define MUL_Z_1  mzd_mul_v_parity_uint64_192_3
#define MUL_Z_10 mzd_mul_v_parity_uint64_192_30

#undef LOWMC
#define LOWMC lowmc_neon_192
#include "lowmc.c.i"

// L5 using NEON
#undef MUL
#undef ADDMUL
#define MUL SELECT_V_VL(mzd_mul_v_neon_256, mzd_mul_vl_neon_256)
#define ADDMUL SELECT_V_VL(mzd_addmul_v_neon_256, mzd_addmul_vl_neon_256)

#undef LOWMC_INSTANCE_1
#undef LOWMC_INSTANCE_10
#undef LOWMC_N
#undef LOWMC_R_1
#undef LOWMC_R_10
#if defined(WITH_LOWMC_256_256_38)
#define LOWMC_INSTANCE_10 (&lowmc_256_256_38)
#endif
#if defined(WITH_LOWMC_256_256_363)
#define LOWMC_INSTANCE_1 (&lowmc_256_256_363)
#endif
#define LOWMC_N LOWMC_L5_N
#define LOWMC_R_10 LOWMC_L5_R
#define LOWMC_R_1 LOWMC_L5_1_R

#undef MUL_R_1
#undef MUL_R_10
#undef MUL_Z_1
#undef MUL_Z_10
#define MUL_R_1  mzd_addmul_v_neon_3_256
#define MUL_R_10 mzd_addmul_v_neon_30_256
#define MUL_Z_1  mzd_mul_v_parity_uint64_256_3
#define MUL_Z_10 mzd_mul_v_parity_uint64_256_30

#undef LOWMC
#define LOWMC lowmc_neon_256
#include "lowmc.c.i"

#endif
#endif

lowmc_implementation_f lowmc_get_implementation(const lowmc_t* lowmc) {
#if defined(WITH_OPT)
#if defined(WITH_AVX2)
  if (CPU_SUPPORTS_AVX2) {
    if (lowmc->m == 10) {
#if defined(WITH_POPCNT)
      if (CPU_SUPPORTS_POPCNT) {
        switch (lowmc->n) {
#if defined(WITH_LOWMC_128_128_20)
        case 128:
          return lowmc_avx_popcnt_128_10;
#endif
#if defined(WITH_LOWMC_192_192_30)
        case 192:
          return lowmc_avx_popcnt_192_10;
#endif
#if defined(WITH_LOWMC_256_256_38)
        case 256:
          return lowmc_avx_popcnt_256_10;
#endif
        }
      }
#endif
      switch (lowmc->n) {
#if defined(WITH_LOWMC_128_128_20)
      case 128:
        return lowmc_avx_128_10;
#endif
#if defined(WITH_LOWMC_192_192_30)
      case 192:
        return lowmc_avx_192_10;
#endif
#if defined(WITH_LOWMC_256_256_38)
      case 256:
        return lowmc_avx_256_10;
#endif
      }
    }
#if defined(WITH_LOWMC_M1)
    if (lowmc->m == 1) {
#if defined(WITH_POPCNT)
      if (CPU_SUPPORTS_POPCNT) {
        switch (lowmc->n) {
#if defined(WITH_LOWMC_128_128_182)
        case 128:
          return lowmc_avx_popcnt_128_1;
#endif
#if defined(WITH_LOWMC_192_192_284)
        case 192:
          return lowmc_avx_popcnt_192_1;
#endif
#if defined(WITH_LOWMC_256_256_363)
        case 256:
          return lowmc_avx_popcnt_256_1;
#endif
        }
      }
#endif
      switch (lowmc->n) {
#if defined(WITH_LOWMC_128_128_182)
      case 128:
        return lowmc_avx_128_1;
#endif
#if defined(WITH_LOWMC_192_192_284)
      case 192:
        return lowmc_avx_192_1;
#endif
#if defined(WITH_LOWMC_256_256_363)
      case 256:
        return lowmc_avx_256_1;
#endif
      }
    }
#endif
  }
#endif
#if defined(WITH_SSE2)
  if (CPU_SUPPORTS_SSE2) {
    if (lowmc->m == 10) {
#if defined(WITH_POPCNT)
      if (CPU_SUPPORTS_POPCNT) {
        switch (lowmc->n) {
#if defined(WITH_LOWMC_128_128_20)
        case 128:
          return lowmc_sse_popcnt_128_10;
#endif
#if defined(WITH_LOWMC_192_192_30)
        case 192:
          return lowmc_sse_popcnt_192_10;
#endif
#if defined(WITH_LOWMC_256_256_38)
        case 256:
          return lowmc_sse_popcnt_256_10;
#endif
        }
      }
#endif
      switch (lowmc->n) {
#if defined(WITH_LOWMC_128_128_20)
      case 128:
        return lowmc_sse_128_10;
#endif
#if defined(WITH_LOWMC_192_192_30)
      case 192:
        return lowmc_sse_192_10;
#endif
#if defined(WITH_LOWMC_256_256_38)
      case 256:
        return lowmc_sse_256_10;
#endif
      }
    }
#if defined(WITH_LOWMC_M1)
    if (lowmc->m == 1) {
#if defined(WITH_POPCNT)
      if (CPU_SUPPORTS_POPCNT) {
        switch (lowmc->n) {
#if defined(WITH_LOWMC_128_128_182)
        case 128:
          return lowmc_sse_popcnt_128_1;
#endif
#if defined(WITH_LOWMC_192_192_284)
        case 192:
          return lowmc_sse_popcnt_192_1;
#endif
#if defined(WITH_LOWMC_256_256_363)
        case 256:
          return lowmc_sse_popcnt_256_1;
#endif
        }
      }
#endif
      switch (lowmc->n) {
#if defined(WITH_LOWMC_128_128_182)
      case 128:
        return lowmc_sse_128_1;
#endif
#if defined(WITH_LOWMC_192_192_284)
      case 192:
        return lowmc_sse_192_1;
#endif
#if defined(WITH_LOWMC_256_256_363)
      case 256:
        return lowmc_sse_256_1;
#endif
      }
    }
#endif
  }
#endif
#if defined(WITH_NEON)
  if (CPU_SUPPORTS_NEON) {
    if (lowmc->m == 10) {
      switch (lowmc->n) {
#if defined(WITH_LOWMC_128_128_20)
      case 128:
        return lowmc_neon_128_10;
#endif
#if defined(WITH_LOWMC_192_192_30)
      case 192:
        return lowmc_neon_192_10;
#endif
#if defined(WITH_LOWMC_256_256_38)
      case 256:
        return lowmc_neon_256_10;
#endif
      }
    }
#if defined(WITH_LOWMC_M1)
    if (lowmc->m == 1) {
      switch (lowmc->n) {
#if defined(WITH_LOWMC_128_128_182)
      case 128:
        return lowmc_neon_128_1;
#endif
#if defined(WITH_LOWMC_192_192_284)
      case 192:
        return lowmc_neon_192_1;
#endif
#if defined(WITH_LOWMC_256_256_363)
      case 256:
        return lowmc_neon_256_1;
#endif
      }
    }
#endif
  }
#endif
#endif

  if (lowmc->m == 10) {
    switch (lowmc->n) {
#if defined(WITH_LOWMC_128_128_20)
    case 128:
      return lowmc_uint64_128_10;
#endif
#if defined(WITH_LOWMC_192_192_30)
    case 192:
      return lowmc_uint64_192_10;
#endif
#if defined(WITH_LOWMC_256_256_38)
    case 256:
      return lowmc_uint64_256_10;
#endif
    }
  }

#if defined(WITH_LOWMC_M1)
  if (lowmc->m == 1) {
    switch (lowmc->n) {
#if defined(WITH_LOWMC_128_128_182)
    case 128:
      return lowmc_uint64_128_1;
#endif
#if defined(WITH_LOWMC_192_192_284)
    case 192:
      return lowmc_uint64_192_1;
#endif
#if defined(WITH_LOWMC_256_256_363)
    case 256:
      return lowmc_uint64_256_1;
#endif
    }
  }
#endif

  return NULL;
}

lowmc_store_implementation_f lowmc_store_get_implementation(const lowmc_t* lowmc) {
#if defined(WITH_OPT)
#if defined(WITH_AVX2)
  if (CPU_SUPPORTS_AVX2) {
    if (lowmc->m == 10) {
#if defined(WITH_POPCNT)
      if (CPU_SUPPORTS_POPCNT) {
        switch (lowmc->n) {
#if defined(WITH_LOWMC_128_128_20)
        case 128:
          return lowmc_avx_popcnt_128_store_10;
#endif
#if defined(WITH_LOWMC_192_192_30)
        case 192:
          return lowmc_avx_popcnt_192_store_10;
#endif
#if defined(WITH_LOWMC_256_256_38)
        case 256:
          return lowmc_avx_popcnt_256_store_10;
#endif
        }
      }
#endif
      switch (lowmc->n) {
#if defined(WITH_LOWMC_128_128_20)
      case 128:
        return lowmc_avx_128_store_10;
#endif
#if defined(WITH_LOWMC_192_192_30)
      case 192:
        return lowmc_avx_192_store_10;
#endif
#if defined(WITH_LOWMC_256_256_38)
      case 256:
        return lowmc_avx_256_store_10;
#endif
      }
    }
#if defined(WITH_LOWMC_M1)
    if (lowmc->m == 1) {
#if defined(WITH_POPCNT)
      if (CPU_SUPPORTS_POPCNT) {
        switch (lowmc->n) {
#if defined(WITH_LOWMC_128_128_182)
        case 128:
          return lowmc_avx_popcnt_128_store_1;
#endif
#if defined(WITH_LOWMC_192_192_284)
        case 192:
          return lowmc_avx_popcnt_192_store_1;
#endif
#if defined(WITH_LOWMC_256_256_363)
        case 256:
          return lowmc_avx_popcnt_256_store_1;
#endif
        }
      }
#endif
      switch (lowmc->n) {
#if defined(WITH_LOWMC_128_128_182)
      case 128:
        return lowmc_avx_128_store_1;
#endif
#if defined(WITH_LOWMC_192_192_284)
      case 192:
        return lowmc_avx_192_store_1;
#endif
#if defined(WITH_LOWMC_256_256_363)
      case 256:
        return lowmc_avx_256_store_1;
#endif
      }
    }
#endif
  }
#endif
#if defined(WITH_SSE2)
  if (CPU_SUPPORTS_SSE2) {
    if (lowmc->m == 10) {
#if defined(WITH_POPCNT)
      if (CPU_SUPPORTS_POPCNT) {
        switch (lowmc->n) {
#if defined(WITH_LOWMC_128_128_20)
        case 128:
          return lowmc_sse_popcnt_128_store_10;
#endif
#if defined(WITH_LOWMC_192_192_30)
        case 192:
          return lowmc_sse_popcnt_192_store_10;
#endif
#if defined(WITH_LOWMC_256_256_38)
        case 256:
          return lowmc_sse_popcnt_256_store_10;
#endif
        }
      }
#endif
      switch (lowmc->n) {
#if defined(WITH_LOWMC_128_128_20)
      case 128:
        return lowmc_sse_128_store_10;
#endif
#if defined(WITH_LOWMC_192_192_30)
      case 192:
        return lowmc_sse_192_store_10;
#endif
#if defined(WITH_LOWMC_256_256_38)
      case 256:
        return lowmc_sse_256_store_10;
#endif
      }
    }
#if defined(WITH_LOWMC_M1)
    if (lowmc->m == 1) {
#if defined(WITH_POPCNT)
      if (CPU_SUPPORTS_POPCNT) {
        switch (lowmc->n) {
#if defined(WITH_LOWMC_128_128_182)
        case 128:
          return lowmc_sse_popcnt_128_store_1;
#endif
#if defined(WITH_LOWMC_192_192_284)
        case 192:
          return lowmc_sse_popcnt_192_store_1;
#endif
#if defined(WITH_LOWMC_256_256_363)
        case 256:
          return lowmc_sse_popcnt_256_store_1;
#endif
        }
      }
#endif
      switch (lowmc->n) {
#if defined(WITH_LOWMC_128_128_182)
      case 128:
        return lowmc_sse_128_store_1;
#endif
#if defined(WITH_LOWMC_192_192_284)
      case 192:
        return lowmc_sse_192_store_1;
#endif
#if defined(WITH_LOWMC_256_256_363)
      case 256:
        return lowmc_sse_256_store_1;
#endif
      }
    }
#endif
  }
#endif
#if defined(WITH_NEON)
  if (CPU_SUPPORTS_NEON) {
    if (lowmc->m == 10) {
      switch (lowmc->n) {
#if defined(WITH_LOWMC_128_128_20)
      case 128:
        return lowmc_neon_128_store_10;
#endif
#if defined(WITH_LOWMC_192_192_30)
      case 192:
        return lowmc_neon_192_store_10;
#endif
#if defined(WITH_LOWMC_256_256_38)
      case 256:
        return lowmc_neon_256_store_10;
#endif
      }
    }
#if defined(WITH_LOWMC_M1)
    if (lowmc->m == 1) {
      switch (lowmc->n) {
#if defined(WITH_LOWMC_128_128_182)
      case 128:
        return lowmc_neon_128_store_1;
#endif
#if defined(WITH_LOWMC_192_192_284)
      case 192:
        return lowmc_neon_192_store_1;
#endif
#if defined(WITH_LOWMC_256_256_363)
      case 256:
        return lowmc_neon_256_store_1;
#endif
      }
    }
#endif
  }
#endif
#endif

  if (lowmc->m == 10) {
    switch (lowmc->n) {
#if defined(WITH_LOWMC_128_128_20)
    case 128:
      return lowmc_uint64_128_store_10;
#endif
#if defined(WITH_LOWMC_192_192_30)
    case 192:
      return lowmc_uint64_192_store_10;
#endif
#if defined(WITH_LOWMC_256_256_38)
    case 256:
      return lowmc_uint64_256_store_10;
#endif
    }
  }

#if defined(WITH_LOWMC_M1)
  if (lowmc->m == 1) {
    switch (lowmc->n) {
#if defined(WITH_LOWMC_128_128_182)
    case 128:
      return lowmc_uint64_128_store_1;
#endif
#if defined(WITH_LOWMC_192_192_284)
    case 192:
      return lowmc_uint64_192_store_1;
#endif
#if defined(WITH_LOWMC_256_256_363)
    case 256:
      return lowmc_uint64_256_store_1;
#endif
    }
  }
#endif

  return NULL;
}
