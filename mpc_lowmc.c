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

#include "lowmc_pars.h"
#include "mpc.h"
#include "mpc_lowmc.h"
#include "mzd_additional.h"

#if !defined(_MSC_VER)
#include <stdalign.h>
#endif
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#if defined(WITH_OPT)
#include "simd.h"
#endif

#define bitsliced_step_1_uint64_10(sc)                                                             \
  uint64_t r0m[sc];                                                                                \
  uint64_t r0s[sc];                                                                                \
  uint64_t r1m[sc];                                                                                \
  uint64_t r1s[sc];                                                                                \
  uint64_t r2m[sc];                                                                                \
  uint64_t x0s[sc];                                                                                \
  uint64_t x1s[sc];                                                                                \
  uint64_t x2m[sc];                                                                                \
  do {                                                                                             \
    for (unsigned int m = 0; m < (sc); ++m) {                                                      \
      const uint64_t inm   = in[m];                                                                \
      const uint64_t rvecm = rvec[m];                                                              \
                                                                                                   \
      x0s[m] = (inm & MASK_X0I) << 2;                                                              \
      x1s[m] = (inm & MASK_X1I) << 1;                                                              \
      x2m[m] = inm & MASK_X2I;                                                                     \
                                                                                                   \
      r0m[m] = rvecm & MASK_X0I;                                                                   \
      r1m[m] = rvecm & MASK_X1I;                                                                   \
      r2m[m] = rvecm & MASK_X2I;                                                                   \
                                                                                                   \
      r0s[m] = r0m[m] << 2;                                                                        \
      r1s[m] = r1m[m] << 1;                                                                        \
    }                                                                                              \
  } while (0)

#define bitsliced_step_2_uint64_10(sc)                                                             \
  do {                                                                                             \
    for (unsigned int m = 0; m < sc; ++m) {                                                        \
      const uint64_t tmp1 = r2m[m] ^ x0s[m];                                                       \
      const uint64_t tmp2 = x0s[m] ^ x1s[m];                                                       \
      const uint64_t tmp3 = tmp2 ^ r1m[m];                                                         \
      const uint64_t tmp4 = tmp2 ^ r0m[m] ^ x2m[m];                                                \
                                                                                                   \
      in[m] = (in[m] & MASK_MASK) ^ (tmp4) ^ (tmp1 >> 2) ^ (tmp3 >> 1);                            \
    }                                                                                              \
  } while (0)

#define bitsliced_step_1_uint64_1(sc)                                                              \
  uint64_t r0m[sc];                                                                                \
  uint64_t r0s[sc];                                                                                \
  uint64_t r1m[sc];                                                                                \
  uint64_t r1s[sc];                                                                                \
  uint64_t r2m[sc];                                                                                \
  uint64_t x0s[sc];                                                                                \
  uint64_t x1s[sc];                                                                                \
  uint64_t x2m[sc];                                                                                \
  do {                                                                                             \
    for (unsigned int m = 0; m < (sc); ++m) {                                                      \
      const uint64_t inm   = in[m];                                                                \
      const uint64_t rvecm = rvec[m];                                                              \
                                                                                                   \
      x0s[m] = (inm & MASK_X0I_1) << 2;                                                            \
      x1s[m] = (inm & MASK_X1I_1) << 1;                                                            \
      x2m[m] = inm & MASK_X2I_1;                                                                   \
                                                                                                   \
      r0m[m] = rvecm & MASK_X0I_1;                                                                 \
      r1m[m] = rvecm & MASK_X1I_1;                                                                 \
      r2m[m] = rvecm & MASK_X2I_1;                                                                 \
                                                                                                   \
      r0s[m] = r0m[m] << 2;                                                                        \
      r1s[m] = r1m[m] << 1;                                                                        \
    }                                                                                              \
  } while (0)

#define bitsliced_step_2_uint64_1(sc)                                                              \
  do {                                                                                             \
    for (unsigned int m = 0; m < sc; ++m) {                                                        \
      const uint64_t tmp1 = r2m[m] ^ x0s[m];                                                       \
      const uint64_t tmp2 = x0s[m] ^ x1s[m];                                                       \
      const uint64_t tmp3 = tmp2 ^ r1m[m];                                                         \
      const uint64_t tmp4 = tmp2 ^ r0m[m] ^ x2m[m];                                                \
                                                                                                   \
      in[m] = (in[m] & MASK_MASK_1) ^ (tmp4) ^ (tmp1 >> 2) ^ (tmp3 >> 1);                          \
    }                                                                                              \
  } while (0)

static void mpc_sbox_layer_bitsliced_uint64_10(uint64_t* in, view_t* view, uint64_t const* rvec) {
  bitsliced_step_1_uint64_10(SC_PROOF);

  mpc_and_uint64(r0m, x0s, x1s, r2m, view, 0);
  mpc_and_uint64(r2m, x1s, x2m, r1s, view, 1);
  mpc_and_uint64(r1m, x0s, x2m, r0s, view, 2);

  bitsliced_step_2_uint64_10(SC_PROOF);
}

static void mpc_sbox_layer_bitsliced_verify_uint64_10(uint64_t* in, view_t* view,
                                                      uint64_t const* rvec) {
  bitsliced_step_1_uint64_10(SC_VERIFY);

  mpc_and_verify_uint64(r0m, x0s, x1s, r2m, view, MASK_X2I, 0);
  mpc_and_verify_uint64(r2m, x1s, x2m, r1s, view, MASK_X2I, 1);
  mpc_and_verify_uint64(r1m, x0s, x2m, r0s, view, MASK_X2I, 2);

  bitsliced_step_2_uint64_10(SC_VERIFY);
}

static void mpc_sbox_layer_bitsliced_uint64_1(uint64_t* in, view_t* view, uint64_t const* rvec) {
  bitsliced_step_1_uint64_1(SC_PROOF);

  mpc_and_uint64(r0m, x0s, x1s, r2m, view, 0);
  mpc_and_uint64(r2m, x1s, x2m, r1s, view, 1);
  mpc_and_uint64(r1m, x0s, x2m, r0s, view, 2);

  bitsliced_step_2_uint64_1(SC_PROOF);
}

static void mpc_sbox_layer_bitsliced_verify_uint64_1(uint64_t* in, view_t* view,
                                                   uint64_t const* rvec) {
  bitsliced_step_1_uint64_1(SC_VERIFY);

  mpc_and_verify_uint64(r0m, x0s, x1s, r2m, view, MASK_X2I_1, 0);
  mpc_and_verify_uint64(r2m, x1s, x2m, r1s, view, MASK_X2I_1, 1);
  mpc_and_verify_uint64(r1m, x0s, x2m, r0s, view, MASK_X2I_1, 2);

  bitsliced_step_2_uint64_1(SC_VERIFY);
}

#define SBOX_mzd(X, sbox, y, x, views, r, lowmcmask, vars, n, shares)                              \
  CONCAT(SBOX_mzd, X)(sbox, y, x, views, r, lowmcmask, vars, n)

#define SBOX_mzd_5(sbox, y, x, views, r, lowmcmask, vars, n) sbox(y, x, views, r, lowmcmask)
#define SBOX_mzd_6(sbox, y, x, views, r, lowmcmask, vars, n) sbox(y, x, views, r, lowmcmask, vars)

#define SBOX_uint64(X, sbox, y, x, views, r, lowmcmask, vars, n, shares)                           \
  do {                                                                                             \
    uint64_t in[shares];                                                                           \
    for (unsigned int count = 0; count < shares; ++count) {                                        \
      in[count] = CONST_FIRST_ROW(x[count])[(n) / (sizeof(word) * 8) - 1];                         \
    }                                                                                              \
    sbox(in, views, r);                                                                            \
    for (unsigned int count = 0; count < shares; ++count) {                                        \
      memcpy(FIRST_ROW(y[count]), CONST_FIRST_ROW(x[count]),                                       \
             ((n) / (sizeof(word) * 8) - 1) * sizeof(word));                                       \
      FIRST_ROW(y[count])[(n) / (sizeof(word) * 8) - 1] = in[count];                               \
    }                                                                                              \
  } while (0)

#define R_mzd mzd_local_t** r = rvec[i].s
#define R_uint64 const uint64_t* r = rvec[i].t

#define VARS_5(shares, n)
#define VARS_6(shares, n)                                                                          \
  sbox_vars_t vars;                                                                                \
  sbox_vars_init(&vars, n, shares)

#define VARS_FREE_5
#define VARS_FREE_6 sbox_vars_clear(&vars)

// uint64 based implementation
#define XOR mzd_xor_uint64
#define MUL SELECT_V_VL(mzd_mul_v_uint64, mzd_mul_vl_uint64)
#define ADDMUL SELECT_V_VL(mzd_addmul_v_uint64, mzd_addmul_vl_uint64)
#define XOR_MC mzd_xor_uint64
#define MUL_MC SELECT_V_VL(mzd_mul_v_uint64, mzd_mul_vl_uint64)

#define LOWMC_N lowmc->n
#define LOWMC_R_10 lowmc->r
#define LOWMC_R_1 lowmc->r

#define SIGN_SBOX mpc_sbox_layer_bitsliced
#define VERIFY_SBOX mpc_sbox_layer_bitsliced_verify
#define SBOX_NUM_ARGS 6

#undef MUL_R_1
#undef MUL_R_10
#undef MUL_Z_1
#undef MUL_Z_10
#define MUL_R_1  mzd_mul_v_uint64_3
#define MUL_R_10 mzd_mul_v_uint64_30
#define MUL_Z_1  mzd_mul_v_3_popcnt
#define MUL_Z_10 mzd_mul_v_30_popcnt

#define SIGN mpc_lowmc_call
#define VERIFY mpc_lowmc_call_verify
#include "mpc_lowmc.c.i"

#if defined(WITH_OPT)
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

#undef SBOX_NUM_ARGS
#define SBOX_NUM_ARGS 5

#if defined(WITH_SSE2)
#undef XOR_MC
#undef MUL_MC
#define XOR_MC mzd_xor_sse
#define MUL_MC SELECT_V_VL(mzd_mul_v_sse, mzd_mul_vl_sse)

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

#undef SIGN_SBOX
#undef VERIFY_SBOX
#define SIGN_SBOX mpc_sbox_layer_bitsliced_128_sse
#define VERIFY_SBOX mpc_sbox_layer_bitsliced_verify_128_sse

#undef MUL_R_1
#undef MUL_R_10
#undef MUL_Z_1
#undef MUL_Z_10
#define MUL_R_1  mzd_mul_v_sse_3_128
#define MUL_R_10 mzd_mul_v_sse_30_128
#define MUL_Z_1  mzd_mul_v_125_3_popcnt
#define MUL_Z_10 mzd_mul_v_98_30_popcnt

#define SIGN mpc_lowmc_call_128_sse
#define VERIFY mpc_lowmc_call_verify_128_sse
#include "mpc_lowmc.c.i"

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

#undef SIGN_SBOX
#undef VERIFY_SBOX
#define SIGN_SBOX mpc_sbox_layer_bitsliced_256_sse
#define VERIFY_SBOX mpc_sbox_layer_bitsliced_verify_256_sse

#undef MUL_R_1
#undef MUL_R_10
#undef MUL_Z_1
#undef MUL_Z_10
#define MUL_R_1  mzd_mul_v_sse_3_192
#define MUL_R_10 mzd_mul_v_sse_30_192
#define MUL_Z_1  mzd_mul_v_189_3_popcnt
#define MUL_Z_10 mzd_mul_v_162_30_popcnt

#define SIGN mpc_lowmc_call_192_sse
#define VERIFY mpc_lowmc_call_verify_192_sse
#include "mpc_lowmc.c.i"

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
#define MUL_R_1  mzd_mul_v_sse_3_256
#define MUL_R_10 mzd_mul_v_sse_30_256
#define MUL_Z_1  mzd_mul_v_253_3_popcnt
#define MUL_Z_10 mzd_mul_v_226_30_popcnt

#define SIGN mpc_lowmc_call_256_sse
#define VERIFY mpc_lowmc_call_verify_256_sse
#include "mpc_lowmc.c.i"

#endif

#if defined(WITH_AVX2)
#undef XOR_MC
#undef MUL_MC
#define XOR_MC mzd_xor_avx
#define MUL_MC SELECT_V_VL(mzd_mul_v_avx, mzd_mul_vl_avx)

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

#undef SIGN_SBOX
#undef VERIFY_SBOX
#define SIGN_SBOX mpc_sbox_layer_bitsliced_128_sse
#define VERIFY_SBOX mpc_sbox_layer_bitsliced_verify_128_sse

#undef MUL_R_1
#undef MUL_R_10
#undef MUL_Z_1
#undef MUL_Z_10
#define MUL_R_1  mzd_mul_v_sse_3_128
#define MUL_R_10 mzd_mul_v_avx_30_128
#define MUL_Z_1  mzd_mul_v_125_3_popcnt
#define MUL_Z_10 mzd_mul_v_98_30_popcnt

#define SIGN mpc_lowmc_call_128_avx
#define VERIFY mpc_lowmc_call_verify_128_avx
#include "mpc_lowmc.c.i"

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

#undef SIGN_SBOX
#undef VERIFY_SBOX
#define SIGN_SBOX mpc_sbox_layer_bitsliced_256_avx
#define VERIFY_SBOX mpc_sbox_layer_bitsliced_verify_256_avx

#undef MUL_R_1
#undef MUL_R_10
#undef MUL_Z_1
#undef MUL_Z_10
#define MUL_R_1  mzd_mul_v_avx_3_192
#define MUL_R_10 mzd_mul_v_avx_30_192
#define MUL_Z_1  mzd_mul_v_189_3_popcnt
#define MUL_Z_10 mzd_mul_v_162_30_popcnt

#define SIGN mpc_lowmc_call_192_avx
#define VERIFY mpc_lowmc_call_verify_192_avx
#include "mpc_lowmc.c.i"

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
#define MUL_R_1  mzd_mul_v_avx_3_256
#define MUL_R_10 mzd_mul_v_avx_30_256
#define MUL_Z_1  mzd_mul_v_253_3_popcnt
#define MUL_Z_10 mzd_mul_v_226_30_popcnt

#define SIGN mpc_lowmc_call_256_avx
#define VERIFY mpc_lowmc_call_verify_256_avx
#include "mpc_lowmc.c.i"

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

#undef SIGN_SBOX
#undef VERIFY_SBOX
#define SIGN_SBOX mpc_sbox_layer_bitsliced_128_neon
#define VERIFY_SBOX mpc_sbox_layer_bitsliced_verify_128_neon

#undef MUL_R_1
#undef MUL_R_10
#undef MUL_Z_1
#undef MUL_Z_10
#define MUL_R_1  mzd_mul_v_neon_3_128
#define MUL_R_10 mzd_mul_v_neon_30_128
#define MUL_Z_1  mzd_mul_v_125_3_popcnt
#define MUL_Z_10 mzd_mul_v_98_30_popcnt

#define SIGN mpc_lowmc_call_128_neon
#define VERIFY mpc_lowmc_call_verify_128_neon
#include "mpc_lowmc.c.i"

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

#undef SIGN_SBOX
#undef VERIFY_SBOX
#define SIGN_SBOX mpc_sbox_layer_bitsliced_256_neon
#define VERIFY_SBOX mpc_sbox_layer_bitsliced_verify_256_neon

#undef MUL_R_1
#undef MUL_R_10
#undef MUL_Z_1
#undef MUL_Z_10
#define MUL_R_1  mzd_mul_v_neon_3_192
#define MUL_R_10 mzd_mul_v_neon_30_192
#define MUL_Z_1  mzd_mul_v_189_3_popcnt
#define MUL_Z_10 mzd_mul_v_162_30_popcnt

#define SIGN mpc_lowmc_call_192_neon
#define VERIFY mpc_lowmc_call_verify_192_neon
#include "mpc_lowmc.c.i"

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
#define MUL_R_1  mzd_mul_v_neon_3_256
#define MUL_R_10 mzd_mul_v_neon_30_256
#define MUL_Z_1  mzd_mul_v_253_3_popcnt
#define MUL_Z_10 mzd_mul_v_226_30_popcnt

#define SIGN mpc_lowmc_call_256_neon
#define VERIFY mpc_lowmc_call_verify_256_neon
#include "mpc_lowmc.c.i"

#endif
#endif

zkbpp_lowmc_implementation_f get_zkbpp_lowmc_implementation(const lowmc_t* lowmc) {
#if defined(WITH_OPT)
#if defined(WITH_AVX2)
  if (CPU_SUPPORTS_AVX2) {
    if (lowmc->m == 10) {
      switch (lowmc->n) {
        case 128:
          return mpc_lowmc_call_128_avx_10;
        case 192:
          return mpc_lowmc_call_192_avx_10;
        case 256:
          return mpc_lowmc_call_256_avx_10;
      }
    }
    if (lowmc->m == 1) {
      switch (lowmc->n) {
        case 128:
          return mpc_lowmc_call_128_avx_1;
        case 192:
          return mpc_lowmc_call_192_avx_1;
        case 256:
          return mpc_lowmc_call_256_avx_1;
      }
    }
  }
#endif
#if defined(WITH_SSE2)
  if (CPU_SUPPORTS_SSE2) {
    if (lowmc->m == 10) {
      switch (lowmc->n) {
        case 128:
          return mpc_lowmc_call_128_sse_10;
        case 192:
          return mpc_lowmc_call_192_sse_10;
        case 256:
          return mpc_lowmc_call_256_sse_10;
      }
    }
    if (lowmc->m == 1) {
      switch (lowmc->n) {
        case 128:
          return mpc_lowmc_call_128_sse_1;
        case 192:
          return mpc_lowmc_call_192_sse_1;
        case 256:
          return mpc_lowmc_call_256_sse_1;
      }
    }
  }
#endif
#if defined(WITH_NEON)
  if (CPU_SUPPORTS_NEON) {
    if (lowmc->m == 10) {
      switch (lowmc->n) {
      case 128:
        return mpc_lowmc_call_128_neon_10;
      case 192:
        return mpc_lowmc_call_192_neon_10;
      case 256:
        return mpc_lowmc_call_256_neon_10;
      }
    }
    if (lowmc->m == 1) {
      switch (lowmc->n) {
      case 128:
        return mpc_lowmc_call_128_neon_1;
      case 192:
        return mpc_lowmc_call_192_neon_1;
      case 256:
        return mpc_lowmc_call_256_neon_1;
      }
    }
  }
#endif
#endif

  if (lowmc->m == 10)
    return mpc_lowmc_call_10;
  if (lowmc->m == 1)
    return mpc_lowmc_call_1;
  return NULL;
}

zkbpp_lowmc_verify_implementation_f get_zkbpp_lowmc_verify_implementation(const lowmc_t* lowmc) {
#if defined(WITH_OPT)
#if defined(WITH_AVX2)
  if (CPU_SUPPORTS_AVX2) {
    if (lowmc->m == 10) {
      switch (lowmc->n) {
        case 128:
          return mpc_lowmc_call_verify_128_avx_10;
        case 192:
          return mpc_lowmc_call_verify_192_avx_10;
        case 256:
          return mpc_lowmc_call_verify_256_avx_10;
      }
    }
    if (lowmc->m == 1) {
      switch (lowmc->n) {
        case 128:
          return mpc_lowmc_call_verify_128_avx_1;
        case 192:
          return mpc_lowmc_call_verify_192_avx_1;
        case 256:
          return mpc_lowmc_call_verify_256_avx_1;
      }
    }
  }
#endif
#if defined(WITH_SSE2)
  if (CPU_SUPPORTS_SSE2) {
    if (lowmc->m == 10) {
      switch (lowmc->n) {
        case 128:
          return mpc_lowmc_call_verify_128_sse_10;
        case 192:
          return mpc_lowmc_call_verify_192_sse_10;
        case 256:
          return mpc_lowmc_call_verify_256_sse_10;
      }
    }
    if (lowmc->m == 1) {
      switch (lowmc->n) {
        case 128:
          return mpc_lowmc_call_verify_128_sse_1;
        case 192:
          return mpc_lowmc_call_verify_192_sse_1;
        case 256:
          return mpc_lowmc_call_verify_256_sse_1;
      }
    }
  }
#endif
#if defined(WITH_NEON)
  if (CPU_SUPPORTS_NEON) {
    if (lowmc->m == 10) {
      switch (lowmc->n) {
      case 128:
        return mpc_lowmc_call_verify_128_neon_10;
      case 192:
        return mpc_lowmc_call_verify_192_neon_10;
      case 256:
        return mpc_lowmc_call_verify_256_neon_10;
      }
    }
    if (lowmc->m == 1) {
      switch (lowmc->n) {
      case 128:
        return mpc_lowmc_call_verify_128_neon_1;
      case 192:
        return mpc_lowmc_call_verify_192_neon_1;
      case 256:
        return mpc_lowmc_call_verify_256_neon_1;
      }
    }
  }
#endif
#endif

  if (lowmc->m == 10)
    return mpc_lowmc_call_verify_10;
  if (lowmc->m == 1)
    return mpc_lowmc_call_verify_1;
  return NULL;
}
