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

#if defined(WITH_LOWMC_M1)
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

#define SBOX_uint64(sbox, y, x, views, r, n, shares)                                               \
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

#define R_uint64 const uint64_t* r = rvec[i].t

// uint64 based implementation
#define XOR mzd_xor_uint64
#define MUL SELECT_V_VL(mzd_mul_v_uint64, mzd_mul_vl_uint64)
#define ADDMUL SELECT_V_VL(mzd_addmul_v_uint64, mzd_addmul_vl_uint64)
#define XOR_MC mzd_xor_uint64
#define MUL_MC SELECT_V_VL(mzd_mul_v_uint64, mzd_mul_vl_uint64)
#define SHUFFLE mzd_shuffle

#define SIGN_SBOX mpc_sbox_layer_bitsliced
#define VERIFY_SBOX mpc_sbox_layer_bitsliced_verify

#define MUL_R_1  mzd_addmul_v_uint64_3
#define MUL_R_10 mzd_addmul_v_uint64_30
#define MUL_Z_1  mzd_mul_v_parity_uint64_128_3
#define MUL_Z_10 mzd_mul_v_parity_uint64_128_30

#define LOWMC_N LOWMC_L1_N
#define LOWMC_R_10 LOWMC_L1_R
#define LOWMC_R_1 LOWMC_L1_1_R
#if defined(WITH_LOWMC_128_128_20)
#define LOWMC_INSTANCE_10 lowmc_128_128_20
#endif
#if defined(WITH_LOWMC_128_128_182)
#define LOWMC_INSTANCE_1 lowmc_128_128_182
#endif
#define SIGN mpc_lowmc_call_uint64_128
#define VERIFY mpc_lowmc_call_verify_uint64_128
#include "mpc_lowmc.c.i"

#undef MUL_Z_1
#undef MUL_Z_10
#define MUL_Z_1  mzd_mul_v_parity_uint64_192_3
#define MUL_Z_10 mzd_mul_v_parity_uint64_192_30

#undef LOWMC_N
#undef LOWMC_R_10
#undef LOWMC_R_1
#undef LOWMC_INSTANCE_10
#undef LOWMC_INSTANCE_1

#define LOWMC_N LOWMC_L3_N
#define LOWMC_R_10 LOWMC_L3_R
#define LOWMC_R_1 LOWMC_L3_1_R
#if defined(WITH_LOWMC_192_192_30)
#define LOWMC_INSTANCE_10 lowmc_192_192_30
#endif
#if defined(WITH_LOWMC_192_192_284)
#define LOWMC_INSTANCE_1 lowmc_192_192_284
#endif
#define SIGN mpc_lowmc_call_uint64_192
#define VERIFY mpc_lowmc_call_verify_uint64_192
#include "mpc_lowmc.c.i"

#undef MUL_Z_1
#undef MUL_Z_10
#define MUL_Z_1  mzd_mul_v_parity_uint64_256_3
#define MUL_Z_10 mzd_mul_v_parity_uint64_256_30

#undef LOWMC_N
#undef LOWMC_R_10
#undef LOWMC_R_1
#undef LOWMC_INSTANCE_10
#undef LOWMC_INSTANCE_1

#define LOWMC_N LOWMC_L5_N
#define LOWMC_R_10 LOWMC_L5_R
#define LOWMC_R_1 LOWMC_L5_1_R
#if defined(WITH_LOWMC_256_256_38)
#define LOWMC_INSTANCE_10 lowmc_256_256_38
#endif
#if defined(WITH_LOWMC_256_256_363)
#define LOWMC_INSTANCE_1 lowmc_256_256_363
#endif
#define SIGN mpc_lowmc_call_uint64_256
#define VERIFY mpc_lowmc_call_verify_uint64_256
#include "mpc_lowmc.c.i"

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
#define LOWMC_INSTANCE_10 lowmc_128_128_20
#endif
#if defined(WITH_LOWMC_128_128_182)
#define LOWMC_INSTANCE_1 lowmc_128_128_182
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
#define MUL_R_1  mzd_addmul_v_sse_3_128
#define MUL_R_10 mzd_addmul_v_sse_30_128
#define MUL_Z_1  mzd_mul_v_parity_uint64_128_3
#define MUL_Z_10 mzd_mul_v_parity_uint64_128_30

#define SIGN mpc_lowmc_call_sse_128
#define VERIFY mpc_lowmc_call_verify_sse_128
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
#define LOWMC_INSTANCE_10 lowmc_192_192_30
#endif
#if defined(WITH_LOWMC_192_192_284)
#define LOWMC_INSTANCE_1 lowmc_192_192_284
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
#define MUL_R_1  mzd_addmul_v_sse_3_192
#define MUL_R_10 mzd_addmul_v_sse_30_192
#define MUL_Z_1  mzd_mul_v_parity_uint64_192_3
#define MUL_Z_10 mzd_mul_v_parity_uint64_192_30

#define SIGN mpc_lowmc_call_sse_192
#define VERIFY mpc_lowmc_call_verify_sse_192
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
#define LOWMC_INSTANCE_10 lowmc_256_256_38
#endif
#if defined(WITH_LOWMC_256_256_363)
#define LOWMC_INSTANCE_1 lowmc_256_256_363
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

#define SIGN mpc_lowmc_call_sse_256
#define VERIFY mpc_lowmc_call_verify_sse_256
#include "mpc_lowmc.c.i"

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
#define LOWMC_INSTANCE_10 lowmc_128_128_20
#endif
#if defined(WITH_LOWMC_128_128_182)
#define LOWMC_INSTANCE_1 lowmc_128_128_182
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
#define MUL_R_1  mzd_addmul_v_sse_3_128
#define MUL_R_10 mzd_addmul_v_sse_30_128
#define MUL_Z_1  mzd_mul_v_parity_popcnt_128_3
#define MUL_Z_10 mzd_mul_v_parity_popcnt_128_30

#define SIGN mpc_lowmc_call_sse_popcnt_128
#define VERIFY mpc_lowmc_call_verify_sse_popcnt_128
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
#define LOWMC_INSTANCE_10 lowmc_192_192_30
#endif
#if defined(WITH_LOWMC_192_192_284)
#define LOWMC_INSTANCE_1 lowmc_192_192_284
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
#define MUL_R_1  mzd_addmul_v_sse_3_192
#define MUL_R_10 mzd_addmul_v_sse_30_192
#define MUL_Z_1  mzd_mul_v_parity_popcnt_192_3
#define MUL_Z_10 mzd_mul_v_parity_popcnt_192_30

#define SIGN mpc_lowmc_call_sse_popcnt_192
#define VERIFY mpc_lowmc_call_verify_sse_popcnt_192
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
#define LOWMC_INSTANCE_10 lowmc_256_256_38
#endif
#if defined(WITH_LOWMC_256_256_363)
#define LOWMC_INSTANCE_1 lowmc_256_256_363
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

#define SIGN mpc_lowmc_call_sse_popcnt_256
#define VERIFY mpc_lowmc_call_verify_sse_popcnt_256
#include "mpc_lowmc.c.i"

#undef FN_ATTR
#endif

#if defined(WITH_AVX2)
#undef XOR_MC
#undef MUL_MC
#undef SHUFFLE
#define XOR_MC mzd_xor_avx
#define MUL_MC SELECT_V_VL(mzd_mul_v_avx, mzd_mul_vl_avx)
#define FN_ATTR ATTR_TARGET("avx2,bmi2")
#define SHUFFLE mzd_shuffle_pext

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
#define LOWMC_INSTANCE_10 lowmc_128_128_20
#endif
#if defined(WITH_LOWMC_128_128_182)
#define LOWMC_INSTANCE_1 lowmc_128_128_182
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
#define MUL_R_1  mzd_addmul_v_sse_3_128
#define MUL_R_10 mzd_addmul_v_avx_30_128
#define MUL_Z_1  mzd_mul_v_parity_uint64_128_3
#define MUL_Z_10 mzd_mul_v_parity_uint64_128_30

#define SIGN mpc_lowmc_call_avx_128
#define VERIFY mpc_lowmc_call_verify_avx_128
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
#define LOWMC_INSTANCE_10 lowmc_192_192_30
#endif
#if defined(WITH_LOWMC_192_192_284)
#define LOWMC_INSTANCE_1 lowmc_192_192_284
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
#define MUL_R_1  mzd_addmul_v_avx_3_192
#define MUL_R_10 mzd_addmul_v_avx_30_192
#define MUL_Z_1  mzd_mul_v_parity_uint64_192_3
#define MUL_Z_10 mzd_mul_v_parity_uint64_192_30

#define SIGN mpc_lowmc_call_avx_192
#define VERIFY mpc_lowmc_call_verify_avx_192
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
#define LOWMC_INSTANCE_10 lowmc_256_256_38
#endif
#if defined(WITH_LOWMC_256_256_363)
#define LOWMC_INSTANCE_1 lowmc_256_256_363
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

#define SIGN mpc_lowmc_call_avx_256
#define VERIFY mpc_lowmc_call_verify_avx_256
#include "mpc_lowmc.c.i"

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
#define LOWMC_INSTANCE_10 lowmc_128_128_20
#endif
#if defined(WITH_LOWMC_128_128_182)
#define LOWMC_INSTANCE_1 lowmc_128_128_182
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
#define MUL_R_1  mzd_addmul_v_sse_3_128
#define MUL_R_10 mzd_addmul_v_avx_30_128
#define MUL_Z_1  mzd_mul_v_parity_popcnt_128_3
#define MUL_Z_10 mzd_mul_v_parity_popcnt_128_30

#define SIGN mpc_lowmc_call_avx_popcnt_128
#define VERIFY mpc_lowmc_call_verify_avx_popcnt_128
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
#define LOWMC_INSTANCE_10 lowmc_192_192_30
#endif
#if defined(WITH_LOWMC_192_192_284)
#define LOWMC_INSTANCE_1 lowmc_192_192_284
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
#define MUL_R_1  mzd_addmul_v_avx_3_192
#define MUL_R_10 mzd_addmul_v_avx_30_192
#define MUL_Z_1  mzd_mul_v_parity_popcnt_192_3
#define MUL_Z_10 mzd_mul_v_parity_popcnt_192_30

#define SIGN mpc_lowmc_call_avx_popcnt_192
#define VERIFY mpc_lowmc_call_verify_avx_popcnt_192
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
#define LOWMC_INSTANCE_10 lowmc_256_256_38
#endif
#if defined(WITH_LOWMC_256_256_363)
#define LOWMC_INSTANCE_1 lowmc_256_256_363
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

#define SIGN mpc_lowmc_call_avx_popcnt_256
#define VERIFY mpc_lowmc_call_verify_avx_popcnt_256
#include "mpc_lowmc.c.i"

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
#define LOWMC_INSTANCE_10 lowmc_128_128_20
#endif
#if defined(WITH_LOWMC_128_128_182)
#define LOWMC_INSTANCE_1 lowmc_128_128_182
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
#define MUL_R_1  mzd_addmul_v_neon_3_128
#define MUL_R_10 mzd_addmul_v_neon_30_128
#define MUL_Z_1  mzd_mul_v_parity_uint64_128_3
#define MUL_Z_10 mzd_mul_v_parity_uint64_128_30

#define SIGN mpc_lowmc_call_neon_128
#define VERIFY mpc_lowmc_call_verify_neon_128
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
#define LOWMC_INSTANCE_10 lowmc_192_192_30
#endif
#if defined(WITH_LOWMC_192_192_284)
#define LOWMC_INSTANCE_1 lowmc_192_192_284
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
#define MUL_R_1  mzd_addmul_v_neon_3_192
#define MUL_R_10 mzd_addmul_v_neon_30_192
#define MUL_Z_1  mzd_mul_v_parity_uint64_192_3
#define MUL_Z_10 mzd_mul_v_parity_uint64_192_30

#define SIGN mpc_lowmc_call_neon_192
#define VERIFY mpc_lowmc_call_verify_neon_192
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
#define LOWMC_INSTANCE_10 lowmc_256_256_38
#endif
#if defined(WITH_LOWMC_256_256_363)
#define LOWMC_INSTANCE_1 lowmc_256_256_363
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

#define SIGN mpc_lowmc_call_neon_256
#define VERIFY mpc_lowmc_call_verify_neon_256
#include "mpc_lowmc.c.i"

#endif
#endif

zkbpp_lowmc_implementation_f get_zkbpp_lowmc_implementation(const lowmc_t* lowmc) {
#if defined(WITH_OPT)
#if defined(WITH_AVX2)
  if (CPU_SUPPORTS_AVX2) {
    if (lowmc->m == 10) {
#if defined(WITH_POPCNT)
      if (CPU_SUPPORTS_POPCNT) {
        switch (lowmc->n) {
#if defined(WITH_LOWMC_128_128_20)
        case 128:
          return mpc_lowmc_call_avx_popcnt_128_10;
#endif
#if defined(WITH_LOWMC_192_192_30)
        case 192:
          return mpc_lowmc_call_avx_popcnt_192_10;
#endif
#if defined(WITH_LOWMC_256_256_38)
        case 256:
          return mpc_lowmc_call_avx_popcnt_256_10;
#endif
        }
      }
#endif
      switch (lowmc->n) {
#if defined(WITH_LOWMC_128_128_20)
      case 128:
        return mpc_lowmc_call_avx_128_10;
#endif
#if defined(WITH_LOWMC_192_192_30)
      case 192:
        return mpc_lowmc_call_avx_192_10;
#endif
#if defined(WITH_LOWMC_256_256_38)
      case 256:
        return mpc_lowmc_call_avx_256_10;
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
          return mpc_lowmc_call_avx_popcnt_128_1;
#endif
#if defined(WITH_LOWMC_192_192_284)
        case 192:
          return mpc_lowmc_call_avx_popcnt_192_1;
#endif
#if defined(WITH_LOWMC_256_256_363)
        case 256:
          return mpc_lowmc_call_avx_popcnt_256_1;
#endif
        }
      }
#endif
      switch (lowmc->n) {
#if defined(WITH_LOWMC_128_128_182)
      case 128:
        return mpc_lowmc_call_avx_128_1;
#endif
#if defined(WITH_LOWMC_192_192_284)
      case 192:
        return mpc_lowmc_call_avx_192_1;
#endif
#if defined(WITH_LOWMC_256_256_363)
      case 256:
        return mpc_lowmc_call_avx_256_1;
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
          return mpc_lowmc_call_sse_popcnt_128_10;
#endif
#if defined(WITH_LOWMC_192_192_30)
        case 192:
          return mpc_lowmc_call_sse_popcnt_192_10;
#endif
#if defined(WITH_LOWMC_256_256_38)
        case 256:
          return mpc_lowmc_call_sse_popcnt_256_10;
#endif
        }
      }
#endif
      switch (lowmc->n) {
#if defined(WITH_LOWMC_128_128_20)
      case 128:
        return mpc_lowmc_call_sse_128_10;
#endif
#if defined(WITH_LOWMC_192_192_30)
      case 192:
        return mpc_lowmc_call_sse_192_10;
#endif
#if defined(WITH_LOWMC_256_256_38)
      case 256:
        return mpc_lowmc_call_sse_256_10;
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
          return mpc_lowmc_call_sse_popcnt_128_1;
#endif
#if defined(WITH_LOWMC_192_192_284)
        case 192:
          return mpc_lowmc_call_sse_popcnt_192_1;
#endif
#if defined(WITH_LOWMC_256_256_363)
        case 256:
          return mpc_lowmc_call_sse_popcnt_256_1;
#endif
        }
      }
#endif
      switch (lowmc->n) {
#if defined(WITH_LOWMC_128_128_182)
      case 128:
        return mpc_lowmc_call_sse_128_1;
#endif
#if defined(WITH_LOWMC_192_192_284)
      case 192:
        return mpc_lowmc_call_sse_192_1;
#endif
#if defined(WITH_LOWMC_256_256_363)
      case 256:
        return mpc_lowmc_call_sse_256_1;
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
        return mpc_lowmc_call_neon_128_10;
#endif
#if defined(WITH_LOWMC_192_192_30)
      case 192:
        return mpc_lowmc_call_neon_192_10;
#endif
#if defined(WITH_LOWMC_256_256_38)
      case 256:
        return mpc_lowmc_call_neon_256_10;
#endif
      }
    }
#if defined(WITH_LOWMC_M1)
    if (lowmc->m == 1) {
      switch (lowmc->n) {
#if defined(WITH_LOWMC_128_128_182)
      case 128:
        return mpc_lowmc_call_neon_128_1;
#endif
#if defined(WITH_LOWMC_192_192_284)
      case 192:
        return mpc_lowmc_call_neon_192_1;
#endif
#if defined(WITH_LOWMC_256_256_363)
      case 256:
        return mpc_lowmc_call_neon_256_1;
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
      return mpc_lowmc_call_uint64_128_10;
#endif
#if defined(WITH_LOWMC_192_192_30)
    case 192:
      return mpc_lowmc_call_uint64_192_10;
#endif
#if defined(WITH_LOWMC_256_256_38)
    case 256:
      return mpc_lowmc_call_uint64_256_10;
#endif
    }
  }

#if defined(WITH_LOWMC_M1)
  if (lowmc->m == 1) {
    switch (lowmc->n) {
#if defined(WITH_LOWMC_128_128_182)
    case 128:
      return mpc_lowmc_call_uint64_128_1;
#endif
#if defined(WITH_LOWMC_192_192_284)
    case 192:
      return mpc_lowmc_call_uint64_192_1;
#endif
#if defined(WITH_LOWMC_256_256_363)
    case 256:
      return mpc_lowmc_call_uint64_256_1;
#endif
    }
  }
#endif

  return NULL;
}

zkbpp_lowmc_verify_implementation_f get_zkbpp_lowmc_verify_implementation(const lowmc_t* lowmc) {
#if defined(WITH_OPT)
#if defined(WITH_AVX2)
  if (CPU_SUPPORTS_AVX2) {
    if (lowmc->m == 10) {
#if defined(WITH_POPCNT)
      if (CPU_SUPPORTS_POPCNT) {
        switch (lowmc->n) {
#if defined(WITH_LOWMC_128_128_20)
        case 128:
          return mpc_lowmc_call_verify_avx_popcnt_128_10;
#endif
#if defined(WITH_LOWMC_192_192_30)
        case 192:
          return mpc_lowmc_call_verify_avx_popcnt_192_10;
#endif
#if defined(WITH_LOWMC_256_256_38)
        case 256:
          return mpc_lowmc_call_verify_avx_popcnt_256_10;
#endif
        }
      }
#endif
      switch (lowmc->n) {
#if defined(WITH_LOWMC_128_128_20)
      case 128:
        return mpc_lowmc_call_verify_avx_128_10;
#endif
#if defined(WITH_LOWMC_192_192_30)
      case 192:
        return mpc_lowmc_call_verify_avx_192_10;
#endif
#if defined(WITH_LOWMC_256_256_38)
      case 256:
        return mpc_lowmc_call_verify_avx_256_10;
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
          return mpc_lowmc_call_verify_avx_popcnt_128_1;
#endif
#if defined(WITH_LOWMC_192_192_284)
        case 192:
          return mpc_lowmc_call_verify_avx_popcnt_192_1;
#endif
#if defined(WITH_LOWMC_256_256_363)
        case 256:
          return mpc_lowmc_call_verify_avx_popcnt_256_1;
#endif
        }
      }
#endif
      switch (lowmc->n) {
#if defined(WITH_LOWMC_128_128_182)
      case 128:
        return mpc_lowmc_call_verify_avx_128_1;
#endif
#if defined(WITH_LOWMC_192_192_284)
      case 192:
        return mpc_lowmc_call_verify_avx_192_1;
#endif
#if defined(WITH_LOWMC_256_256_363)
      case 256:
        return mpc_lowmc_call_verify_avx_256_1;
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
          return mpc_lowmc_call_verify_sse_popcnt_128_10;
#endif
#if defined(WITH_LOWMC_192_192_30)
        case 192:
          return mpc_lowmc_call_verify_sse_popcnt_192_10;
#endif
#if defined(WITH_LOWMC_256_256_38)
        case 256:
          return mpc_lowmc_call_verify_sse_popcnt_256_10;
#endif
        }
      }
#endif
      switch (lowmc->n) {
#if defined(WITH_LOWMC_128_128_20)
      case 128:
        return mpc_lowmc_call_verify_sse_128_10;
#endif
#if defined(WITH_LOWMC_192_192_30)
      case 192:
        return mpc_lowmc_call_verify_sse_192_10;
#endif
#if defined(WITH_LOWMC_256_256_38)
      case 256:
        return mpc_lowmc_call_verify_sse_256_10;
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
          return mpc_lowmc_call_verify_sse_popcnt_128_1;
#endif
#if defined(WITH_LOWMC_192_192_284)
        case 192:
          return mpc_lowmc_call_verify_sse_popcnt_192_1;
#endif
#if defined(WITH_LOWMC_256_256_363)
        case 256:
          return mpc_lowmc_call_verify_sse_popcnt_256_1;
#endif
        }
      }
#endif
      switch (lowmc->n) {
#if defined(WITH_LOWMC_128_128_182)
      case 128:
        return mpc_lowmc_call_verify_sse_128_1;
#endif
#if defined(WITH_LOWMC_192_192_284)
      case 192:
        return mpc_lowmc_call_verify_sse_192_1;
#endif
#if defined(WITH_LOWMC_256_256_363)
      case 256:
        return mpc_lowmc_call_verify_sse_256_1;
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
        return mpc_lowmc_call_verify_neon_128_10;
#endif
#if defined(WITH_LOWMC_192_192_30)
      case 192:
        return mpc_lowmc_call_verify_neon_192_10;
#endif
#if defined(WITH_LOWMC_256_256_38)
      case 256:
        return mpc_lowmc_call_verify_neon_256_10;
#endif
      }
    }
#if defined(WITH_LOWMC_M1)
    if (lowmc->m == 1) {
      switch (lowmc->n) {
#if defined(WITH_LOWMC_128_128_182)
      case 128:
        return mpc_lowmc_call_verify_neon_128_1;
#endif
#if defined(WITH_LOWMC_192_192_284)
      case 192:
        return mpc_lowmc_call_verify_neon_192_1;
#endif
#if defined(WITH_LOWMC_256_256_363)
      case 256:
        return mpc_lowmc_call_verify_neon_256_1;
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
      return mpc_lowmc_call_verify_uint64_128_10;
#endif
#if defined(WITH_LOWMC_192_192_30)
    case 192:
      return mpc_lowmc_call_verify_uint64_192_10;
#endif
#if defined(WITH_LOWMC_256_256_38)
    case 256:
      return mpc_lowmc_call_verify_uint64_256_10;
#endif
    }
  }

#if defined(WITH_LOWMC_M1)
  if (lowmc->m == 1) {
    switch (lowmc->n) {
#if defined(WITH_LOWMC_128_128_182)
    case 128:
      return mpc_lowmc_call_verify_uint64_128_1;
#endif
#if defined(WITH_LOWMC_192_192_284)
    case 192:
      return mpc_lowmc_call_verify_uint64_192_1;
#endif
#if defined(WITH_LOWMC_256_256_363)
    case 256:
      return mpc_lowmc_call_verify_uint64_256_1;
#endif
    }
  }
#endif

  return NULL;
}
