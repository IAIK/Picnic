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

#if defined(WITH_CUSTOM_INSTANCES)
static void sbox_layer_bitsliced(mzd_local_t* in, mask_t const* mask) {
  mzd_local_t* buffer[6] = {NULL};
  mzd_local_init_multiple_ex(buffer, 6, 1, in->ncols, false);

  mzd_local_t* x0m = buffer[0];
  mzd_local_t* x1m = buffer[1];
  mzd_local_t* x2m = buffer[2];
  mzd_local_t* t0  = buffer[3];
  mzd_local_t* t1  = buffer[4];
  mzd_local_t* t2  = buffer[5];

  // a
  mzd_and(x0m, mask->x0, in);
  // b
  mzd_and(x1m, mask->x1, in);
  // c
  mzd_and(x2m, mask->x2, in);

  mzd_shift_left(x0m, x0m, 2);
  mzd_shift_left(x1m, x1m, 1);

  // b & c
  mzd_and(t0, x1m, x2m);
  // c & a
  mzd_and(t1, x0m, x2m);
  // a & b
  mzd_and(t2, x0m, x1m);

  // (b & c) ^ a
  mzd_xor(t0, t0, x0m);

  // (c & a) ^ a ^ b
  mzd_xor(t1, t1, x0m);
  mzd_xor(t1, t1, x1m);

  // (a & b) ^ a ^ b ^c
  mzd_xor(t2, t2, x0m);
  mzd_xor(t2, t2, x1m);
  mzd_xor(t2, t2, x2m);

  mzd_shift_right(t0, t0, 2);
  mzd_shift_right(t1, t1, 1);

  mzd_and(in, in, mask->mask);
  mzd_xor(in, in, t2);
  mzd_xor(in, in, t0);
  mzd_xor(in, in, t1);

  mzd_local_free_multiple(buffer);
}

#if defined(WITH_OPT)
#if defined(WITH_SSE2)
ATTR_TARGET("sse") static void sbox_layer_sse(mzd_local_t* in, mask_t const* mask) {
  __m128i* ip       = (__m128i*)ASSUME_ALIGNED(CONST_FIRST_ROW(in), alignof(__m128i));
  __m128i const min = *ip;

  __m128i const* x0p = (__m128i const*)ASSUME_ALIGNED(CONST_FIRST_ROW(mask->x0), alignof(__m128i));
  __m128i const* x1p = (__m128i const*)ASSUME_ALIGNED(CONST_FIRST_ROW(mask->x1), alignof(__m128i));
  __m128i const* x2p = (__m128i const*)ASSUME_ALIGNED(CONST_FIRST_ROW(mask->x2), alignof(__m128i));

  __m128i x0m = _mm_and_si128(min, *x0p);
  __m128i x1m = _mm_and_si128(min, *x1p);
  __m128i x2m = _mm_and_si128(min, *x2p);

  x0m = mm128_shift_left(x0m, 2);
  x1m = mm128_shift_left(x1m, 1);

  __m128i t0 = _mm_and_si128(x1m, x2m);
  __m128i t1 = _mm_and_si128(x0m, x2m);
  __m128i t2 = _mm_and_si128(x0m, x1m);

  t0 = _mm_xor_si128(t0, x0m);

  x0m = _mm_xor_si128(x0m, x1m);
  t1  = _mm_xor_si128(t1, x0m);

  t2 = _mm_xor_si128(t2, x0m);
  t2 = _mm_xor_si128(t2, x2m);

  t0 = mm128_shift_right(t0, 2);
  t1 = mm128_shift_right(t1, 1);

  __m128i const* xmp =
      (__m128i const*)ASSUME_ALIGNED(CONST_FIRST_ROW(mask->mask), alignof(__m128i));

  __m128i mout = _mm_and_si128(min, *xmp);

  mout = _mm_xor_si128(mout, t2);
  mout = _mm_xor_si128(mout, t1);
  *ip  = _mm_xor_si128(mout, t0);
}
#endif

#if defined(WITH_AVX2)
ATTR_TARGET("avx2") static void sbox_layer_avx(mzd_local_t* in, mask_t const* mask) {
  __m256i* ip       = (__m256i*)ASSUME_ALIGNED(CONST_FIRST_ROW(in), alignof(__m256i));
  __m256i const min = *ip;

  __m256i const* x0p = (__m256i const*)ASSUME_ALIGNED(CONST_FIRST_ROW(mask->x0), alignof(__m256i));
  __m256i const* x1p = (__m256i const*)ASSUME_ALIGNED(CONST_FIRST_ROW(mask->x1), alignof(__m256i));
  __m256i const* x2p = (__m256i const*)ASSUME_ALIGNED(CONST_FIRST_ROW(mask->x2), alignof(__m256i));

  __m256i x0m = _mm256_and_si256(min, *x0p);
  __m256i x1m = _mm256_and_si256(min, *x1p);
  __m256i x2m = _mm256_and_si256(min, *x2p);

  x0m = mm256_shift_left(x0m, 2);
  x1m = mm256_shift_left(x1m, 1);

  __m256i t0 = _mm256_and_si256(x1m, x2m);
  __m256i t1 = _mm256_and_si256(x0m, x2m);
  __m256i t2 = _mm256_and_si256(x0m, x1m);

  t0 = _mm256_xor_si256(t0, x0m);

  x0m = _mm256_xor_si256(x0m, x1m);
  t1  = _mm256_xor_si256(t1, x0m);

  t2 = _mm256_xor_si256(t2, x0m);
  t2 = _mm256_xor_si256(t2, x2m);

  t0 = mm256_shift_right(t0, 2);
  t1 = mm256_shift_right(t1, 1);

  __m256i const* xmp =
      (__m256i const*)ASSUME_ALIGNED(CONST_FIRST_ROW(mask->mask), alignof(__m256i));

  __m256i mout = _mm256_and_si256(min, *xmp);

  mout = _mm256_xor_si256(mout, t2);
  mout = _mm256_xor_si256(mout, t1);
  *ip  = _mm256_xor_si256(mout, t0);
}
#endif

#if defined(WITH_NEON)
static void sbox_layer_neon(mzd_local_t* in, mask_t const* mask) {
  uint32x4_t* ip       = (uint32x4_t*)ASSUME_ALIGNED(CONST_FIRST_ROW(in), alignof(uint32x4_t));
  uint32x4_t const min = *ip;

  uint32x4_t const* x0p =
      (uint32x4_t const*)ASSUME_ALIGNED(CONST_FIRST_ROW(mask->x0), alignof(uint32x4_t));
  uint32x4_t const* x1p =
      (uint32x4_t const*)ASSUME_ALIGNED(CONST_FIRST_ROW(mask->x1), alignof(uint32x4_t));
  uint32x4_t const* x2p =
      (uint32x4_t const*)ASSUME_ALIGNED(CONST_FIRST_ROW(mask->x2), alignof(uint32x4_t));

  uint32x4_t x0m = vandq_u32(min, *x0p);
  uint32x4_t x1m = vandq_u32(min, *x1p);
  uint32x4_t x2m = vandq_u32(min, *x2p);

  x0m = mm128_shift_left(x0m, 2);
  x1m = mm128_shift_left(x1m, 1);

  uint32x4_t t0 = vandq_u32(x1m, x2m);
  uint32x4_t t1 = vandq_u32(x0m, x2m);
  uint32x4_t t2 = vandq_u32(x0m, x1m);

  t0 = veorq_u32(t0, x0m);

  x0m = veorq_u32(x0m, x1m);
  t1  = veorq_u32(t1, x0m);

  t2 = veorq_u32(t2, x0m);
  t2 = veorq_u32(t2, x2m);

  t0 = mm128_shift_right(t0, 2);
  t1 = mm128_shift_right(t1, 1);

  uint32x4_t const* xmp =
      (uint32x4_t const*)ASSUME_ALIGNED(CONST_FIRST_ROW(mask->mask), alignof(uint32x4_t));

  uint32x4_t mout = vandq_u32(min, *xmp);

  mout = veorq_u32(mout, t2);
  mout = veorq_u32(mout, t1);
  *ip  = veorq_u32(mout, t0);
}
#endif
#endif
#endif

// uint64 based implementation
#define XOR mzd_xor_uint64
#define MUL SELECT_V_VL(mzd_mul_v_uint64, mzd_mul_vl_uint64)
#define ADDMUL SELECT_V_VL(mzd_addmul_v_uint64, mzd_addmul_vl_uint64)
#define XOR_MC mzd_xor_uint64
#define MUL_MC SELECT_V_VL(mzd_mul_v_uint64, mzd_mul_vl_uint64)

#undef MUL_Z_1
#undef MUL_Z_10
#undef MUL_A_1
#undef MUL_A_10
#define MUL_Z_1  mzd_mul_v_uint64_3
#define MUL_Z_10 mzd_mul_v_uint64_30
#define MUL_A_1  mzd_mul_v_3_popcnt
#define MUL_A_10 mzd_mul_v_3_popcnt

#define LOWMC_N lowmc->n
#define LOWMC_R lowmc->r
#define LOWMC_M lowmc->m

#define SBOX_IMPL sbox_layer_bitsliced
#define LOWMC lowmc_uint64
#include "lowmc.c.i"

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
#if defined(WITH_LOWMC_128_128_20)
#include "lowmc_128_128_20.h"
#endif
#if defined(WITH_LOWMC_192_192_30)
#include "lowmc_192_192_30.h"
#endif
#if defined(WITH_LOWMC_256_256_38)
#include "lowmc_256_256_38.h"
#endif

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

#undef LOWMC_INSTANCE
#undef LOWMC_N
#undef LOWMC_R
#if defined(WITH_LOWMC_128_128_20)
#define LOWMC_INSTANCE (&lowmc_128_128_20)
#endif
#define LOWMC_N LOWMC_L1_N
#define LOWMC_R LOWMC_L1_R

#undef MUL_Z_1
#undef MUL_Z_10
#undef MUL_A_1
#undef MUL_A_10
#define MUL_Z_1  mzd_mul_v_sse_3_128
#define MUL_Z_10 mzd_mul_v_sse_30_128
#define MUL_A_1  mzd_mul_v_125_3_popcnt
#define MUL_A_10 mzd_mul_v_98_30_popcnt

#undef SBOX_IMPL
#undef LOWMC
#define SBOX_IMPL sbox_layer_sse
#define LOWMC lowmc_sse_128
#include "lowmc.c.i"

// L3 using SSE2
#undef XOR
#undef MUL
#undef ADDMUL
#define XOR mzd_xor_sse_256
#define MUL SELECT_V_VL(mzd_mul_v_sse_192, mzd_mul_vl_sse_192)
#define ADDMUL SELECT_V_VL(mzd_addmul_v_sse_192, mzd_addmul_vl_sse_192)

#undef LOWMC_INSTANCE
#undef LOWMC_N
#undef LOWMC_R
#if defined(WITH_LOWMC_192_192_30)
#define LOWMC_INSTANCE (&lowmc_192_192_30)
#endif
#define LOWMC_N LOWMC_L3_N
#define LOWMC_R LOWMC_L3_R

#undef MUL_Z_1
#undef MUL_Z_10
#undef MUL_A_1
#undef MUL_A_10
#define MUL_Z_1  mzd_mul_v_sse_3_192
#define MUL_Z_10 mzd_mul_v_sse_30_192
#define MUL_A_1  mzd_mul_v_189_3_popcnt
#define MUL_A_10 mzd_mul_v_162_30_popcnt

#undef SBOX_IMPL
#undef LOWMC
#define SBOX_IMPL sbox_layer_bitsliced
#define LOWMC lowmc_sse_192
#include "lowmc.c.i"

// L5 using SSE2
#undef MUL
#undef ADDMUL
#define MUL SELECT_V_VL(mzd_mul_v_sse_256, mzd_mul_vl_sse_256)
#define ADDMUL SELECT_V_VL(mzd_addmul_v_sse_256, mzd_addmul_vl_sse_256)

#undef LOWMC_INSTANCE
#undef LOWMC_N
#undef LOWMC_R
#if defined(WITH_LOWMC_256_256_38)
#define LOWMC_INSTANCE (&lowmc_256_256_38)
#endif
#define LOWMC_N LOWMC_L5_N
#define LOWMC_R LOWMC_L5_R

#undef MUL_Z_1
#undef MUL_Z_10
#undef MUL_A_1
#undef MUL_A_10
#define MUL_Z_1  mzd_mul_v_sse_3_256
#define MUL_Z_10 mzd_mul_v_sse_30_256
#define MUL_A_1  mzd_mul_v_253_3_popcnt
#define MUL_A_10 mzd_mul_v_226_30_popcnt

#undef LOWMC
#define LOWMC lowmc_sse_256
#include "lowmc.c.i"

#if defined(WITH_CUSTOM_INSTANCES)
#undef XOR
#undef MUL
#undef ADDMUL
#define XOR mzd_xor_sse
#define MUL SELECT_V_VL(mzd_mul_v_sse, mzd_mul_vl_sse)
#define ADDMUL SELECT_V_VL(mzd_addmul_v_sse, mzd_addmul_vl_sse)

#undef LOWMC_INSTANCE
#undef LOWMC_N
#undef LOWMC_R
#define LOWMC_N lowmc->n
#define LOWMC_R lowmc->r

// generic using SSE2
#undef LOWMC
#define LOWMC lowmc_sse
#include "lowmc.c.i"
#endif
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

#undef LOWMC_INSTANCE
#undef LOWMC_N
#undef LOWMC_R
#if defined(WITH_LOWMC_128_128_20)
#define LOWMC_INSTANCE (&lowmc_128_128_20)
#endif
#define LOWMC_N LOWMC_L1_N
#define LOWMC_R LOWMC_L1_R

#undef MUL_Z_1
#undef MUL_Z_10
#undef MUL_A_1
#undef MUL_A_10
#define MUL_Z_1  mzd_mul_v_sse_3_128
#define MUL_Z_10 mzd_mul_v_avx_30_128
#define MUL_A_1  mzd_mul_v_125_3_popcnt
#define MUL_A_10 mzd_mul_v_98_30_popcnt

#undef SBOX_IMPL
#undef LOWMC
#define SBOX_IMPL sbox_layer_sse
#define LOWMC lowmc_avx_128
#include "lowmc.c.i"

// L3 using AVX2
#undef XOR
#undef MUL
#undef ADDMUL
#define XOR mzd_xor_avx_256
#define MUL SELECT_V_VL(mzd_mul_v_avx_192, mzd_mul_vl_avx_192)
#define ADDMUL SELECT_V_VL(mzd_addmul_v_avx_192, mzd_addmul_vl_avx_192)

#undef LOWMC_INSTANCE
#undef LOWMC_N
#undef LOWMC_R
#if defined(WITH_LOWMC_192_192_30)
#define LOWMC_INSTANCE (&lowmc_192_192_30)
#endif
#define LOWMC_N LOWMC_L3_N
#define LOWMC_R LOWMC_L3_R

#undef MUL_Z_1
#undef MUL_Z_10
#undef MUL_A_1
#undef MUL_A_10
#define MUL_Z_1  mzd_mul_v_avx_3_192
#define MUL_Z_10 mzd_mul_v_avx_30_192
#define MUL_A_1  mzd_mul_v_189_3_popcnt
#define MUL_A_10 mzd_mul_v_162_30_popcnt

#undef SBOX_IMPL
#undef LOWMC
#define SBOX_IMPL sbox_layer_avx
#define LOWMC lowmc_avx_192
#include "lowmc.c.i"

// L5 using AVX2
#undef MUL
#undef ADDMUL
#define MUL SELECT_V_VL(mzd_mul_v_avx_256, mzd_mul_vl_avx_256)
#define ADDMUL SELECT_V_VL(mzd_addmul_v_avx_256, mzd_addmul_vl_avx_256)

#undef LOWMC_INSTANCE
#undef LOWMC_N
#undef LOWMC_R
#if defined(WITH_LOWMC_256_256_38)
#define LOWMC_INSTANCE (&lowmc_256_256_38)
#endif
#define LOWMC_N LOWMC_L5_N
#define LOWMC_R LOWMC_L5_R

#undef MUL_Z_1
#undef MUL_Z_10
#undef MUL_A_1
#undef MUL_A_10
#define MUL_Z_1  mzd_mul_v_avx_3_256
#define MUL_Z_10 mzd_mul_v_avx_30_256
#define MUL_A_1  mzd_mul_v_253_3_popcnt
#define MUL_A_10 mzd_mul_v_226_30_popcnt

#undef LOWMC
#define LOWMC lowmc_avx_256
#include "lowmc.c.i"

#if defined(WITH_CUSTOM_INSTANCES)
#undef XOR
#undef MUL
#undef ADDMUL
#define XOR mzd_xor_avx
#define MUL SELECT_V_VL(mzd_mul_v_avx, mzd_mul_vl_avx)
#define ADDMUL SELECT_V_VL(mzd_addmul_v_avx, mzd_addmul_vl_avx)

#undef LOWMC_INSTANCE
#undef LOWMC_N
#undef LOWMC_R
#define LOWMC_N lowmc->n
#define LOWMC_R lowmc->r

// generic using AVX2
#undef SBOX_IMPL
#undef LOWMC
#define SBOX_IMPL sbox_layer_bitsliced
#define LOWMC lowmc_avx
#include "lowmc.c.i"
#endif
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

#undef LOWMC_N
#undef LOWMC_R
#define LOWMC_N LOWMC_L1_N
#define LOWMC_R LOWMC_L1_R

#undef SBOX_IMPL
#undef LOWMC
#define SBOX_IMPL sbox_layer_neon
#define LOWMC lowmc_neon_128
#include "lowmc.c.i"

// L3 using NEON
#undef XOR
#undef MUL
#undef ADDMUL
#define XOR mzd_xor_neon_256
#define MUL SELECT_V_VL(mzd_mul_v_neon_192, mzd_mul_vl_neon_192)
#define ADDMUL SELECT_V_VL(mzd_addmul_v_neon_192, mzd_addmul_vl_neon_192)

#undef LOWMC_N
#undef LOWMC_R
#define LOWMC_N LOWMC_L3_N
#define LOWMC_R LOWMC_L3_R

#undef SBOX_IMPL
#undef LOWMC
#define SBOX_IMPL sbox_layer_bitsliced
#define LOWMC lowmc_neon_192
#include "lowmc.c.i"

// L5 using NEON
#undef MUL
#undef ADDMUL
#define MUL SELECT_V_VL(mzd_mul_v_neon_256, mzd_mul_vl_neon_256)
#define ADDMUL SELECT_V_VL(mzd_addmul_v_neon_256, mzd_addmul_vl_neon_256)

#undef LOWMC_N
#undef LOWMC_R
#define LOWMC_N LOWMC_L5_N
#define LOWMC_R LOWMC_L5_R

#undef SBOX_IMPL
#undef LOWMC
#define SBOX_IMPL sbox_layer_bitsliced
#define LOWMC lowmc_neon_256
#include "lowmc.c.i"

#if defined(WITH_CUSTOM_INSTANCES)
#undef XOR
#undef MUL
#undef ADDMUL
#define XOR mzd_xor_neon
#define MUL SELECT_V_VL(mzd_mul_v_neon, mzd_mul_vl_neon)
#define ADDMUL SELECT_V_VL(mzd_addmul_v_neon, mzd_addmul_vl_neon)

#undef LOWMC_INSTANCE
#undef LOWMC_N
#undef LOWMC_R
#define LOWMC_N lowmc->n
#define LOWMC_R lowmc->r

// generic using NEON
#undef SBOX_IMPL
#undef LOWMC
#define SBOX_IMPL sbox_layer_bitsliced
#define LOWMC lowmc_neon
#include "lowmc.c.i"
#endif
#endif
#endif

lowmc_implementation_f lowmc_get_implementation(const lowmc_t* lowmc) {
#if defined(WITH_OPT)
#if defined(WITH_AVX2)
  if (CPU_SUPPORTS_AVX2) {
    if(lowmc->m == 10) {
      switch (lowmc->n) {
        case 128:
          return general_or_10(lowmc, lowmc_avx_128);
        case 192:
          return general_or_10(lowmc, lowmc_avx_192);
        case 256:
          return general_or_10(lowmc, lowmc_avx_256);
      }
    }
    if(lowmc->m == 1) {
      switch (lowmc->n) {
        case 128:
          return general_or_1(lowmc, lowmc_avx_128);
        case 192:
          return general_or_1(lowmc, lowmc_avx_192);
        case 256:
          return general_or_1(lowmc, lowmc_avx_256);
      }
    }
#if defined(WITH_CUSTOM_INSTANCES)
    if (lowmc->n > 256) {
      return general_or_10(lowmc, lowmc_avx);
    }
#endif
  }
#endif
#if defined(WITH_SSE2)
  if (CPU_SUPPORTS_SSE2) {
    if(lowmc->m == 10) {
      switch (lowmc->n) {
        case 128:
          return general_or_10(lowmc, lowmc_sse_128);
        case 192:
          return general_or_10(lowmc, lowmc_sse_192);
        case 256:
          return general_or_10(lowmc, lowmc_sse_256);
      }
    }
    if(lowmc->m == 1) {
      switch (lowmc->n) {
        case 128:
          return general_or_1(lowmc, lowmc_sse_128);
        case 192:
          return general_or_1(lowmc, lowmc_sse_192);
        case 256:
          return general_or_1(lowmc, lowmc_sse_256);
      }
    }
#if defined(WITH_CUSTOM_INSTANCES)
    if (lowmc->n > 256) {
      return general_or_10(lowmc, lowmc_sse);
    }
#endif
  }
#endif
#if defined(WITH_NEON)
  if (CPU_SUPPORTS_NEON) {
    if(lowmc->m == 10) {
      switch (lowmc->n) {
        case 128:
          return general_or_10(lowmc, lowmc_neon_128);
        case 192:
          return general_or_10(lowmc, lowmc_neon_192);
        case 256:
          return general_or_10(lowmc, lowmc_neon_256);
      }
    }
    if(lowmc->m == 1) {
      switch (lowmc->n) {
        case 128:
          return general_or_1(lowmc, lowmc_neon_128);
        case 192:
          return general_or_1(lowmc, lowmc_neon_192);
        case 256:
          return general_or_1(lowmc, lowmc_neon_256);
      }
    }
#if defined(WITH_CUSTOM_INSTANCES)
    if (lowmc->n > 256) {
      return general_or_10(lowmc, lowmc_neon);
    }
#endif
  }
#endif
#endif

  (void)lowmc;
  if(lowmc->m == 10)
    return general_or_10(lowmc, lowmc_uint64);
  else if (lowmc->m == 1)
    return general_or_1(lowmc, lowmc_uint64);
  else
    return NULL;
}

lowmc_store_implementation_f lowmc_store_get_implementation(const lowmc_t* lowmc) {
#if defined(WITH_OPT)
#if defined(WITH_AVX2)
  if (CPU_SUPPORTS_AVX2) {
    switch (lowmc->n) {
    case 128:
      return general_or_10(lowmc_store, lowmc_avx_128_store);
    case 192:
      return general_or_10(lowmc_store, lowmc_avx_192_store);
    case 256:
      return general_or_10(lowmc_store, lowmc_avx_256_store);
    }
#if defined(WITH_CUSTOM_INSTANCES)
    if (lowmc->n > 256) {
      return general_or_10(lowmc_store, lowmc_avx_store);
    }
#endif
  }
#endif
#if defined(WITH_SSE2)
  if (CPU_SUPPORTS_SSE2) {
    switch (lowmc->n) {
    case 128:
      return general_or_10(lowmc_store, lowmc_sse_128_store);
    case 192:
      return general_or_10(lowmc_store, lowmc_sse_192_store);
    case 256:
      return general_or_10(lowmc_store, lowmc_sse_256_store);
    }
#if defined(WITH_CUSTOM_INSTANCES)
    if (lowmc->n > 256) {
      return general_or_10(lowmc_store, lowmc_sse_store);
    }
#endif
  }
#endif
#if defined(WITH_NEON)
  if (CPU_SUPPORTS_NEON) {
    switch (lowmc->n) {
    case 128:
      return general_or_10(lowmc_store, lowmc_neon_128_store);
    case 192:
      return general_or_10(lowmc_store, lowmc_neon_192_store);
    case 256:
      return general_or_10(lowmc_store, lowmc_neon_256_store);
    }
#if defined(WITH_CUSTOM_INSTANCES)
    if (lowmc->n > 256) {
      return general_or_10(lowmc_store, lowmc_neon_store);
    }
#endif
  }
#endif
#endif

  (void)lowmc;
  return general_or_10(lowmc_store, lowmc_uint64_store);
}

mzd_local_t* lowmc_call(lowmc_t const* lowmc, lowmc_key_t const* lowmc_key, mzd_local_t const* p) {
  lowmc_implementation_f impl = lowmc_get_implementation(lowmc);
  return impl(lowmc, lowmc_key, p);
}
