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

#ifdef WITH_OPT
#include "simd.h"
#endif

#if !defined(_MSC_VER)
#include <stdalign.h>
#endif
#include <string.h>

static uint64_t sbox_layer_bitsliced_uint64(uint64_t in) {
  // a, b, c
  const uint64_t x0m = (in & MASK_X0I) << 2;
  const uint64_t x1m = (in & MASK_X1I) << 1;
  const uint64_t x2m = in & MASK_X2I;

  // (b & c) ^ a
  const uint64_t t0 = (x1m & x2m) ^ x0m;
  // (c & a) ^ a ^ b
  const uint64_t t1 = (x0m & x2m) ^ x0m ^ x1m;
  // (a & b) ^ a ^ b ^c
  const uint64_t t2 = (x0m & x1m) ^ x0m ^ x1m ^ x2m;

  return (in & MASK_MASK) ^ (t0 >> 2) ^ (t1 >> 1) ^ t2;
}

/**
 * S-box for m = 10
 */
static void sbox_layer_uint64(mzd_local_t* x, mask_t const* mask) {
  (void)mask;

  uint64_t* d = &FIRST_ROW(x)[x->width - 1];
  *d = sbox_layer_bitsliced_uint64(*d);
}

#ifdef WITH_CUSTOM_INSTANCES
static void sbox_layer_bitsliced(mzd_local_t* in, mask_t const* mask) {
  mzd_local_t* buffer[6] = {NULL};
  mzd_local_init_multiple_ex(buffer, 6, 1, in->ncols, false);

  mzd_local_t* x0m = buffer[0];
  mzd_local_t* x1m = buffer[1];
  mzd_local_t* x2m = buffer[2];
  mzd_local_t* t0 = buffer[3];
  mzd_local_t* t1 = buffer[4];
  mzd_local_t* t2 = buffer[5];

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

#ifdef WITH_OPT
#ifdef WITH_SSE2
ATTR_TARGET("sse") static void sbox_layer_sse(mzd_local_t* in,
                                                           mask_t const* mask) {
  __m128i* ip = (__m128i*)ASSUME_ALIGNED(CONST_FIRST_ROW(in), alignof(__m128i));
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

#ifdef WITH_AVX2
/**
 * AVX2 version of LowMC. It assumes that mzd_local_t's row[0] is always 32 byte
 * aligned.
 */
ATTR_TARGET("avx2") static void sbox_layer_avx(mzd_local_t* in, mask_t const* mask) {
  __m256i* ip = (__m256i*)ASSUME_ALIGNED(CONST_FIRST_ROW(in), alignof(__m256i));
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

#ifdef WITH_NEON
static void sbox_layer_neon(mzd_local_t* in, mask_t const* mask) {
  uint32x4_t* ip =
      (uint32x4_t*)ASSUME_ALIGNED(CONST_FIRST_ROW(in), alignof(uint32x4_t));
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

typedef void (*sbox_layer_impl)(mzd_local_t*, mask_t const*);

static sbox_layer_impl get_sbox_layer(const lowmc_t* lowmc) {
  if (lowmc->m == 10) {
    return sbox_layer_uint64;
  }
#ifdef WITH_CUSTOM_INSTANCES
#ifdef WITH_OPT
#ifdef WITH_SSE2
  if (CPU_SUPPORTS_SSE2 && lowmc->n == 128) {
    return sbox_layer_sse;
  }
#endif
#ifdef WITH_AVX2
  if (CPU_SUPPORTS_AVX2 && lowmc->n == 256) {
    return sbox_layer_avx;
  }
#endif
#ifdef WITH_NEON
  if (CPU_SUPPORTS_NEON && lowmc->n == 128) {
    return sbox_layer_neon;
  }
#endif
#endif
  return sbox_layer_bitsliced;
#else
  return NULL;
#endif
}

#if defined(REDUCED_LINEAR_LAYER)
static mzd_local_t* lowmc_reduced_linear_layer(lowmc_t const* lowmc, lowmc_key_t const* lowmc_key,
                                               mzd_local_t const* p) {
  mzd_local_t* x       = mzd_local_init_ex(1, lowmc->n, false);
  mzd_local_t* y       = mzd_local_init_ex(1, lowmc->n, false);
  mzd_local_t* nl_part = mzd_local_init_ex(1, lowmc->r * 32, false);

  mzd_local_copy(x, p);
#if defined(MUL_M4RI)
  mzd_addmul_vl(x, lowmc_key, lowmc->k0_lookup);
  mzd_mul_vl(nl_part, lowmc_key, lowmc->precomputed_non_linear_part_lookup);
#else
  mzd_addmul_v(x, lowmc_key, lowmc->k0_matrix);
  mzd_mul_v(nl_part, lowmc_key, lowmc->precomputed_non_linear_part_matrix);
#endif

  lowmc_round_t const* round = lowmc->rounds;
  for (unsigned i = 0; i < lowmc->r; ++i, ++round) {
    sbox_layer_uint64(x, NULL);

    const word mask          = (i & 1) ? WORD_C(0xFFFFFFFF00000000) : WORD_C(0x00000000FFFFFFFF);
    const unsigned int shift = (i & 1) ? 2 : 34;

    FIRST_ROW(x)[x->width - 1] ^= (CONST_FIRST_ROW(nl_part)[i >> 1] & mask) << shift;

#if defined(MUL_M4RI)
    mzd_mul_vl(y, x, round->l_lookup);
#else
    mzd_mul_v(y, x, round->l_matrix);
#endif
    mzd_xor(x, y, round->constant);
  }

  mzd_local_free(y);
  mzd_local_free(nl_part);
  return x;
}
#else
static mzd_local_t* lowmc_plain(lowmc_t const* lowmc, lowmc_key_t const* lowmc_key,
                                mzd_local_t const* p, sbox_layer_impl sbox_layer) {
  (void)sbox_layer;

  mzd_local_t* x = mzd_local_init_ex(1, lowmc->n, false);
  mzd_local_t* y = mzd_local_init_ex(1, lowmc->n, false);

  mzd_local_copy(x, p);
#if defined(MUL_M4RI)
  mzd_addmul_vl(x, lowmc_key, lowmc->k0_lookup);
#else
  mzd_addmul_v(x, lowmc_key, lowmc->k0_matrix);
#endif

  lowmc_round_t const* round = lowmc->rounds;
  for (unsigned int i = lowmc->r; i; --i, ++round) {
#if defined(WITH_CUSTOM_INSTANCE)
    sbox_layer(x, &lowmc->mask);
#else
    sbox_layer_uint64(x, NULL);
#endif

#if defined(MUL_M4RI)
    mzd_mul_vl(y, x, round->l_lookup);
#else
    mzd_mul_v(y, x, round->l_matrix);
#endif
    mzd_xor(x, y, round->constant);
#if defined(MUL_M4RI) && !defined(REDUCED_LINEAR_LAYER)
    mzd_addmul_vl(x, lowmc_key, round->k_lookup);
#else
    mzd_addmul_v(x, lowmc_key, round->k_matrix);
#endif
  }

  mzd_local_free(y);
  return x;
}
#endif

mzd_local_t* lowmc_call(lowmc_t const* lowmc, lowmc_key_t const* lowmc_key, mzd_local_t const* p) {
  sbox_layer_impl sbox_layer = get_sbox_layer(lowmc);
  if (!sbox_layer) {
    return NULL;
  }

#if defined(REDUCED_LINEAR_LAYER)
  if (lowmc->m == 10) {
    return lowmc_reduced_linear_layer(lowmc, lowmc_key, p);
  }
  return NULL;
#else
  return lowmc_plain(lowmc, lowmc_key, p, sbox_layer);
#endif
}
