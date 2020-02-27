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
#include "mzd_additional.h"
#if defined(WITH_KKW)
#include "picnic2_impl.h"
#include "picnic2_simulate_mul.h"
#endif
#if defined(WITH_OPT)
#include "simd.h"
#endif

#if !defined(_MSC_VER)
#include <stdalign.h>
#endif
#include <string.h>
#include <assert.h>

/**
 * S-box for m = 42
 */
static void sbox_uint64_lowmc_126_126_4(mzd_local_t* in) {
  mzd_local_t x0m[1], x1m[1], x2m[1];
  // a
  mzd_and_uint64_128(x0m, mask_126_126_42_a, in);
  // b
  mzd_and_uint64_128(x1m, mask_126_126_42_b, in);
  // c
  mzd_and_uint64_128(x2m, mask_126_126_42_c, in);

  mzd_shift_left_uint64_128(x0m, x0m, 2);
  mzd_shift_left_uint64_128(x1m, x1m, 1);

  mzd_local_t t0[1], t1[1], t2[1];
  // b & c
  mzd_and_uint64_128(t0, x1m, x2m);
  // c & a
  mzd_and_uint64_128(t1, x0m, x2m);
  // a & b
  mzd_and_uint64_128(t2, x0m, x1m);

  // (b & c) ^ a
  mzd_xor_uint64_128(t0, t0, x0m);

  // (c & a) ^ a ^ b
  mzd_xor_uint64_128(t1, t1, x0m);
  mzd_xor_uint64_128(t1, t1, x1m);

  // (a & b) ^ a ^ b ^c
  mzd_xor_uint64_128(t2, t2, x0m);
  mzd_xor_uint64_128(t2, t2, x1m);
  mzd_xor_uint64_128(t2, t2, x2m);

  mzd_shift_right_uint64_128(t0, t0, 2);
  mzd_shift_right_uint64_128(t1, t1, 1);

  mzd_xor_uint64_128(t2, t2, t1);
  mzd_xor_uint64_128(in, t2, t0);
}

/**
 * S-box for m = 64
 */
static void sbox_uint64_lowmc_192_192_4(mzd_local_t* in) {
  mzd_local_t x0m[1], x1m[1], x2m[1];
  // a
  mzd_and_uint64_192(x0m, mask_192_192_64_a, in);
  // b
  mzd_and_uint64_192(x1m, mask_192_192_64_b, in);
  // c
  mzd_and_uint64_192(x2m, mask_192_192_64_c, in);

  mzd_rotate_left_uint64_192(x0m, x0m, 2);
  mzd_rotate_left_uint64_192(x1m, x1m, 1);

  mzd_local_t t0[1], t1[1], t2[1];
  // b & c
  mzd_and_uint64_192(t0, x1m, x2m);
  // c & a
  mzd_and_uint64_192(t1, x0m, x2m);
  // a & b
  mzd_and_uint64_192(t2, x0m, x1m);

  // (b & c) ^ a
  mzd_xor_uint64_192(t0, t0, x0m);

  // (c & a) ^ a ^ b
  mzd_xor_uint64_192(t1, t1, x0m);
  mzd_xor_uint64_192(t1, t1, x1m);

  // (a & b) ^ a ^ b ^c
  mzd_xor_uint64_192(t2, t2, x0m);
  mzd_xor_uint64_192(t2, t2, x1m);
  mzd_xor_uint64_192(t2, t2, x2m);

  mzd_rotate_right_uint64_192(t0, t0, 2);
  mzd_rotate_right_uint64_192(t1, t1, 1);

  mzd_xor_uint64_192(t2, t2, t1);
  mzd_xor_uint64_192(in, t2, t0);
}

/**
 * S-box for m = 85
 */
static void sbox_uint64_lowmc_255_255_4(mzd_local_t* in) {
  mzd_local_t x0m[1], x1m[1], x2m[1];
  // a
  mzd_and_uint64_256(x0m, mask_255_255_85_a, in);
  // b
  mzd_and_uint64_256(x1m, mask_255_255_85_b, in);
  // c
  mzd_and_uint64_256(x2m, mask_255_255_85_c, in);

  mzd_rotate_left_uint64_256(x0m, x0m, 2);
  mzd_rotate_left_uint64_256(x1m, x1m, 1);

  mzd_local_t t0[1], t1[1], t2[1];
  // b & c
  mzd_and_uint64_256(t0, x1m, x2m);
  // c & a
  mzd_and_uint64_256(t1, x0m, x2m);
  // a & b
  mzd_and_uint64_256(t2, x0m, x1m);

  // (b & c) ^ a
  mzd_xor_uint64_256(t0, t0, x0m);

  // (c & a) ^ a ^ b
  mzd_xor_uint64_256(t1, t1, x0m);
  mzd_xor_uint64_256(t1, t1, x1m);

  // (a & b) ^ a ^ b ^c
  mzd_xor_uint64_256(t2, t2, x0m);
  mzd_xor_uint64_256(t2, t2, x1m);
  mzd_xor_uint64_256(t2, t2, x2m);

  mzd_rotate_right_uint64_256(t0, t0, 2);
  mzd_rotate_right_uint64_256(t1, t1, 1);

  mzd_xor_uint64_256(t2, t2, t1);
  mzd_xor_uint64_256(in, t2, t0);
}

#if defined(WITH_OPT)
#if defined(WITH_SSE2) || defined(WITH_NEON)
ATTR_TARGET_S128
static inline void sbox_s128_lowmc_126_126_4(mzd_local_t* in) {
  const word128 min = in->w128[0];

  word128 x0m = mm128_and(min, mask_126_126_42_a->w128[0]);
  word128 x1m = mm128_and(min, mask_126_126_42_b->w128[0]);
  word128 x2m = mm128_and(min, mask_126_126_42_c->w128[0]);

  x0m = mm128_shift_left(x0m, 2);
  x1m = mm128_shift_left(x1m, 1);

  word128 t0 = mm128_and(x1m, x2m);
  word128 t1 = mm128_and(x0m, x2m);
  word128 t2 = mm128_and(x0m, x1m);

  t0 = mm128_xor(t0, x0m);

  x0m = mm128_xor(x0m, x1m);
  t1  = mm128_xor(t1, x0m);

  t2 = mm128_xor(t2, x0m);
  t2 = mm128_xor(t2, x2m);

  t0 = mm128_shift_right(t0, 2);
  t1 = mm128_shift_right(t1, 1);

  in->w128[0] = mm128_xor(mm128_xor(t0, t1), t2);
}

#define sbox_s128_lowmc_192_192_4 sbox_uint64_lowmc_192_192_4
#define sbox_s128_lowmc_255_255_4 sbox_uint64_lowmc_255_255_4
#endif

#if defined(WITH_AVX2)
ATTR_TARGET_AVX2
static inline word256 sbox_s256_lowmc_full(const word256 min, const word256 mask_a,
                                           const word256 mask_b, const word256 mask_c) {
  word256 x0m = mm256_and(min, mask_a);
  word256 x1m = mm256_and(min, mask_b);
  word256 x2m = mm256_and(min, mask_c);

  x0m = mm256_rotate_left(x0m, 2);
  x1m = mm256_rotate_left(x1m, 1);

  word256 t0 = mm256_and(x1m, x2m);
  word256 t1 = mm256_and(x0m, x2m);
  word256 t2 = mm256_and(x0m, x1m);

  t0 = mm256_xor(t0, x0m);

  x0m = mm256_xor(x0m, x1m);
  t1  = mm256_xor(t1, x0m);

  t2 = mm256_xor(t2, x0m);
  t2 = mm256_xor(t2, x2m);

  t0 = mm256_rotate_right(t0, 2);
  t1 = mm256_rotate_right(t1, 1);

  return mm256_xor(mm256_xor(t0, t1), t2);
}

/* TODO: implement with SSE? */
ATTR_TARGET_AVX2
static inline void sbox_s256_lowmc_126_126_4(mzd_local_t* in) {
  in->w256 = sbox_s256_lowmc_full(in->w256, mask_126_126_42_a->w256, mask_126_126_42_b->w256,
                                  mask_126_126_42_c->w256);
}

ATTR_TARGET_AVX2
static inline void sbox_s256_lowmc_192_192_4(mzd_local_t* in) {
  in->w256 = sbox_s256_lowmc_full(in->w256, mask_192_192_64_a->w256, mask_192_192_64_b->w256,
                                  mask_192_192_64_c->w256);
}

ATTR_TARGET_AVX2
static inline void sbox_s256_lowmc_255_255_4(mzd_local_t* in) {
  in->w256 = sbox_s256_lowmc_full(in->w256, mask_255_255_85_a->w256, mask_255_255_85_b->w256,
                                  mask_255_255_85_c->w256);
}
#endif /* WITH_AVX2 */
#endif /* WITH_OPT */

#if defined(WITH_LOWMC_126_126_4)
#include "lowmc_126_126_4.h"
#endif
#if defined(WITH_LOWMC_192_192_4)
#include "lowmc_192_192_4.h"
#endif
#if defined(WITH_LOWMC_255_255_4)
#include "lowmc_255_255_4.h"
#endif

#if !defined(NO_UINT64_FALLBACK)
// uint64 based implementation
#define IMPL uint64
#include "lowmc_fns_uint64_L1.h"
#include "lowmc.c.i"

#include "lowmc_fns_uint64_L3.h"
#include "lowmc.c.i"

#include "lowmc_fns_uint64_L5.h"
#include "lowmc.c.i"
#endif

#if defined(WITH_OPT)
#if defined(WITH_SSE2) || defined(WITH_NEON)
#if defined(WITH_SSE2)
#define FN_ATTR ATTR_TARGET_SSE2
#endif
#undef IMPL
#define IMPL s128

// L1 using SSE2/NEON
#include "lowmc_fns_s128_L1.h"
#include "lowmc.c.i"

// L3 using SSE2/NEON
#include "lowmc_fns_s128_L3.h"
#include "lowmc.c.i"

// L5 using SSE2/NEON
#include "lowmc_fns_s128_L5.h"
#include "lowmc.c.i"

#undef FN_ATTR
#endif

#if defined(WITH_AVX2)
#define FN_ATTR ATTR_TARGET_AVX2
#undef IMPL
#define IMPL s256

// L1 using AVX2
#include "lowmc_fns_s256_L1.h"
#include "lowmc.c.i"

// L3 using AVX2
#include "lowmc_fns_s256_L3.h"
#include "lowmc.c.i"

// L5 using AVX2
#include "lowmc_fns_s256_L5.h"
#include "lowmc.c.i"

#undef FN_ATTR

#endif
#endif

lowmc_implementation_f lowmc_get_implementation(const lowmc_t* lowmc) {
  assert(lowmc->m == 42 || lowmc->m == 64 || lowmc->m == 85);
  assert(lowmc->n == 126 || lowmc->n == 192 || lowmc->n == 255);

#if defined(WITH_OPT)
#if defined(WITH_AVX2)
  if (CPU_SUPPORTS_AVX2) {
#if defined(WITH_LOWMC_126_126_4)
    if (lowmc->n == 126 && lowmc->m == 42)
      return lowmc_s256_lowmc_126_126_4;
#endif
#if defined(WITH_LOWMC_192_192_4)
    if (lowmc->n == 192 && lowmc->m == 64)
      return lowmc_s256_lowmc_192_192_4;
#endif
#if defined(WITH_LOWMC_255_255_4)
    if (lowmc->n == 255 && lowmc->m == 85)
      return lowmc_s256_lowmc_255_255_4;
#endif
  }
#endif
#if defined(WITH_SSE2) || defined(WITH_NEON)
  if (CPU_SUPPORTS_SSE2 || CPU_SUPPORTS_NEON) {
#if defined(WITH_LOWMC_126_126_4)
    if (lowmc->n == 126 && lowmc->m == 42)
      return lowmc_s128_lowmc_126_126_4;
#endif
#if defined(WITH_LOWMC_192_192_4)
    if (lowmc->n == 192 && lowmc->m == 64)
      return lowmc_s128_lowmc_192_192_4;
#endif
#if defined(WITH_LOWMC_255_255_4)
    if (lowmc->n == 255 && lowmc->m == 85)
      return lowmc_s128_lowmc_255_255_4;
#endif
  }
#endif
#endif

#if !defined(NO_UINT64_FALLBACK)
#if defined(WITH_LOWMC_126_126_4)
  if (lowmc->n == 126 && lowmc->m == 42)
    return lowmc_uint64_lowmc_126_126_4;
#endif
#if defined(WITH_LOWMC_192_192_4)
  if (lowmc->n == 192 && lowmc->m == 64)
    return lowmc_uint64_lowmc_192_192_4;
#endif
#if defined(WITH_LOWMC_255_255_4)
  if (lowmc->n == 255 && lowmc->m == 85)
    return lowmc_uint64_lowmc_255_255_4;
#endif
#endif

  return NULL;
}

#if defined(WITH_ZKBPP)
lowmc_store_implementation_f lowmc_store_get_implementation(const lowmc_t* lowmc) {
  assert(lowmc->m == 42 || lowmc->m == 64 || lowmc->m == 85);
  assert(lowmc->n == 126 || lowmc->n == 192 || lowmc->n == 255);

#if defined(WITH_OPT)
#if defined(WITH_AVX2)
  if (CPU_SUPPORTS_AVX2) {
#if defined(WITH_LOWMC_126_126_4)
    if (lowmc->n == 126 && lowmc->m == 42)
      return lowmc_store_s256_lowmc_126_126_4;
#endif
#if defined(WITH_LOWMC_192_192_4)
    if (lowmc->n == 192 && lowmc->m == 64)
      return lowmc_store_s256_lowmc_192_192_4;
#endif
#if defined(WITH_LOWMC_255_255_4)
    if (lowmc->n == 255 && lowmc->m == 85)
      return lowmc_store_s256_lowmc_255_255_4;
#endif
  }
#endif
#if defined(WITH_SSE2) || defined(WITH_NEON)
  if (CPU_SUPPORTS_SSE2 || CPU_SUPPORTS_NEON) {
#if defined(WITH_LOWMC_126_126_4)
    if (lowmc->n == 126 && lowmc->m == 42)
      return lowmc_store_s128_lowmc_126_126_4;
#endif
#if defined(WITH_LOWMC_192_192_4)
    if (lowmc->n == 192 && lowmc->m == 64)
      return lowmc_store_s128_lowmc_192_192_4;
#endif
#if defined(WITH_LOWMC_255_255_4)
    if (lowmc->n == 255 && lowmc->m == 85)
      return lowmc_store_s128_lowmc_255_255_4;
#endif
  }
#endif
#endif

#if !defined(NO_UINT64_FALLBACK)
#if defined(WITH_LOWMC_126_126_4)
  if (lowmc->n == 126 && lowmc->m == 42)
    return lowmc_store_uint64_lowmc_126_126_4;
#endif
#if defined(WITH_LOWMC_192_192_4)
  if (lowmc->n == 192 && lowmc->m == 64)
    return lowmc_store_uint64_lowmc_192_192_4;
#endif
#if defined(WITH_LOWMC_255_255_4)
  if (lowmc->n == 255 && lowmc->m == 85)
    return lowmc_store_uint64_lowmc_255_255_4;
#endif
#endif

  return NULL;
}
#endif

#if defined(WITH_KKW)
lowmc_compute_aux_implementation_f lowmc_compute_aux_get_implementation(const lowmc_t* lowmc) {
  assert(lowmc->m == 42 || lowmc->m == 64 || lowmc->m == 85);
  assert(lowmc->n == 126 || lowmc->n == 192 || lowmc->n == 255);

#if defined(WITH_OPT)
#if defined(WITH_AVX2)
  if (CPU_SUPPORTS_AVX2) {
#if defined(WITH_LOWMC_126_126_4)
    if (lowmc->n == 126 && lowmc->m == 42)
      return lowmc_compute_aux_s256_lowmc_126_126_4;
#endif
#if defined(WITH_LOWMC_192_192_4)
    if (lowmc->n == 192 && lowmc->m == 64)
      return lowmc_compute_aux_s256_lowmc_192_192_4;
#endif
#if defined(WITH_LOWMC_255_255_4)
    if (lowmc->n == 255 && lowmc->m == 85)
      return lowmc_compute_aux_s256_lowmc_255_255_4;
#endif
  }
#endif
#if defined(WITH_SSE2) || defined(WITH_NEON)
  if (CPU_SUPPORTS_SSE2 || CPU_SUPPORTS_NEON) {
#if defined(WITH_LOWMC_126_126_4)
    if (lowmc->n == 126 && lowmc->m == 42)
      return lowmc_compute_aux_s128_lowmc_126_126_4;
#endif
#if defined(WITH_LOWMC_192_192_4)
    if (lowmc->n == 192 && lowmc->m == 64)
      return lowmc_compute_aux_s128_lowmc_192_192_4;
#endif
#if defined(WITH_LOWMC_255_255_4)
    if (lowmc->n == 255 && lowmc->m == 85)
      return lowmc_compute_aux_s128_lowmc_255_255_4;
#endif
  }
#endif
#endif

#if !defined(NO_UINT64_FALLBACK)
#if defined(WITH_LOWMC_126_126_4)
  if (lowmc->n == 126 && lowmc->m == 42)
    return lowmc_compute_aux_uint64_lowmc_126_126_4;
#endif
#if defined(WITH_LOWMC_192_192_4)
  if (lowmc->n == 192 && lowmc->m == 64)
    return lowmc_compute_aux_uint64_lowmc_192_192_4;
#endif
#if defined(WITH_LOWMC_255_255_4)
  if (lowmc->n == 255 && lowmc->m == 85)
    return lowmc_compute_aux_uint64_lowmc_255_255_4;
#endif
#endif

  return NULL;
}
#endif
