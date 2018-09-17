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

#if defined(WITH_CUSTOM_INSTANCES)
typedef struct {
  mzd_local_t* x0m[SC_PROOF]; // a
  mzd_local_t* x1m[SC_PROOF]; // b
  mzd_local_t* x2m[SC_PROOF]; // c
  mzd_local_t* r0m[SC_PROOF];
  mzd_local_t* r1m[SC_PROOF];
  mzd_local_t* r2m[SC_PROOF];
  mzd_local_t* x0s[SC_PROOF];
  mzd_local_t* r0s[SC_PROOF];
  mzd_local_t* x1s[SC_PROOF];
  mzd_local_t* r1s[SC_PROOF];
  mzd_local_t* v[SC_PROOF];

  mzd_local_t** storage;
} sbox_vars_t;

static void sbox_vars_clear(sbox_vars_t* vars) {
  if (vars->storage) {
    mzd_local_free_multiple(vars->storage);
    free(vars->storage);
    memset(vars, 0, sizeof(*vars));
  }
}

static sbox_vars_t* sbox_vars_init(sbox_vars_t* vars, uint32_t n, unsigned sc) {
  vars->storage = calloc(11 * sc, sizeof(mzd_local_t*));
  mzd_local_init_multiple_ex(vars->storage, 11 * sc, 1, n, false);

  for (unsigned int i = 0; i < sc; ++i) {
    vars->x0m[i] = vars->storage[11 * i + 0];
    vars->x1m[i] = vars->storage[11 * i + 1];
    vars->x2m[i] = vars->storage[11 * i + 2];
    vars->r0m[i] = vars->storage[11 * i + 3];
    vars->r1m[i] = vars->storage[11 * i + 4];
    vars->r2m[i] = vars->storage[11 * i + 5];
    vars->x0s[i] = vars->storage[11 * i + 6];
    vars->x1s[i] = vars->storage[11 * i + 7];
    vars->r0s[i] = vars->storage[11 * i + 8];
    vars->r1s[i] = vars->storage[11 * i + 9];
    vars->v[i]   = vars->storage[11 * i + 10];
  }

  return vars;
}

#define bitsliced_step_1(sc)                                                                       \
  mpc_and_const(out, in, mask->mask, sc);                                                          \
                                                                                                   \
  mpc_and_const(vars->x0m, in, mask->x0, sc);                                                      \
  mpc_and_const(vars->x1m, in, mask->x1, sc);                                                      \
  mpc_and_const(vars->x2m, in, mask->x2, sc);                                                      \
  mpc_and_const(vars->r0m, rvec, mask->x0, sc);                                                    \
  mpc_and_const(vars->r1m, rvec, mask->x1, sc);                                                    \
  mpc_and_const(vars->r2m, rvec, mask->x2, sc);                                                    \
                                                                                                   \
  mpc_shift_left(vars->x0s, vars->x0m, 2, sc);                                                     \
  mpc_shift_left(vars->r0s, vars->r0m, 2, sc);                                                     \
                                                                                                   \
  mpc_shift_left(vars->x1s, vars->x1m, 1, sc);                                                     \
  mpc_shift_left(vars->r1s, vars->r1m, 1, sc)

#define bitsliced_step_2(sc)                                                                       \
  /* (b & c) ^ a */                                                                                \
  mpc_xor(vars->r2m, vars->r2m, vars->x0s, sc);                                                    \
  /* a ^ b */                                                                                      \
  mpc_xor(vars->x0s, vars->x0s, vars->x1s, sc);                                                    \
  /* (c & a) ^ a ^ b */                                                                            \
  mpc_xor(vars->r1m, vars->r1m, vars->x0s, sc);                                                    \
  /* (a & b) ^ a ^ b ^ c */                                                                        \
  mpc_xor(vars->r0m, vars->r0m, vars->x0s, sc);                                                    \
  mpc_xor(vars->r0m, vars->r0m, vars->x2m, sc);                                                    \
                                                                                                   \
  mpc_shift_right(vars->x0s, vars->r2m, 2, sc);                                                    \
  mpc_shift_right(vars->x1s, vars->r1m, 1, sc);                                                    \
                                                                                                   \
  mpc_xor(out, out, vars->r0m, sc);                                                                \
  mpc_xor(out, out, vars->x0s, sc);                                                                \
  mpc_xor(out, out, vars->x1s, sc)

static void mpc_sbox_layer_bitsliced(mzd_local_t** out, mzd_local_t* const* in, view_t* view,
                                     mzd_local_t* const* rvec, mask_t const* mask,
                                     sbox_vars_t const* vars) {
  bitsliced_step_1(SC_PROOF);

  mpc_clear(view->s, SC_PROOF);
  // a & b
  mpc_and(vars->r0m, vars->x0s, vars->x1s, vars->r2m, view, 0, vars->v);
  // b & c
  mpc_and(vars->r2m, vars->x1s, vars->x2m, vars->r1s, view, 1, vars->v);
  // c & a
  mpc_and(vars->r1m, vars->x0s, vars->x2m, vars->r0s, view, 2, vars->v);

  bitsliced_step_2(SC_PROOF);
}

static void mpc_sbox_layer_bitsliced_verify(mzd_local_t** out, mzd_local_t* const* in, view_t* view,
                                            mzd_local_t* const* rvec, mask_t const* mask,
                                            sbox_vars_t const* vars) {
  bitsliced_step_1(SC_VERIFY);

  mzd_local_clear(view->s[0]);
  // a & b
  mpc_and_verify(vars->r0m, vars->x0s, vars->x1s, vars->r2m, view, mask->x2, 0, vars->v);
  // b & c
  mpc_and_verify(vars->r2m, vars->x1s, vars->x2m, vars->r1s, view, mask->x2, 1, vars->v);
  // c & a
  mpc_and_verify(vars->r1m, vars->x0s, vars->x2m, vars->r0s, view, mask->x2, 2, vars->v);

  bitsliced_step_2(SC_VERIFY);
}
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

#if defined(WITH_OPT) && defined(WITH_CUSTOM_INSTANCES)
#define bitsliced_mm_step_1(sc, type, and, shift_left)                                             \
  type r0m[sc] ATTR_ALIGNED(alignof(type));                                                        \
  type r0s[sc] ATTR_ALIGNED(alignof(type));                                                        \
  type r1m[sc] ATTR_ALIGNED(alignof(type));                                                        \
  type r1s[sc] ATTR_ALIGNED(alignof(type));                                                        \
  type r2m[sc] ATTR_ALIGNED(alignof(type));                                                        \
  type x0s[sc] ATTR_ALIGNED(alignof(type));                                                        \
  type x1s[sc] ATTR_ALIGNED(alignof(type));                                                        \
  type x2m[sc] ATTR_ALIGNED(alignof(type));                                                        \
  const type mx2 ATTR_ALIGNED(alignof(type)) =                                                     \
      *((const type*)ASSUME_ALIGNED(CONST_FIRST_ROW(mask->x2), alignof(type)));                    \
  do {                                                                                             \
    const type mx0 ATTR_ALIGNED(alignof(type)) =                                                   \
        *((const type*)ASSUME_ALIGNED(CONST_FIRST_ROW(mask->x0), alignof(type)));                  \
    const type mx1 ATTR_ALIGNED(alignof(type)) =                                                   \
        *((const type*)ASSUME_ALIGNED(CONST_FIRST_ROW(mask->x1), alignof(type)));                  \
                                                                                                   \
    for (unsigned int m = 0; m < (sc); ++m) {                                                      \
      const type inm ATTR_ALIGNED(alignof(type)) =                                                 \
          *((const type*)ASSUME_ALIGNED(CONST_FIRST_ROW(in[m]), alignof(type)));                   \
      const type rvecm ATTR_ALIGNED(alignof(type)) =                                               \
          *((const type*)ASSUME_ALIGNED(CONST_FIRST_ROW(rvec[m]), alignof(type)));                 \
                                                                                                   \
      type tmp1 = (and)(inm, mx0);                                                                 \
      type tmp2 = (and)(inm, mx1);                                                                 \
      x2m[m]    = (and)(inm, mx2);                                                                 \
                                                                                                   \
      x0s[m] = (shift_left)(tmp1, 2);                                                              \
      x1s[m] = (shift_left)(tmp2, 1);                                                              \
                                                                                                   \
      r0m[m] = tmp1 = (and)(rvecm, mx0);                                                           \
      r1m[m] = tmp2 = (and)(rvecm, mx1);                                                           \
      r2m[m]        = (and)(rvecm, mx2);                                                           \
                                                                                                   \
      r0s[m] = (shift_left)(tmp1, 2);                                                              \
      r1s[m] = (shift_left)(tmp2, 1);                                                              \
    }                                                                                              \
  } while (0)

#define bitsliced_mm_step_2(sc, type, and, xor, shift_right)                                       \
  do {                                                                                             \
    const type maskm ATTR_ALIGNED(alignof(type)) =                                                 \
        *((const type*)ASSUME_ALIGNED(CONST_FIRST_ROW(mask->mask), alignof(type)));                \
    for (unsigned int m = 0; m < sc; ++m) {                                                        \
      const type inm ATTR_ALIGNED(alignof(type)) =                                                 \
          *((const type*)ASSUME_ALIGNED(CONST_FIRST_ROW(in[m]), alignof(type)));                   \
      type* outm = (type*)ASSUME_ALIGNED(CONST_FIRST_ROW(out[m]), alignof(type));                  \
                                                                                                   \
      type tmp1 = (xor)(r2m[m], x0s[m]);                                                           \
      type tmp2 = (xor)(x0s[m], x1s[m]);                                                           \
      type tmp3 = (xor)(tmp2, r1m[m]);                                                             \
                                                                                                   \
      type mout = (and)(maskm, inm);                                                               \
                                                                                                   \
      type tmp4 = (xor)(tmp2, r0m[m]);                                                             \
      tmp4      = (xor)(tmp4, x2m[m]);                                                             \
      mout      = (xor)(mout, tmp4);                                                               \
                                                                                                   \
      tmp2 = (shift_right)(tmp1, 2);                                                               \
      mout = (xor)(mout, tmp2);                                                                    \
                                                                                                   \
      tmp1  = (shift_right)(tmp3, 1);                                                              \
      *outm = (xor)(mout, tmp1);                                                                   \
    }                                                                                              \
  } while (0)

#define bitsliced_mm_step_1_multiple_of_128(sc, type, and, shift_left, size)                       \
  type r0m[sc][size] ATTR_ALIGNED(alignof(type));                                                  \
  type r0s[sc][size] ATTR_ALIGNED(alignof(type));                                                  \
  type r1m[sc][size] ATTR_ALIGNED(alignof(type));                                                  \
  type r1s[sc][size] ATTR_ALIGNED(alignof(type));                                                  \
  type r2m[sc][size] ATTR_ALIGNED(alignof(type));                                                  \
  type x0s[sc][size] ATTR_ALIGNED(alignof(type));                                                  \
  type x1s[sc][size] ATTR_ALIGNED(alignof(type));                                                  \
  type x2m[sc][size] ATTR_ALIGNED(alignof(type));                                                  \
  const type* mx2 ATTR_ALIGNED(alignof(type)) =                                                    \
      ((const type*)ASSUME_ALIGNED(CONST_FIRST_ROW(mask->x2), alignof(type)));                     \
  do {                                                                                             \
    const type* mx0 ATTR_ALIGNED(alignof(type)) =                                                  \
        ((const type*)ASSUME_ALIGNED(CONST_FIRST_ROW(mask->x0), alignof(type)));                   \
    const type* mx1 ATTR_ALIGNED(alignof(type)) =                                                  \
        ((const type*)ASSUME_ALIGNED(CONST_FIRST_ROW(mask->x1), alignof(type)));                   \
                                                                                                   \
    for (unsigned int m = 0; m < (sc); ++m) {                                                      \
      const type* inm ATTR_ALIGNED(alignof(type)) =                                                \
          ((const type*)ASSUME_ALIGNED(CONST_FIRST_ROW(in[m]), alignof(type)));                    \
      const type* rvecm ATTR_ALIGNED(alignof(type)) =                                              \
          ((const type*)ASSUME_ALIGNED(CONST_FIRST_ROW(rvec[m]), alignof(type)));                  \
                                                                                                   \
      type tmp1[size] ATTR_ALIGNED(alignof(type));                                                 \
      type tmp2[size] ATTR_ALIGNED(alignof(type));                                                 \
      (and)(tmp1, inm, mx0);                                                                       \
      (and)(tmp2, inm, mx1);                                                                       \
      (and)(x2m[m], inm, mx2);                                                                     \
                                                                                                   \
      (shift_left)(x0s[m], tmp1, 2);                                                               \
      (shift_left)(x1s[m], tmp2, 1);                                                               \
                                                                                                   \
      (and)(tmp1, rvecm, mx0);                                                                     \
      memcpy(r0m[m], tmp1, size * sizeof(type));                                                   \
                                                                                                   \
      (and)(tmp2, rvecm, mx1);                                                                     \
      memcpy(r1m[m], tmp2, size * sizeof(type));                                                   \
                                                                                                   \
      (and)(r2m[m], rvecm, mx2);                                                                   \
                                                                                                   \
      (shift_left)(r0s[m], tmp1, 2);                                                               \
      (shift_left)(r1s[m], tmp2, 1);                                                               \
    }                                                                                              \
  } while (0)

#define bitsliced_mm_step_2_multiple_of_128(sc, type, and, xor, shift_right, size)                 \
  do {                                                                                             \
    const type* maskm ATTR_ALIGNED(alignof(type)) =                                                \
        ((const type*)ASSUME_ALIGNED(CONST_FIRST_ROW(mask->mask), alignof(type)));                 \
    for (unsigned int m = 0; m < sc; ++m) {                                                        \
      const type* inm ATTR_ALIGNED(alignof(type)) =                                                \
          ((const type*)ASSUME_ALIGNED(CONST_FIRST_ROW(in[m]), alignof(type)));                    \
      type* outm = (type*)ASSUME_ALIGNED(CONST_FIRST_ROW(out[m]), alignof(type));                  \
                                                                                                   \
      type tmp1[size], tmp2[size], tmp3[size];                                                     \
      (xor)(tmp1, r2m[m], x0s[m]);                                                                 \
      (xor)(tmp2, x0s[m], x1s[m]);                                                                 \
      (xor)(tmp3, tmp2, r1m[m]);                                                                   \
                                                                                                   \
      type mout[size];                                                                             \
      (and)(mout, maskm, inm);                                                                     \
                                                                                                   \
      type tmp4[size];                                                                             \
      (xor)(tmp4, tmp2, r0m[m]);                                                                   \
      (xor)(tmp4, tmp4, x2m[m]);                                                                   \
      (xor)(mout, mout, tmp4);                                                                     \
                                                                                                   \
      (shift_right)(tmp2, tmp1, 2);                                                                \
      (xor)(mout, mout, tmp2);                                                                     \
      (shift_right)(tmp1, tmp3, 1);                                                                \
      (xor)(outm, mout, tmp1);                                                                     \
    }                                                                                              \
  } while (0)

#if defined(WITH_SSE2)
ATTR_TARGET("sse2")
static void mpc_sbox_layer_bitsliced_128_sse(mzd_local_t** out, mzd_local_t* const* in,
                                             view_t* view, mzd_local_t** rvec, mask_t const* mask) {
  bitsliced_mm_step_1(SC_PROOF, __m128i, _mm_and_si128, mm128_shift_left);

  mpc_clear(view->s, SC_PROOF);
  mpc_and_sse(r0m, x0s, x1s, r2m, view, 0);
  mpc_and_sse(r2m, x1s, x2m, r1s, view, 1);
  mpc_and_sse(r1m, x0s, x2m, r0s, view, 2);

  bitsliced_mm_step_2(SC_PROOF, __m128i, _mm_and_si128, _mm_xor_si128, mm128_shift_right);
}

ATTR_TARGET("sse2")
static void mpc_sbox_layer_bitsliced_verify_128_sse(mzd_local_t** out, mzd_local_t* const* in,
                                                    view_t* view, mzd_local_t** rvec,
                                                    mask_t const* mask) {
  bitsliced_mm_step_1(SC_VERIFY, __m128i, _mm_and_si128, mm128_shift_left);

  mzd_local_clear(view->s[0]);
  mpc_and_verify_sse(r0m, x0s, x1s, r2m, view, mx2, 0);
  mpc_and_verify_sse(r2m, x1s, x2m, r1s, view, mx2, 1);
  mpc_and_verify_sse(r1m, x0s, x2m, r0s, view, mx2, 2);

  bitsliced_mm_step_2(SC_VERIFY, __m128i, _mm_and_si128, _mm_xor_si128, mm128_shift_right);
}

ATTR_TARGET("sse2")
static void mpc_sbox_layer_bitsliced_256_sse(mzd_local_t** out, mzd_local_t* const* in,
                                             view_t* view, mzd_local_t** rvec, mask_t const* mask) {
  bitsliced_mm_step_1_multiple_of_128(SC_PROOF, __m128i, mm256_and_sse, mm256_shift_left_sse, 2);

  mpc_clear(view->s, SC_PROOF);
  mpc_and_256_sse(r0m, x0s, x1s, r2m, view, 0);
  mpc_and_256_sse(r2m, x1s, x2m, r1s, view, 1);
  mpc_and_256_sse(r1m, x0s, x2m, r0s, view, 2);

  bitsliced_mm_step_2_multiple_of_128(SC_PROOF, __m128i, mm256_and_sse, mm256_xor_sse,
                                      mm256_shift_right_sse, 2);
}

ATTR_TARGET("sse2")
static void mpc_sbox_layer_bitsliced_verify_256_sse(mzd_local_t** out, mzd_local_t* const* in,
                                                    view_t* view, mzd_local_t** rvec,
                                                    mask_t const* mask) {
  bitsliced_mm_step_1_multiple_of_128(SC_VERIFY, __m128i, mm256_and_sse, mm256_shift_left_sse, 2);

  mzd_local_clear(view->s[0]);
  mpc_and_verify_256_sse(r0m, x0s, x1s, r2m, view, mx2, 0);
  mpc_and_verify_256_sse(r2m, x1s, x2m, r1s, view, mx2, 1);
  mpc_and_verify_256_sse(r1m, x0s, x2m, r0s, view, mx2, 2);

  bitsliced_mm_step_2_multiple_of_128(SC_VERIFY, __m128i, mm256_and_sse, mm256_xor_sse,
                                      mm256_shift_right_sse, 2);
}

ATTR_TARGET("sse2")
static void mpc_sbox_layer_bitsliced_384_sse(mzd_local_t** out, mzd_local_t* const* in,
                                             view_t* view, mzd_local_t** rvec, mask_t const* mask) {
  bitsliced_mm_step_1_multiple_of_128(SC_PROOF, __m128i, mm384_and_sse, mm384_shift_left_sse, 3);

  mpc_clear(view->s, SC_PROOF);
  mpc_and_384_sse(r0m, x0s, x1s, r2m, view, 0);
  mpc_and_384_sse(r2m, x1s, x2m, r1s, view, 1);
  mpc_and_384_sse(r1m, x0s, x2m, r0s, view, 2);

  bitsliced_mm_step_2_multiple_of_128(SC_PROOF, __m128i, mm384_and_sse, mm384_xor_sse,
                                      mm384_shift_right_sse, 3);
}

ATTR_TARGET("sse2")
static void mpc_sbox_layer_bitsliced_verify_384_sse(mzd_local_t** out, mzd_local_t* const* in,
                                                    view_t* view, mzd_local_t** rvec,
                                                    mask_t const* mask) {
  bitsliced_mm_step_1_multiple_of_128(SC_VERIFY, __m128i, mm384_and_sse, mm384_shift_left_sse, 3);

  mzd_local_clear(view->s[0]);
  mpc_and_verify_384_sse(r0m, x0s, x1s, r2m, view, mx2, 0);
  mpc_and_verify_384_sse(r2m, x1s, x2m, r1s, view, mx2, 1);
  mpc_and_verify_384_sse(r1m, x0s, x2m, r0s, view, mx2, 2);

  bitsliced_mm_step_2_multiple_of_128(SC_VERIFY, __m128i, mm384_and_sse, mm384_xor_sse,
                                      mm384_shift_right_sse, 3);
}

ATTR_TARGET("sse2")
static void mpc_sbox_layer_bitsliced_512_sse(mzd_local_t** out, mzd_local_t* const* in,
                                             view_t* view, mzd_local_t** rvec, mask_t const* mask) {
  bitsliced_mm_step_1_multiple_of_128(SC_PROOF, __m128i, mm512_and_sse, mm512_shift_left_sse, 4);

  mpc_clear(view->s, SC_PROOF);
  mpc_and_512_sse(r0m, x0s, x1s, r2m, view, 0);
  mpc_and_512_sse(r2m, x1s, x2m, r1s, view, 1);
  mpc_and_512_sse(r1m, x0s, x2m, r0s, view, 2);

  bitsliced_mm_step_2_multiple_of_128(SC_PROOF, __m128i, mm512_and_sse, mm512_xor_sse,
                                      mm512_shift_right_sse, 4);
}

ATTR_TARGET("sse2")
static void mpc_sbox_layer_bitsliced_verify_512_sse(mzd_local_t** out, mzd_local_t* const* in,
                                                    view_t* view, mzd_local_t** rvec,
                                                    mask_t const* mask) {
  bitsliced_mm_step_1_multiple_of_128(SC_VERIFY, __m128i, mm512_and_sse, mm512_shift_left_sse, 4);

  mzd_local_clear(view->s[0]);
  mpc_and_verify_512_sse(r0m, x0s, x1s, r2m, view, mx2, 0);
  mpc_and_verify_512_sse(r2m, x1s, x2m, r1s, view, mx2, 1);
  mpc_and_verify_512_sse(r1m, x0s, x2m, r0s, view, mx2, 2);

  bitsliced_mm_step_2_multiple_of_128(SC_VERIFY, __m128i, mm512_and_sse, mm512_xor_sse,
                                      mm512_shift_right_sse, 4);
}
#endif

#if defined(WITH_AVX2)
ATTR_TARGET("avx2")
static void mpc_sbox_layer_bitsliced_256_avx(mzd_local_t** out, mzd_local_t* const* in,
                                             view_t* view, mzd_local_t** rvec, mask_t const* mask) {
  bitsliced_mm_step_1(SC_PROOF, __m256i, _mm256_and_si256, mm256_shift_left);

  mpc_clear(view->s, SC_PROOF);
  mpc_and_avx(r0m, x0s, x1s, r2m, view, 0);
  mpc_and_avx(r2m, x1s, x2m, r1s, view, 1);
  mpc_and_avx(r1m, x0s, x2m, r0s, view, 2);

  bitsliced_mm_step_2(SC_PROOF, __m256i, _mm256_and_si256, _mm256_xor_si256, mm256_shift_right);
}

ATTR_TARGET("avx2")
static void mpc_sbox_layer_bitsliced_verify_256_avx(mzd_local_t** out, mzd_local_t** in,
                                                    view_t* view, mzd_local_t* const* rvec,
                                                    mask_t const* mask) {
  bitsliced_mm_step_1(SC_VERIFY, __m256i, _mm256_and_si256, mm256_shift_left);

  mzd_local_clear(view->s[0]);
  mpc_and_verify_avx(r0m, x0s, x1s, r2m, view, mx2, 0);
  mpc_and_verify_avx(r2m, x1s, x2m, r1s, view, mx2, 1);
  mpc_and_verify_avx(r1m, x0s, x2m, r0s, view, mx2, 2);

  bitsliced_mm_step_2(SC_VERIFY, __m256i, _mm256_and_si256, _mm256_xor_si256, mm256_shift_right);
}

ATTR_TARGET("avx2")
static void mpc_sbox_layer_bitsliced_512_avx(mzd_local_t** out, mzd_local_t* const* in,
                                             view_t* view, mzd_local_t** rvec, mask_t const* mask) {
  bitsliced_mm_step_1_multiple_of_128(SC_PROOF, __m256i, mm512_and_avx, mm512_shift_left_avx, 2);

  mpc_clear(view->s, SC_PROOF);
  mpc_and_512_avx(r0m, x0s, x1s, r2m, view, 0);
  mpc_and_512_avx(r2m, x1s, x2m, r1s, view, 1);
  mpc_and_512_avx(r1m, x0s, x2m, r0s, view, 2);

  bitsliced_mm_step_2_multiple_of_128(SC_PROOF, __m256i, mm512_and_avx, mm512_xor_avx,
                                      mm512_shift_right_avx, 2);
}

ATTR_TARGET("avx2")
static void mpc_sbox_layer_bitsliced_verify_512_avx(mzd_local_t** out, mzd_local_t** in,
                                                    view_t* view, mzd_local_t* const* rvec,
                                                    mask_t const* mask) {
  bitsliced_mm_step_1_multiple_of_128(SC_VERIFY, __m256i, mm512_and_avx, mm512_shift_left_avx, 2);

  mzd_local_clear(view->s[0]);
  mpc_and_verify_512_avx(r0m, x0s, x1s, r2m, view, mx2, 0);
  mpc_and_verify_512_avx(r2m, x1s, x2m, r1s, view, mx2, 1);
  mpc_and_verify_512_avx(r1m, x0s, x2m, r0s, view, mx2, 2);

  bitsliced_mm_step_2_multiple_of_128(SC_VERIFY, __m256i, mm512_and_avx, mm512_xor_avx,
                                      mm512_shift_right_avx, 2);
}
#endif

#if defined(WITH_NEON)
static void mpc_sbox_layer_bitsliced_128_neon(mzd_local_t** out, mzd_local_t* const* in,
                                              view_t* view, mzd_local_t** rvec,
                                              mask_t const* mask) {
  bitsliced_mm_step_1(SC_PROOF, uint32x4_t, vandq_u32, mm128_shift_left);

  mpc_clear(view->s, SC_PROOF);
  mpc_and_neon(r0m, x0s, x1s, r2m, view, 0);
  mpc_and_neon(r2m, x1s, x2m, r1s, view, 1);
  mpc_and_neon(r1m, x0s, x2m, r0s, view, 2);

  bitsliced_mm_step_2(SC_PROOF, uint32x4_t, vandq_u32, veorq_u32, mm128_shift_right);
}

static void mpc_sbox_layer_bitsliced_verify_128_neon(mzd_local_t** out, mzd_local_t* const* in,
                                                     view_t* view, mzd_local_t** rvec,
                                                     mask_t const* mask) {
  bitsliced_mm_step_1(SC_VERIFY, uint32x4_t, vandq_u32, mm128_shift_left);

  mzd_local_clear(view->s[0]);
  mpc_and_verify_neon(r0m, x0s, x1s, r2m, view, mx2, 0);
  mpc_and_verify_neon(r2m, x1s, x2m, r1s, view, mx2, 1);
  mpc_and_verify_neon(r1m, x0s, x2m, r0s, view, mx2, 2);

  bitsliced_mm_step_2(SC_VERIFY, uint32x4_t, vandq_u32, veorq_u32, mm128_shift_right);
}

static void mpc_sbox_layer_bitsliced_256_neon(mzd_local_t** out, mzd_local_t* const* in,
                                              view_t* view, mzd_local_t** rvec,
                                              mask_t const* mask) {
  bitsliced_mm_step_1_multiple_of_128(SC_PROOF, uint32x4_t, mm256_and, mm256_shift_left, 2);

  mpc_clear(view->s, SC_PROOF);
  mpc_and_256_neon(r0m, x0s, x1s, r2m, view, 0);
  mpc_and_256_neon(r2m, x1s, x2m, r1s, view, 1);
  mpc_and_256_neon(r1m, x0s, x2m, r0s, view, 2);

  bitsliced_mm_step_2_multiple_of_128(SC_PROOF, uint32x4_t, mm256_and, mm256_xor, mm256_shift_right,
                                      2);
}

static void mpc_sbox_layer_bitsliced_verify_256_neon(mzd_local_t** out, mzd_local_t* const* in,
                                                     view_t* view, mzd_local_t** rvec,
                                                     mask_t const* mask) {
  bitsliced_mm_step_1_multiple_of_128(SC_VERIFY, uint32x4_t, mm256_and, mm256_shift_left, 2);

  mzd_local_clear(view->s[0]);
  mpc_and_verify_256_neon(r0m, x0s, x1s, r2m, view, mx2, 0);
  mpc_and_verify_256_neon(r2m, x1s, x2m, r1s, view, mx2, 1);
  mpc_and_verify_256_neon(r1m, x0s, x2m, r0s, view, mx2, 2);

  bitsliced_mm_step_2_multiple_of_128(SC_VERIFY, uint32x4_t, mm256_and, mm256_xor,
                                      mm256_shift_right, 2);
}

static void mpc_sbox_layer_bitsliced_384_neon(mzd_local_t** out, mzd_local_t* const* in,
                                              view_t* view, mzd_local_t** rvec,
                                              mask_t const* mask) {
  bitsliced_mm_step_1_multiple_of_128(SC_PROOF, uint32x4_t, mm384_and, mm384_shift_left, 3);

  mpc_clear(view->s, SC_PROOF);
  mpc_and_384_neon(r0m, x0s, x1s, r2m, view, 0);
  mpc_and_384_neon(r2m, x1s, x2m, r1s, view, 1);
  mpc_and_384_neon(r1m, x0s, x2m, r0s, view, 2);

  bitsliced_mm_step_2_multiple_of_128(SC_PROOF, uint32x4_t, mm384_and, mm384_xor, mm384_shift_right,
                                      3);
}

static void mpc_sbox_layer_bitsliced_verify_384_neon(mzd_local_t** out, mzd_local_t* const* in,
                                                     view_t* view, mzd_local_t** rvec,
                                                     mask_t const* mask) {
  bitsliced_mm_step_1_multiple_of_128(SC_VERIFY, uint32x4_t, mm384_and, mm384_shift_left, 3);

  mzd_local_clear(view->s[0]);
  mpc_and_verify_384_neon(r0m, x0s, x1s, r2m, view, mx2, 0);
  mpc_and_verify_384_neon(r2m, x1s, x2m, r1s, view, mx2, 1);
  mpc_and_verify_384_neon(r1m, x0s, x2m, r0s, view, mx2, 2);

  bitsliced_mm_step_2_multiple_of_128(SC_VERIFY, uint32x4_t, mm384_and, mm384_xor,
                                      mm384_shift_right, 3);
}

static void mpc_sbox_layer_bitsliced_512_neon(mzd_local_t** out, mzd_local_t* const* in,
                                              view_t* view, mzd_local_t** rvec,
                                              mask_t const* mask) {
  bitsliced_mm_step_1_multiple_of_128(SC_PROOF, uint32x4_t, mm512_and, mm512_shift_left, 4);

  mpc_clear(view->s, SC_PROOF);
  mpc_and_512_neon(r0m, x0s, x1s, r2m, view, 0);
  mpc_and_512_neon(r2m, x1s, x2m, r1s, view, 1);
  mpc_and_512_neon(r1m, x0s, x2m, r0s, view, 2);

  bitsliced_mm_step_2_multiple_of_128(SC_PROOF, uint32x4_t, mm512_and, mm512_xor, mm512_shift_right,
                                      4);
}

static void mpc_sbox_layer_bitsliced_verify_512_neon(mzd_local_t** out, mzd_local_t* const* in,
                                                     view_t* view, mzd_local_t** rvec,
                                                     mask_t const* mask) {
  bitsliced_mm_step_1_multiple_of_128(SC_VERIFY, uint32x4_t, mm512_and, mm512_shift_left, 4);

  mzd_local_clear(view->s[0]);
  mpc_and_verify_512_neon(r0m, x0s, x1s, r2m, view, mx2, 0);
  mpc_and_verify_512_neon(r2m, x1s, x2m, r1s, view, mx2, 1);
  mpc_and_verify_512_neon(r1m, x0s, x2m, r0s, view, mx2, 2);

  bitsliced_mm_step_2_multiple_of_128(SC_VERIFY, uint32x4_t, mm512_and, mm512_xor,
                                      mm512_shift_right, 4);
}
#endif
#endif

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

// 384 bit using SSE2
#undef SIGN_SBOX
#undef VERIFY_SBOX
#define SIGN_SBOX mpc_sbox_layer_bitsliced_384_sse
#define VERIFY_SBOX mpc_sbox_layer_bitsliced_verify_384_sse

#define SIGN mpc_lowmc_call_384_sse
#define VERIFY mpc_lowmc_call_verify_384_sse
#include "mpc_lowmc.c.i"

// 512 bit using SSE2
#undef SIGN_SBOX
#undef VERIFY_SBOX
#define SIGN_SBOX mpc_sbox_layer_bitsliced_512_sse
#define VERIFY_SBOX mpc_sbox_layer_bitsliced_verify_512_sse

#define SIGN mpc_lowmc_call_512_sse
#define VERIFY mpc_lowmc_call_verify_512_sse
#include "mpc_lowmc.c.i"
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

// 384 bit using AVX2
#undef SIGN_SBOX
#undef VERIFY_SBOX
#define SIGN_SBOX mpc_sbox_layer_bitsliced_512_avx
#define VERIFY_SBOX mpc_sbox_layer_bitsliced_verify_512_avx

#define SIGN mpc_lowmc_call_384_avx
#define VERIFY mpc_lowmc_call_verify_384_avx
#include "mpc_lowmc.c.i"

// 512 bit using AVX2
#undef SIGN_SBOX
#undef VERIFY_SBOX
#define SIGN_SBOX mpc_sbox_layer_bitsliced_512_avx
#define VERIFY_SBOX mpc_sbox_layer_bitsliced_verify_512_avx

#define SIGN mpc_lowmc_call_512_avx
#define VERIFY mpc_lowmc_call_verify_512_avx
#include "mpc_lowmc.c.i"
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
#if defined(WITH_LOWMC_192_192_30)
#define LOWMC_INSTANCE_10 (&lowmc_192_192_30)
#endif
#if defined(WITH_LOWMC_192_192_284)
#define LOWMC_INSTANCE_1 (&lowmc_192_192_284)
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

// 384 bit using NEON
#undef SIGN_SBOX
#undef VERIFY_SBOX
#define SIGN_SBOX mpc_sbox_layer_bitsliced_384_neon
#define VERIFY_SBOX mpc_sbox_layer_bitsliced_verify_384_neon

#define SIGN mpc_lowmc_call_384_neon
#define VERIFY mpc_lowmc_call_verify_384_neon
#include "mpc_lowmc.c.i"

// 512 bit using NEON
#undef SIGN_SBOX
#undef VERIFY_SBOX
#define SIGN_SBOX mpc_sbox_layer_bitsliced_512_neon
#define VERIFY_SBOX mpc_sbox_layer_bitsliced_verify_512_neon

#define SIGN mpc_lowmc_call_512_neon
#define VERIFY mpc_lowmc_call_verify_512_neon
#include "mpc_lowmc.c.i"
#endif
#endif
#endif

zkbpp_lowmc_implementation_f get_zkbpp_lowmc_implementation(const lowmc_t* lowmc) {
#if defined(WITH_OPT)
#if defined(WITH_AVX2)
  if (CPU_SUPPORTS_AVX2) {
    if (lowmc->m == 10) {
      switch (lowmc->n) {
        case 128:
          return general_or_10(lowmc, mpc_lowmc_call_128_avx);
        case 192:
          return general_or_10(lowmc, mpc_lowmc_call_192_avx);
        case 256:
          return general_or_10(lowmc, mpc_lowmc_call_256_avx);
#if defined(WITH_CUSTOM_INSTANCES)
        case 384:
          return general_or_10(lowmc, mpc_lowmc_call_384_avx);
        case 512:
          return general_or_10(lowmc, mpc_lowmc_call_512_avx);
#endif
      }
    }
    if (lowmc->m == 1) {
      switch (lowmc->n) {
        case 128:
          return general_or_1(lowmc, mpc_lowmc_call_128_avx);
        case 192:
          return general_or_1(lowmc, mpc_lowmc_call_192_avx);
        case 256:
          return general_or_1(lowmc, mpc_lowmc_call_256_avx);
#if defined(WITH_CUSTOM_INSTANCES)
        case 384:
          return general_or_1(lowmc, mpc_lowmc_call_384_avx);
        case 512:
          return general_or_1(lowmc, mpc_lowmc_call_512_avx);
#endif
      }
    }
  }
#endif
#if defined(WITH_SSE2)
  if (CPU_SUPPORTS_SSE2) {
    if(lowmc->m == 10) {
      switch (lowmc->n) {
        case 128:
          return general_or_10(lowmc, mpc_lowmc_call_128_sse);
        case 192:
          return general_or_10(lowmc, mpc_lowmc_call_192_sse);
        case 256:
          return general_or_10(lowmc, mpc_lowmc_call_256_sse);
#if defined(WITH_CUSTOM_INSTANCES)
        case 384:
          return general_or_10(lowmc, mpc_lowmc_call_384_sse);
        case 512:
          return general_or_10(lowmc, mpc_lowmc_call_512_sse);
#endif
      }
    }
    if(lowmc->m == 1) {
      switch (lowmc->n) {
        case 128:
          return general_or_1(lowmc, mpc_lowmc_call_128_sse);
        case 192:
          return general_or_1(lowmc, mpc_lowmc_call_192_sse);
        case 256:
          return general_or_1(lowmc, mpc_lowmc_call_256_sse);
#if defined(WITH_CUSTOM_INSTANCES)
        case 384:
          return general_or_1(lowmc, mpc_lowmc_call_384_sse);
        case 512:
          return general_or_1(lowmc, mpc_lowmc_call_512_sse);
#endif
      }
    }
  }
#endif
#if defined(WITH_NEON)
  if (CPU_SUPPORTS_NEON)
    switch (lowmc->n) {
    case 128:
      return general_or_10(lowmc, mpc_lowmc_call_128_neon);
    case 192:
      return general_or_10(lowmc, mpc_lowmc_call_192_neon);
    case 256:
      return general_or_10(lowmc, mpc_lowmc_call_256_neon);
#if defined(WITH_CUSTOM_INSTANCES)
    case 384:
      return general_or_10(lowmc, mpc_lowmc_call_384_neon);
    case 512:
      return general_or_10(lowmc, mpc_lowmc_call_512_neon);
#endif
    }
#endif
#endif

  (void)lowmc;
  if(lowmc->m == 10)
    return general_or_10(lowmc, mpc_lowmc_call);
  if(lowmc->m == 1)
    return general_or_1(lowmc, mpc_lowmc_call);
  return NULL;
}

zkbpp_lowmc_verify_implementation_f get_zkbpp_lowmc_verify_implementation(const lowmc_t* lowmc) {
#if defined(WITH_OPT)
#if defined(WITH_AVX2)
  if (CPU_SUPPORTS_AVX2) {
    if(lowmc->m == 10) {
      switch (lowmc->n) {
        case 128:
          return general_or_10(lowmc, mpc_lowmc_call_verify_128_avx);
        case 192:
          return general_or_10(lowmc, mpc_lowmc_call_verify_192_avx);
        case 256:
          return general_or_10(lowmc, mpc_lowmc_call_verify_256_avx);
#if defined(WITH_CUSTOM_INSTANCES)
        case 384:
          return general_or_10(lowmc, mpc_lowmc_call_verify_384_avx);
        case 512:
          return general_or_10(lowmc, mpc_lowmc_call_verify_512_avx);
#endif
      }
    }
    if(lowmc->m == 1) {
      switch (lowmc->n) {
        case 128:
          return general_or_1(lowmc, mpc_lowmc_call_verify_128_avx);
        case 192:
          return general_or_1(lowmc, mpc_lowmc_call_verify_192_avx);
        case 256:
          return general_or_1(lowmc, mpc_lowmc_call_verify_256_avx);
#if defined(WITH_CUSTOM_INSTANCES)
        case 384:
          return general_or_1(lowmc, mpc_lowmc_call_verify_384_avx);
        case 512:
          return general_or_1(lowmc, mpc_lowmc_call_verify_512_avx);
#endif
      }
    }
  }
#endif
#if defined(WITH_SSE2)
  if (CPU_SUPPORTS_SSE2) {
    if(lowmc->m == 10) {
      switch (lowmc->n) {
        case 128:
          return general_or_10(lowmc, mpc_lowmc_call_verify_128_sse);
        case 192:
          return general_or_10(lowmc, mpc_lowmc_call_verify_192_sse);
        case 256:
          return general_or_10(lowmc, mpc_lowmc_call_verify_256_sse);
#if defined(WITH_CUSTOM_INSTANCES)
        case 384:
          return general_or_10(lowmc, mpc_lowmc_call_verify_384_sse);
        case 512:
          return general_or_10(lowmc, mpc_lowmc_call_verify_512_sse);
#endif
      }
    }
    if(lowmc->m == 1) {
      switch (lowmc->n) {
        case 128:
          return general_or_1(lowmc, mpc_lowmc_call_verify_128_sse);
        case 192:
          return general_or_1(lowmc, mpc_lowmc_call_verify_192_sse);
        case 256:
          return general_or_1(lowmc, mpc_lowmc_call_verify_256_sse);
#if defined(WITH_CUSTOM_INSTANCES)
        case 384:
          return general_or_1(lowmc, mpc_lowmc_call_verify_384_sse);
        case 512:
          return general_or_1(lowmc, mpc_lowmc_call_verify_512_sse);
#endif
      }
    }
  }
#endif
#if defined(WITH_NEON)
  if (CPU_SUPPORTS_NEON) {
    if(lowmc->m == 10) {
      switch (lowmc->n) {
      case 128:
        return general_or_10(lowmc, mpc_lowmc_call_verify_128_neon);
      case 192:
        return general_or_10(lowmc, mpc_lowmc_call_verify_192_neon);
      case 256:
        return general_or_10(lowmc, mpc_lowmc_call_verify_256_neon);
#if defined(WITH_CUSTOM_INSTANCES)
      case 384:
        return general_or_10(lowmc, mpc_lowmc_call_verify_384_neon);
      case 512:
        return general_or_10(lowmc, mpc_lowmc_call_verify_512_neon);
#endif
      }
    }
    if(lowmc->m == 1) {
      switch (lowmc->n) {
      case 128:
        return general_or_1(lowmc, mpc_lowmc_call_verify_128_neon);
      case 192:
        return general_or_1(lowmc, mpc_lowmc_call_verify_192_neon);
      case 256:
        return general_or_1(lowmc, mpc_lowmc_call_verify_256_neon);
#if defined(WITH_CUSTOM_INSTANCES)
      case 384:
        return general_or_1(lowmc, mpc_lowmc_call_verify_384_neon);
      case 512:
        return general_or_1(lowmc, mpc_lowmc_call_verify_512_neon);
#endif
      }
    }
  }
#endif
#endif

  (void)lowmc;
  if(lowmc->m == 10)
    return general_or_10(lowmc, mpc_lowmc_call_verify);
  if(lowmc->m == 1)
    return general_or_1(lowmc, mpc_lowmc_call_verify);
  return NULL;
}
