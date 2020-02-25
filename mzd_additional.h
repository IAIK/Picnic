/*
 *  This file is part of the optimized implementation of the Picnic signature scheme.
 *  See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */

/* Inspired by m4ri's mzd implementation, but completely re-written for our use-case. */

#ifndef MZD_ADDITIONAL_H
#define MZD_ADDITIONAL_H

#include "macros.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

typedef uint64_t word;
#define WORD_C(v) UINT64_C(v)

#if defined(WITH_OPT)
#include "simd.h"
#endif

typedef union {
  word w64[4];
#if defined(WITH_OPT)
#if defined(WITH_SSE2) || defined(WITH_NEON)
  word128 w128[2];
#endif
#if defined(WITH_AVX2)
  word256 w256;
#endif
#endif
} block_t ATTR_ALIGNED(32);

/**
 * Representation of matrices and vectors
 *
 * The basic memory unit is a block of 256 bit. Each row is stored in (possible multiple) blocks
 * depending on the number of columns. Matrices with up to 128 columns are the only excpetion. In
 * this case a block actually contains two rows. The row with even index is contained in w64[0] and
 * w61[1], the row with odd index is contained in w64[2] and w64[3].
 */
typedef block_t mzd_local_t;

mzd_local_t* mzd_local_init_ex(uint32_t r, uint32_t c, bool clear) ATTR_ASSUME_ALIGNED(32);

#define mzd_local_init(r, c) mzd_local_init_ex(r, c, true)

void mzd_local_free(mzd_local_t* v);

void mzd_local_init_multiple_ex(mzd_local_t** dst, size_t n, uint32_t r, uint32_t c, bool clear)
    ATTR_NONNULL_ARG(1);

#define mzd_local_init_multiple(dst, n, r, c) mzd_local_init_multiple_ex(dst, n, r, c, true)

/**
 * mzd_local_free for mzd_local_init_multiple.
 */
void mzd_local_free_multiple(mzd_local_t** vs);

void mzd_copy_uint64_128(mzd_local_t* dst, mzd_local_t const* src) ATTR_NONNULL;
void mzd_copy_uint64_192(mzd_local_t* dst, mzd_local_t const* src) ATTR_NONNULL;
void mzd_copy_uint64_256(mzd_local_t* dst, mzd_local_t const* src) ATTR_NONNULL;
void mzd_copy_s128_128(mzd_local_t* dst, mzd_local_t const* src) ATTR_NONNULL;
void mzd_copy_s128_256(mzd_local_t* dst, mzd_local_t const* src) ATTR_NONNULL;
void mzd_copy_s256_128(mzd_local_t* dst, mzd_local_t const* src) ATTR_NONNULL;
void mzd_copy_s256_256(mzd_local_t* dst, mzd_local_t const* src) ATTR_NONNULL;

/**
 * mzd_xor variants
 */
void mzd_xor_uint64_128(mzd_local_t* res, mzd_local_t const* first,
                        mzd_local_t const* second) ATTR_NONNULL;
void mzd_xor_uint64_192(mzd_local_t* res, mzd_local_t const* first,
                        mzd_local_t const* second) ATTR_NONNULL;
void mzd_xor_uint64_256(mzd_local_t* res, mzd_local_t const* first,
                        mzd_local_t const* second) ATTR_NONNULL;
void mzd_xor_s128_128(mzd_local_t* res, mzd_local_t const* first,
                      mzd_local_t const* second) ATTR_NONNULL;
void mzd_xor_s128_256(mzd_local_t* res, mzd_local_t const* first,
                      mzd_local_t const* second) ATTR_NONNULL;
void mzd_xor_s256_128(mzd_local_t* res, mzd_local_t const* first,
                      mzd_local_t const* second) ATTR_NONNULL;
void mzd_xor_s256_256(mzd_local_t* res, mzd_local_t const* first,
                      mzd_local_t const* second) ATTR_NONNULL;

/**
 * mzd_and variants
 */
void mzd_and_uint64_128(mzd_local_t* res, mzd_local_t const* first,
                        mzd_local_t const* second) ATTR_NONNULL;
void mzd_and_uint64_192(mzd_local_t* res, mzd_local_t const* first,
                        mzd_local_t const* second) ATTR_NONNULL;
void mzd_and_uint64_256(mzd_local_t* res, mzd_local_t const* first,
                        mzd_local_t const* second) ATTR_NONNULL;
void mzd_and_s128_128(mzd_local_t* res, mzd_local_t const* first,
                      mzd_local_t const* second) ATTR_NONNULL;
void mzd_and_s128_256(mzd_local_t* res, mzd_local_t const* first,
                      mzd_local_t const* second) ATTR_NONNULL;
void mzd_and_s256_128(mzd_local_t* res, mzd_local_t const* first,
                      mzd_local_t const* second) ATTR_NONNULL;
void mzd_and_s256_256(mzd_local_t* res, mzd_local_t const* first,
                      mzd_local_t const* second) ATTR_NONNULL;

/**
 * mzd_left_shift and mzd_right_shift variants
 */
void mzd_shift_left_uint64_128(mzd_local_t* val, unsigned int count);
void mzd_shift_right_uint64_128(mzd_local_t* val, unsigned int count);

/**
 * Compute v * A optimized for v being a vector.
 */
void mzd_mul_v_uint64_128(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* At) ATTR_NONNULL;
void mzd_mul_v_uint64_192(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* At) ATTR_NONNULL;
void mzd_mul_v_uint64_256(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* At) ATTR_NONNULL;
void mzd_mul_v_s128_128(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_mul_v_s128_192(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_mul_v_s128_256(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_mul_v_s256_128(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_mul_v_s256_192(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_mul_v_s256_256(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;

/**
 * Compute c + v * A optimized for c and v being vectors.
 */
void mzd_addmul_v_uint64_128(mzd_local_t* c, mzd_local_t const* v,
                             mzd_local_t const* A) ATTR_NONNULL;
void mzd_addmul_v_uint64_192(mzd_local_t* c, mzd_local_t const* v,
                             mzd_local_t const* A) ATTR_NONNULL;
void mzd_addmul_v_uint64_256(mzd_local_t* c, mzd_local_t const* v,
                             mzd_local_t const* A) ATTR_NONNULL;
void mzd_addmul_v_s128_128(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_addmul_v_s128_192(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_addmul_v_s128_256(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_addmul_v_s256_128(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_addmul_v_s256_192(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_addmul_v_s256_256(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;

#define BLOCK(v, b) ((block_t*)ASSUME_ALIGNED(&(v)[(b)], 32))
#define CONST_BLOCK(v, b) ((const block_t*)ASSUME_ALIGNED(&(v)[(b)], 32))

#endif
