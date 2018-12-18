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

typedef struct {
  uint32_t nrows, ncols, width, rowstride;
  uint32_t padding[4];
  uint64_t rows[];
} mzd_local_t ATTR_ALIGNED(32);

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
/**
 * Improved mzd_copy for specific memory layouts. Requires that dst already exists.
 */
void mzd_local_copy(mzd_local_t* dst, mzd_local_t const* src) ATTR_NONNULL_ARG(2);

void mzd_local_clear(mzd_local_t* c) ATTR_NONNULL;

void mzd_xor_uint64(mzd_local_t* res, mzd_local_t const* first, mzd_local_t const* second) ATTR_NONNULL;
void mzd_xor_sse(mzd_local_t* res, mzd_local_t const* first, mzd_local_t const* second) ATTR_NONNULL;
void mzd_xor_sse_128(mzd_local_t* res, mzd_local_t const* first, mzd_local_t const* second) ATTR_NONNULL;
void mzd_xor_sse_256(mzd_local_t* res, mzd_local_t const* first, mzd_local_t const* second) ATTR_NONNULL;
void mzd_xor_avx(mzd_local_t* res, mzd_local_t const* first, mzd_local_t const* second) ATTR_NONNULL;
void mzd_xor_avx_128(mzd_local_t* res, mzd_local_t const* first, mzd_local_t const* second) ATTR_NONNULL;
void mzd_xor_avx_256(mzd_local_t* res, mzd_local_t const* first, mzd_local_t const* second) ATTR_NONNULL;
void mzd_xor_neon(mzd_local_t* res, mzd_local_t const* first, mzd_local_t const* second) ATTR_NONNULL;
void mzd_xor_neon_128(mzd_local_t* res, mzd_local_t const* first, mzd_local_t const* second) ATTR_NONNULL;
void mzd_xor_neon_256(mzd_local_t* res, mzd_local_t const* first, mzd_local_t const* second) ATTR_NONNULL;

/**
 * Compare two vectors for equality. Note that this version is optimized for
 * vectors with a multiple of sizeof(word) * 8 columns.
 *
 * \param first
 *          first vector
 * \param second
 *          second vector
 * \returns true if both vectors are equal, false otherwise.
 */
bool mzd_local_equal(mzd_local_t const* first, mzd_local_t const* second) ATTR_NONNULL;

/**
 * Compute v * A optimized for v being a vector.
 */
void mzd_mul_v_uint64(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* At) ATTR_NONNULL;
void mzd_mul_v_sse(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_mul_v_sse_128(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_mul_v_sse_192(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_mul_v_sse_256(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_mul_v_avx(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_mul_v_avx_128(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_mul_v_avx_192(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_mul_v_avx_256(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_mul_v_neon(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_mul_v_neon_128(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_mul_v_neon_192(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_mul_v_neon_256(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_mul_v_parity_uint64(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* At) ATTR_NONNULL;
void mzd_mul_v_parity_popcnt(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* At) ATTR_NONNULL;

/**
 * Compute v * A optimized for v being a vector, for specific sizes depending on instance
 * Only work for specific sizes and RLL_NEXT algorithm using uint64 operations
 */
void mzd_addmul_v_uint64_30(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_addmul_v_uint64_3(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;

/**
 * Use SSE2
 */
void mzd_addmul_v_sse_30_128(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_addmul_v_sse_30_192(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_addmul_v_sse_30_256(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_addmul_v_sse_3_128(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_addmul_v_sse_3_192(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_addmul_v_sse_3_256(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;

/**
 * Use AVX2
 */
void mzd_addmul_v_avx_30_128(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_addmul_v_avx_30_192(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_addmul_v_avx_30_256(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_addmul_v_avx_3_192(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_addmul_v_avx_3_256(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;

/**
 * Use NEON
 */
void mzd_addmul_v_neon_30_128(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_addmul_v_neon_30_192(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_addmul_v_neon_30_256(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_addmul_v_neon_3_128(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_addmul_v_neon_3_192(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_addmul_v_neon_3_256(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;

/**
 * Compute using parity based algorithm
 * */
void mzd_mul_v_parity_uint64_30(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_mul_v_parity_uint64_128_30(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_mul_v_parity_uint64_192_30(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_mul_v_parity_uint64_256_30(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_mul_v_parity_uint64_3(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_mul_v_parity_uint64_128_3(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_mul_v_parity_uint64_192_3(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_mul_v_parity_uint64_256_3(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;

/**
 * Using popcnt
 */
void mzd_mul_v_parity_popcnt_30(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_mul_v_parity_popcnt_128_30(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_mul_v_parity_popcnt_192_30(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_mul_v_parity_popcnt_256_30(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_mul_v_parity_popcnt_3(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_mul_v_parity_popcnt_128_3(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_mul_v_parity_popcnt_192_3(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_mul_v_parity_popcnt_256_3(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;

/**
 * Compute c + v * A optimized for c and v being vectors.
 */
void mzd_addmul_v_uint64(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_addmul_v_sse(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_addmul_v_sse_128(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_addmul_v_sse_192(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_addmul_v_sse_256(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_addmul_v_avx(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_addmul_v_avx_128(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_addmul_v_avx_192(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_addmul_v_avx_256(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_addmul_v_neon(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_addmul_v_neon_128(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_addmul_v_neon_192(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_addmul_v_neon_256(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;

/**
 * Compute v * A optimized for v being a vector.
 */
void mzd_mul_vl_uint64(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_mul_vl_sse_128(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_mul_vl_sse_192(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_mul_vl_sse_256(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_mul_vl_sse(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_mul_vl_avx_128(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_mul_vl_avx_192(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_mul_vl_avx_256(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_mul_vl_avx(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_mul_vl_neon_128(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_mul_vl_neon_192(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_mul_vl_neon_256(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_mul_vl_neon(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
/**
 * Compute c + v * A optimized for c and v being vectors.
 */
void mzd_addmul_vl_sse_128(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_addmul_vl_sse_192(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_addmul_vl_sse_256(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_addmul_vl_avx_128(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_addmul_vl_avx_192(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_addmul_vl_avx_256(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_addmul_vl_sse(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_addmul_vl_avx(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_addmul_vl_uint64(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_addmul_vl_neon_128(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;
void mzd_addmul_vl_neon(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) ATTR_NONNULL;

/**
 * Pre-compute matrices for mzd_{add,}mul_vl computions.
 */
mzd_local_t* mzd_precompute_matrix_lookup(mzd_local_t const* A) ATTR_NONNULL;

/**
 * Shuffle vector x according to info in mask. Needed for OLLE optimiztaions.
 */
void mzd_shuffle_30(mzd_local_t* x, const word mask) ATTR_NONNULL;
void mzd_shuffle_3(mzd_local_t* x, const word mask) ATTR_NONNULL;
void mzd_shuffle_pext_30(mzd_local_t* x, const word mask) ATTR_NONNULL;
void mzd_shuffle_pext_3(mzd_local_t* x, const word mask) ATTR_NONNULL;

#define ROW(v, r) (&(v)->rows[(v)->rowstride * (r)])
#define CONST_ROW(v, r) ((word const*)ROW(v, r))

#define FIRST_ROW(v) ROW(v, 0)
#define CONST_FIRST_ROW(v) CONST_ROW(v, 0)

#endif
