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

#if !defined(_MSC_VER)
#include <stdalign.h>
#endif
#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "mzd_additional.h"

static_assert(((sizeof(mzd_local_t) + 0x1f) & ~0x1f) == 32, "sizeof mzd_local_t not supported");

#if defined(WITH_OPT)
#include "simd.h"
#endif
static const unsigned int align_bound = 128 / (8 * sizeof(word));

static inline size_t calculate_rowstride(size_t width) {
  // As soon as we hit the AVX bound, use 32 byte alignment. Otherwise use 16
  // byte alignment for SSE2 and 128 bit vectors.
  if (width > align_bound) {
    return ((width * sizeof(word) + 31) & ~31) / sizeof(word);
  } else {
    return ((width * sizeof(word) + 15) & ~15) / sizeof(word);
  }
}

static inline size_t calculate_width(size_t c) {
  return (c + sizeof(word) * 8 - 1) / (sizeof(word) * 8);
}

// Notes on the memory layout: mzd_init allocates multiple memory blocks (one
// for mzd_local_t, one for rows and multiple for the buffers). We use one memory
// block for mzd_local_t, rows and the buffer. This improves memory locality and
// requires less calls to malloc.
//
// In mzd_local_init_multiple we do the same, but store n mzd_local_t instances in one
// memory block.

mzd_local_t* mzd_local_init_ex(unsigned int r, unsigned int c, bool clear) {
  const size_t rowstride = calculate_rowstride(calculate_width(c));

  const size_t buffer_size = r * rowstride * sizeof(word);
  const size_t alloc_size  = (buffer_size + 31) & ~31;

  /* We always align mzd_local_ts to 32 bytes. Thus the first row is always
   * aligned to 32 bytes as well. For 128 bit and SSE all other rows are then
   * aligned to 16 bytes. */
  unsigned char* buffer = picnic_aligned_alloc(32, alloc_size);
  if (clear) {
    memset(buffer, 0, alloc_size);
  }

  return (mzd_local_t*)buffer;
}

/* implementation of copy */

void mzd_copy_uint64_128(mzd_local_t* dst, mzd_local_t const* src) {
  const block_t* sblock = CONST_BLOCK(src, 0);
  block_t* dblock       = BLOCK(dst, 0);

  for (unsigned int i = 0; i < 2; ++i) {
    dblock->w64[i] = sblock->w64[i];
  }
}

void mzd_copy_uint64_192(mzd_local_t* dst, mzd_local_t const* src) {
  const block_t* sblock = CONST_BLOCK(src, 0);
  block_t* dblock       = BLOCK(dst, 0);

  for (unsigned int i = 0; i < 3; ++i) {
    dblock->w64[i] = sblock->w64[i];
  }
}

void mzd_copy_uint64_256(mzd_local_t* dst, mzd_local_t const* src) {
  const block_t* sblock = CONST_BLOCK(src, 0);
  block_t* dblock       = BLOCK(dst, 0);

  for (unsigned int i = 0; i < 4; ++i) {
    dblock->w64[i] = sblock->w64[i];
  }
}

#if defined(WITH_OPT)
#if defined(WITH_SSE2) || defined(WITH_NEON)
ATTR_TARGET_S128
void mzd_copy_s128_128(mzd_local_t* dst, mzd_local_t const* src) {
  mm128_store(BLOCK(dst, 0)->w64, mm128_load(CONST_BLOCK(src, 0)->w64));
}

ATTR_TARGET_S128
void mzd_copy_s128_256(mzd_local_t* dst, mzd_local_t const* src) {
  for (unsigned int i = 0; i < 2; ++i) {
    mm128_store(&dst->w64[2 * i], mm128_load(&src->w64[2 * i]));
  }
}
#endif

#if defined(WITH_AVX2)
ATTR_TARGET_AVX2
void mzd_copy_s256_128(mzd_local_t* dst, mzd_local_t const* src) {
  mm128_store(BLOCK(dst, 0)->w64, mm128_load(CONST_BLOCK(src, 0)->w64));
}

ATTR_TARGET_AVX2
void mzd_copy_s256_256(mzd_local_t* dst, mzd_local_t const* src) {
  mm256_store(BLOCK(dst, 0)->w64, mm256_load(CONST_BLOCK(src, 0)->w64));
}
#endif
#endif

/* implementation of mzd_xor and variants */

#if defined(WITH_OPT)
#if defined(WITH_SSE2) || defined(WITH_NEON)
ATTR_TARGET_S128
void mzd_xor_s128_128(mzd_local_t* res, mzd_local_t const* first, mzd_local_t const* second) {
  block_t* rblock       = BLOCK(res, 0);
  const block_t* fblock = CONST_BLOCK(first, 0);
  const block_t* sblock = CONST_BLOCK(second, 0);

  mm128_store(rblock->w64, mm128_xor(mm128_load(fblock->w64), mm128_load(sblock->w64)));
}

ATTR_TARGET_S128
static inline void mzd_xor_s128_blocks(block_t* rblock, const block_t* fblock,
                                       const block_t* sblock, unsigned int count) {
  for (; count; --count, ++rblock, ++fblock, ++sblock) {
    mm128_store(&rblock->w64[0],
                mm128_xor(mm128_load(&fblock->w64[0]), mm128_load(&sblock->w64[0])));
    mm128_store(&rblock->w64[2],
                mm128_xor(mm128_load(&fblock->w64[2]), mm128_load(&sblock->w64[2])));
  }
}

ATTR_TARGET_S128
void mzd_xor_s128_256(mzd_local_t* res, mzd_local_t const* first, mzd_local_t const* second) {
  mzd_xor_s128_blocks(BLOCK(res, 0), CONST_BLOCK(first, 0), CONST_BLOCK(second, 0), 1);
}

#if defined(WITH_LOWMC_128_128_20)
ATTR_TARGET_S128
void mzd_xor_s128_640(mzd_local_t* res, mzd_local_t const* first, mzd_local_t const* second) {
  mzd_xor_s128_blocks(BLOCK(res, 0), CONST_BLOCK(first, 0), CONST_BLOCK(second, 0), 2);
  mm128_store(BLOCK(res, 2)->w64, mm128_xor(mm128_load(CONST_BLOCK(first, 2)->w64),
                                            mm128_load(CONST_BLOCK(second, 2)->w64)));
}
#endif

#if defined(WITH_LOWMC_192_192_30)
ATTR_TARGET_S128
void mzd_xor_s128_1024(mzd_local_t* res, mzd_local_t const* first, mzd_local_t const* second) {
  mzd_xor_s128_blocks(BLOCK(res, 0), CONST_BLOCK(first, 0), CONST_BLOCK(second, 0), 4);
}
#endif

#if defined(WITH_LOWMC_256_256_38)
ATTR_TARGET_S128
void mzd_xor_s128_1280(mzd_local_t* res, mzd_local_t const* first, mzd_local_t const* second) {
  mzd_xor_s128_blocks(BLOCK(res, 0), CONST_BLOCK(first, 0), CONST_BLOCK(second, 0), 5);
}
#endif
#endif

#if defined(WITH_AVX2)
ATTR_TARGET_AVX2
void mzd_xor_s256_128(mzd_local_t* res, mzd_local_t const* first, mzd_local_t const* second) {
  block_t* rblock       = BLOCK(res, 0);
  const block_t* fblock = CONST_BLOCK(first, 0);
  const block_t* sblock = CONST_BLOCK(second, 0);

  mm128_store(rblock->w64, mm128_xor(mm128_load(fblock->w64), mm128_load(sblock->w64)));
}

ATTR_TARGET_AVX2
static inline void mzd_xor_s256_blocks(block_t* rblock, const block_t* fblock,
                                       const block_t* sblock, unsigned int count) {
  for (unsigned int idx = 0; idx != count; ++idx) {
    mm256_store(rblock[idx].w64,
                mm256_xor(mm256_load(fblock[idx].w64), mm256_load(sblock[idx].w64)));
  }
}

ATTR_TARGET_AVX2
void mzd_xor_s256_256(mzd_local_t* res, mzd_local_t const* first, mzd_local_t const* second) {
  mzd_xor_s256_blocks(BLOCK(res, 0), CONST_BLOCK(first, 0), CONST_BLOCK(second, 0), 1);
}

ATTR_TARGET_AVX2
void mzd_xor_s256_768(mzd_local_t* res, mzd_local_t const* first, mzd_local_t const* second) {
  mzd_xor_s256_blocks(BLOCK(res, 0), CONST_BLOCK(first, 0), CONST_BLOCK(second, 0), 3);
}

void mzd_xor_s256_1024(mzd_local_t* res, mzd_local_t const* first, mzd_local_t const* second) {
  mzd_xor_s256_blocks(BLOCK(res, 0), CONST_BLOCK(first, 0), CONST_BLOCK(second, 0), 4);
}

void mzd_xor_s256_1280(mzd_local_t* res, mzd_local_t const* first, mzd_local_t const* second) {
  mzd_xor_s256_blocks(BLOCK(res, 0), CONST_BLOCK(first, 0), CONST_BLOCK(second, 0), 5);
}
#endif
#endif

static inline void mzd_xor_uint64_block(block_t* rblock, const block_t* fblock,
                                        const block_t* sblock, const unsigned int len) {
  for (unsigned int i = 0; i < len; ++i) {
    rblock->w64[i] = fblock->w64[i] ^ sblock->w64[i];
  }
}

// clang-format off
#if defined(WITH_LOWMC_128_128_20) || defined(WITH_LOWMC_192_192_30) || defined(WITH_LOWMC_256_256_38)
// clang-format on
static inline void mzd_xor_uint64_blocks(block_t* rblock, const block_t* fblock,
                                         const block_t* sblock, const unsigned int len) {
  for (unsigned int i = len; i; --i, ++rblock, ++fblock, ++sblock) {
    mzd_xor_uint64_block(rblock, fblock, sblock, 4);
  }
}
#endif

void mzd_xor_uint64_128(mzd_local_t* res, mzd_local_t const* first, mzd_local_t const* second) {
  mzd_xor_uint64_block(BLOCK(res, 0), CONST_BLOCK(first, 0), CONST_BLOCK(second, 0), 2);
}

void mzd_xor_uint64_192(mzd_local_t* res, mzd_local_t const* first, mzd_local_t const* second) {
  mzd_xor_uint64_block(BLOCK(res, 0), CONST_BLOCK(first, 0), CONST_BLOCK(second, 0), 3);
}

void mzd_xor_uint64_256(mzd_local_t* res, mzd_local_t const* first, mzd_local_t const* second) {
  mzd_xor_uint64_block(BLOCK(res, 0), CONST_BLOCK(first, 0), CONST_BLOCK(second, 0), 4);
}

#if defined(WITH_LOWMC_128_128_20)
void mzd_xor_uint64_640(mzd_local_t* res, mzd_local_t const* first, mzd_local_t const* second) {
  mzd_xor_uint64_blocks(BLOCK(res, 0), CONST_BLOCK(first, 0), CONST_BLOCK(second, 0), 2);
  mzd_xor_uint64_block(BLOCK(res, 2), CONST_BLOCK(first, 2), CONST_BLOCK(second, 2), 2);
}
#endif

#if defined(WITH_LOWMC_192_192_30)
void mzd_xor_uint64_960(mzd_local_t* res, mzd_local_t const* first, mzd_local_t const* second) {
  mzd_xor_uint64_blocks(BLOCK(res, 0), CONST_BLOCK(first, 0), CONST_BLOCK(second, 0), 3);
  mzd_xor_uint64_block(BLOCK(res, 3), CONST_BLOCK(first, 3), CONST_BLOCK(second, 3), 3);
}
#endif

#if defined(WITH_LOWMC_256_256_38)
void mzd_xor_uint64_1216(mzd_local_t* res, mzd_local_t const* first, mzd_local_t const* second) {
  mzd_xor_uint64_blocks(BLOCK(res, 0), CONST_BLOCK(first, 0), CONST_BLOCK(second, 0), 4);
  mzd_xor_uint64_block(BLOCK(res, 4), CONST_BLOCK(first, 4), CONST_BLOCK(second, 4), 3);
}
#endif

/* implementation of mzd_and_* and variants */

static inline void mzd_and_uint64_block(block_t* rblock, const block_t* fblock,
                                        const block_t* sblock, const unsigned int len) {
  for (unsigned int i = 0; i < len; ++i) {
    rblock->w64[i] = fblock->w64[i] & sblock->w64[i];
  }
}

void mzd_and_uint64_128(mzd_local_t* res, mzd_local_t const* first, mzd_local_t const* second) {
  mzd_and_uint64_block(BLOCK(res, 0), CONST_BLOCK(first, 0), CONST_BLOCK(second, 0), 2);
}

void mzd_and_uint64_192(mzd_local_t* res, mzd_local_t const* first, mzd_local_t const* second) {
  mzd_and_uint64_block(BLOCK(res, 0), CONST_BLOCK(first, 0), CONST_BLOCK(second, 0), 3);
}

void mzd_and_uint64_256(mzd_local_t* res, mzd_local_t const* first, mzd_local_t const* second) {
  mzd_and_uint64_block(BLOCK(res, 0), CONST_BLOCK(first, 0), CONST_BLOCK(second, 0), 4);
}

/* shifts and rotations */

void mzd_shift_left_uint64_128(mzd_local_t* res, const mzd_local_t* val, unsigned int count) {
  const unsigned int right_count = 8 * sizeof(word) - count;
  const block_t* block           = CONST_BLOCK(val, 0);
  block_t* rblock                = BLOCK(res, 0);

  rblock->w64[1] = (block->w64[1] << count) | (block->w64[0] >> right_count);
  rblock->w64[0] = block->w64[0] << count;
}

void mzd_shift_right_uint64_128(mzd_local_t* res, const mzd_local_t* val, unsigned int count) {
  const unsigned int left_count = 8 * sizeof(word) - count;
  const block_t* block          = CONST_BLOCK(val, 0);
  block_t* rblock               = BLOCK(res, 0);

  rblock->w64[0] = (block->w64[0] >> count) | (block->w64[1] << left_count);
  rblock->w64[1] = block->w64[1] >> count;
}

void mzd_shift_left_uint64_192(mzd_local_t* res, const mzd_local_t* val, unsigned int count) {
  const unsigned int right_count = 8 * sizeof(word) - count;
  const block_t* block           = CONST_BLOCK(val, 0);
  block_t* rblock                = BLOCK(res, 0);

  rblock->w64[2] = (block->w64[2] << count) | (block->w64[1] >> right_count);
  rblock->w64[1] = (block->w64[1] << count) | (block->w64[0] >> right_count);
  rblock->w64[0] = block->w64[0] << count;
}

void mzd_shift_right_uint64_192(mzd_local_t* res, const mzd_local_t* val, unsigned int count) {
  const unsigned int left_count = 8 * sizeof(word) - count;
  const block_t* block          = CONST_BLOCK(val, 0);
  block_t* rblock               = BLOCK(res, 0);

  rblock->w64[0] = (block->w64[0] >> count) | (block->w64[1] << left_count);
  rblock->w64[1] = (block->w64[1] >> count) | (block->w64[2] << left_count);
  rblock->w64[2] = block->w64[2] >> count;
}

void mzd_shift_left_uint64_256(mzd_local_t* res, const mzd_local_t* val, unsigned int count) {
  const unsigned int right_count = 8 * sizeof(word) - count;
  const block_t* block           = CONST_BLOCK(val, 0);
  block_t* rblock                = BLOCK(res, 0);

  rblock->w64[3] = (block->w64[3] << count) | (block->w64[2] >> right_count);
  rblock->w64[2] = (block->w64[2] << count) | (block->w64[1] >> right_count);
  rblock->w64[1] = (block->w64[1] << count) | (block->w64[0] >> right_count);
  rblock->w64[0] = block->w64[0] << count;
}

void mzd_shift_right_uint64_256(mzd_local_t* res, const mzd_local_t* val, unsigned int count) {
  const unsigned int left_count = 8 * sizeof(word) - count;
  const block_t* block          = CONST_BLOCK(val, 0);
  block_t* rblock               = BLOCK(res, 0);

  rblock->w64[0] = (block->w64[0] >> count) | (block->w64[1] << left_count);
  rblock->w64[1] = (block->w64[1] >> count) | (block->w64[2] << left_count);
  rblock->w64[2] = (block->w64[2] >> count) | (block->w64[3] << left_count);
  rblock->w64[3] = block->w64[3] >> count;
}

#if defined(PICNIC_STATIC)
void mzd_rotate_left_uint64_128(mzd_local_t* res, const mzd_local_t* val, unsigned int count) {
  const unsigned int right_count = 8 * sizeof(word) - count;
  const block_t* block           = CONST_BLOCK(val, 0);
  block_t* rblock                = BLOCK(res, 0);

  const word tmp = block->w64[1] >> right_count;
  rblock->w64[1] = (block->w64[1] << count) | (block->w64[0] >> right_count);
  rblock->w64[0] = (block->w64[0] << count) | tmp;
}

void mzd_rotate_right_uint64_128(mzd_local_t* res, const mzd_local_t* val, unsigned int count) {
  const unsigned int left_count = 8 * sizeof(word) - count;
  const block_t* block          = CONST_BLOCK(val, 0);
  block_t* rblock               = BLOCK(res, 0);

  const word tmp = block->w64[0] << left_count;
  rblock->w64[0] = (block->w64[0] >> count) | (block->w64[1] << left_count);
  rblock->w64[1] = (block->w64[1] >> count) | tmp;
}

void mzd_rotate_left_uint64_192(mzd_local_t* res, const mzd_local_t* val, unsigned int count) {
  const unsigned int right_count = 8 * sizeof(word) - count;
  const block_t* block           = CONST_BLOCK(val, 0);
  block_t* rblock                = BLOCK(res, 0);

  const word tmp = block->w64[2] >> right_count;
  rblock->w64[2] = (block->w64[2] << count) | (block->w64[1] >> right_count);
  rblock->w64[1] = (block->w64[1] << count) | (block->w64[0] >> right_count);
  rblock->w64[0] = (block->w64[0] << count) | tmp;
}

void mzd_rotate_right_uint64_192(mzd_local_t* res, const mzd_local_t* val, unsigned int count) {
  const unsigned int left_count = 8 * sizeof(word) - count;
  const block_t* block          = CONST_BLOCK(val, 0);
  block_t* rblock               = BLOCK(res, 0);

  const word tmp = block->w64[0] << left_count;
  rblock->w64[0] = (block->w64[0] >> count) | (block->w64[1] << left_count);
  rblock->w64[1] = (block->w64[1] >> count) | (block->w64[2] << left_count);
  rblock->w64[2] = (block->w64[2] >> count) | tmp;
}

void mzd_rotate_left_uint64_256(mzd_local_t* res, const mzd_local_t* val, unsigned int count) {
  const unsigned int right_count = 8 * sizeof(word) - count;
  const block_t* block           = CONST_BLOCK(val, 0);
  block_t* rblock                = BLOCK(res, 0);

  const word tmp = block->w64[3] >> right_count;
  rblock->w64[3] = (block->w64[3] << count) | (block->w64[2] >> right_count);
  rblock->w64[2] = (block->w64[2] << count) | (block->w64[1] >> right_count);
  rblock->w64[1] = (block->w64[1] << count) | (block->w64[0] >> right_count);
  rblock->w64[0] = (block->w64[0] << count) | tmp;
}

void mzd_rotate_right_uint64_256(mzd_local_t* res, const mzd_local_t* val, unsigned int count) {
  const unsigned int left_count = 8 * sizeof(word) - count;
  const block_t* block          = CONST_BLOCK(val, 0);
  block_t* rblock               = BLOCK(res, 0);

  const word tmp = block->w64[0] << left_count;
  rblock->w64[0] = (block->w64[0] >> count) | (block->w64[1] << left_count);
  rblock->w64[1] = (block->w64[1] >> count) | (block->w64[2] << left_count);
  rblock->w64[2] = (block->w64[2] >> count) | (block->w64[3] << left_count);
  rblock->w64[3] = (block->w64[3] >> count) | tmp;
}
#endif

#if defined(WITH_LOWMC_128_128_20)
void mzd_mul_v_parity_uint64_128_30(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* At) {
  block_t* cblock       = BLOCK(c, 0);
  const block_t* vblock = CONST_BLOCK(v, 0);

  word res = 0;
  for (unsigned int i = 15; i; --i) {
    const block_t* Ablock = CONST_BLOCK(At, 15 - i);
    const word parity1 =
        parity64_uint64((vblock->w64[0] & Ablock->w64[0]) ^ (vblock->w64[1] & Ablock->w64[1]));
    const word parity2 =
        parity64_uint64((vblock->w64[0] & Ablock->w64[2]) ^ (vblock->w64[1] & Ablock->w64[3]));
    res |= (parity1 | (parity2 << 1)) << (64 - (2 * i));
  }
  cblock->w64[0] = 0;
  cblock->w64[1] = res;
}
#endif

#if defined(WITH_LOWMC_192_192_30)
void mzd_mul_v_parity_uint64_192_30(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* At) {
  block_t* cblock       = BLOCK(c, 0);
  const block_t* vblock = CONST_BLOCK(v, 0);

  word res = 0;
  for (unsigned int i = 30; i; --i) {
    const block_t* Ablock = CONST_BLOCK(At, 30 - i);
    const word parity =
        parity64_uint64((vblock->w64[0] & Ablock->w64[0]) ^ (vblock->w64[1] & Ablock->w64[1]) ^
                        (vblock->w64[2] & Ablock->w64[2]));
    res |= parity << (64 - i);
  }
  for (unsigned int j = 0; j < 2; j++) {
    cblock->w64[j] = 0;
  }
  cblock->w64[2] = res;
}
#endif

#if defined(WITH_LOWMC_256_256_38)
void mzd_mul_v_parity_uint64_256_30(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* At) {
  block_t* cblock       = BLOCK(c, 0);
  const block_t* vblock = CONST_BLOCK(v, 0);

  word res = 0;
  for (unsigned int i = 30; i; --i) {
    const block_t* Ablock = CONST_BLOCK(At, 30 - i);
    const word parity =
        parity64_uint64((vblock->w64[0] & Ablock->w64[0]) ^ (vblock->w64[1] & Ablock->w64[1]) ^
                        (vblock->w64[2] & Ablock->w64[2]) ^ (vblock->w64[3] & Ablock->w64[3]));
    res |= parity << (64 - i);
  }
  for (unsigned int j = 0; j < 3; j++) {
    cblock->w64[j] = 0;
  }
  cblock->w64[3] = res;
}
#endif

#if defined(WITH_OPT)
#if defined(WITH_SSE2) || defined(WITH_NEON)
ATTR_TARGET_S128 ATTR_ARTIFICIAL ATTR_CONST static inline word128
mm128_compute_mask(const word idx, const size_t bit) {
  return mm128_broadcast_u64(-((idx >> bit) & 1));
}

ATTR_TARGET_S128
void mzd_mul_v_s128_128(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  block_t* cblock       = BLOCK(c, 0);
  const word* vptr      = CONST_BLOCK(v, 0)->w64;
  const block_t* Ablock = CONST_BLOCK(A, 0);

  word128 cval[4] = {mm128_zero, mm128_zero, mm128_zero, mm128_zero};
  for (unsigned int w = 2; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int i = sizeof(word) * 8; i; i -= 4, idx >>= 4, Ablock += 2) {
      cval[0] = mm128_xor_mask(cval[0], mm128_load(&Ablock[0].w64[0]), mm128_compute_mask(idx, 0));
      cval[1] = mm128_xor_mask(cval[1], mm128_load(&Ablock[0].w64[2]), mm128_compute_mask(idx, 1));
      cval[0] = mm128_xor_mask(cval[0], mm128_load(&Ablock[1].w64[0]), mm128_compute_mask(idx, 2));
      cval[1] = mm128_xor_mask(cval[1], mm128_load(&Ablock[1].w64[2]), mm128_compute_mask(idx, 3));
    }
  }
  mm128_store(cblock->w64, mm128_xor(cval[0], cval[1]));
}

ATTR_TARGET_S128
void mzd_addmul_v_s128_128(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  block_t* cblock       = BLOCK(c, 0);
  const word* vptr      = CONST_BLOCK(v, 0)->w64;
  const block_t* Ablock = CONST_BLOCK(A, 0);

  word128 cval[4] = {mm128_load(cblock->w64), mm128_zero, mm128_zero, mm128_zero};
  for (unsigned int w = 2; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int i = sizeof(word) * 8; i; i -= 4, idx >>= 4, Ablock += 2) {
      cval[0] = mm128_xor_mask(cval[0], mm128_load(&Ablock[0].w64[0]), mm128_compute_mask(idx, 0));
      cval[1] = mm128_xor_mask(cval[1], mm128_load(&Ablock[0].w64[2]), mm128_compute_mask(idx, 1));
      cval[0] = mm128_xor_mask(cval[0], mm128_load(&Ablock[1].w64[0]), mm128_compute_mask(idx, 2));
      cval[1] = mm128_xor_mask(cval[1], mm128_load(&Ablock[1].w64[2]), mm128_compute_mask(idx, 3));
    }
  }
  mm128_store(cblock->w64, mm128_xor(cval[0], cval[1]));
}

ATTR_TARGET_S128
void mzd_mul_v_s128_129(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  block_t* cblock       = BLOCK(c, 0);
  const word* vptr      = CONST_BLOCK(v, 0)->w64;
  const block_t* Ablock = CONST_BLOCK(A, 0);

  word128 cval[4] = {mm128_zero, mm128_zero, mm128_zero, mm128_zero};
  {
    Ablock += 63;
    word idx = (*vptr) >> 63;
    mm128_xor_mask_region(&cval[0], Ablock[0].w64, mm128_compute_mask(idx, 0), 2);
    vptr++;
    Ablock++;
  }
  for (unsigned int w = 2; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int i = sizeof(word) * 8; i; i -= 2, idx >>= 2, Ablock += 2) {
      mm128_xor_mask_region(&cval[0], Ablock[0].w64, mm128_compute_mask(idx, 0), 2);
      mm128_xor_mask_region(&cval[2], Ablock[1].w64, mm128_compute_mask(idx, 1), 2);
    }
  }
  mm128_store(&cblock->w64[0], mm128_xor(cval[0], cval[2]));
  mm128_store(&cblock->w64[2], mm128_xor(cval[1], cval[3]));
}

ATTR_TARGET_S128
void mzd_addmul_v_s128_129(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  block_t* cblock       = BLOCK(c, 0);
  const word* vptr      = CONST_BLOCK(v, 0)->w64;
  const block_t* Ablock = CONST_BLOCK(A, 0);

  word128 cval[4] = {mm128_load(&cblock->w64[0]), mm128_load(&cblock->w64[2]), mm128_zero,
                     mm128_zero};
  {
    Ablock += 63;
    word idx = (*vptr) >> 63;
    mm128_xor_mask_region(&cval[0], Ablock[0].w64, mm128_compute_mask(idx, 0), 2);
    vptr++;
    Ablock++;
  }
  for (unsigned int w = 2; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int i = sizeof(word) * 8; i; i -= 2, idx >>= 2, Ablock += 2) {
      mm128_xor_mask_region(&cval[0], Ablock[0].w64, mm128_compute_mask(idx, 0), 2);
      mm128_xor_mask_region(&cval[2], Ablock[1].w64, mm128_compute_mask(idx, 1), 2);
    }
  }
  mm128_store(&cblock->w64[0], mm128_xor(cval[0], cval[2]));
  mm128_store(&cblock->w64[2], mm128_xor(cval[1], cval[3]));
}

ATTR_TARGET_S128
void mzd_mul_v_s128_192(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  block_t* cblock       = BLOCK(c, 0);
  const word* vptr      = CONST_BLOCK(v, 0)->w64;
  const block_t* Ablock = CONST_BLOCK(A, 0);

  word128 cval[4] = {mm128_zero, mm128_zero, mm128_zero, mm128_zero};
  for (unsigned int w = 3; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int i = sizeof(word) * 8; i; i -= 2, idx >>= 2, Ablock += 2) {
      mm128_xor_mask_region(&cval[0], Ablock[0].w64, mm128_compute_mask(idx, 0), 2);
      mm128_xor_mask_region(&cval[2], Ablock[1].w64, mm128_compute_mask(idx, 1), 2);
    }
  }
  mm128_store(&cblock->w64[0], mm128_xor(cval[0], cval[2]));
  mm128_store(&cblock->w64[2], mm128_xor(cval[1], cval[3]));
}

ATTR_TARGET_S128
void mzd_addmul_v_s128_192(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  block_t* cblock       = BLOCK(c, 0);
  const word* vptr      = CONST_BLOCK(v, 0)->w64;
  const block_t* Ablock = CONST_BLOCK(A, 0);

  word128 cval[4] = {mm128_load(&cblock->w64[0]), mm128_load(&cblock->w64[2]), mm128_zero,
                     mm128_zero};
  for (unsigned int w = 3; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int i = sizeof(word) * 8; i; i -= 2, idx >>= 2, Ablock += 2) {
      mm128_xor_mask_region(&cval[0], Ablock[0].w64, mm128_compute_mask(idx, 0), 2);
      mm128_xor_mask_region(&cval[2], Ablock[1].w64, mm128_compute_mask(idx, 1), 2);
    }
  }
  mm128_store(&cblock->w64[0], mm128_xor(cval[0], cval[2]));
  mm128_store(&cblock->w64[2], mm128_xor(cval[1], cval[3]));
}

ATTR_TARGET_S128
void mzd_mul_v_s128_256(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  block_t* cblock       = BLOCK(c, 0);
  const word* vptr      = CONST_BLOCK(v, 0)->w64;
  const block_t* Ablock = CONST_BLOCK(A, 0);

  word128 cval[4] = {mm128_zero, mm128_zero, mm128_zero, mm128_zero};
  for (unsigned int w = 4; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int i = sizeof(word) * 8; i; i -= 2, idx >>= 2, Ablock += 2) {
      mm128_xor_mask_region(&cval[0], Ablock[0].w64, mm128_compute_mask(idx, 0), 2);
      mm128_xor_mask_region(&cval[2], Ablock[1].w64, mm128_compute_mask(idx, 1), 2);
    }
  }
  mm128_store(&cblock->w64[0], mm128_xor(cval[0], cval[2]));
  mm128_store(&cblock->w64[2], mm128_xor(cval[1], cval[3]));
}

ATTR_TARGET_S128
void mzd_addmul_v_s128_256(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  block_t* cblock       = BLOCK(c, 0);
  const word* vptr      = CONST_BLOCK(v, 0)->w64;
  const block_t* Ablock = CONST_BLOCK(A, 0);

  word128 cval[4] = {mm128_load(&cblock->w64[0]), mm128_load(&cblock->w64[2]), mm128_zero,
                     mm128_zero};
  for (unsigned int w = 4; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int i = sizeof(word) * 8; i; i -= 2, idx >>= 2, Ablock += 2) {
      mm128_xor_mask_region(&cval[0], Ablock[0].w64, mm128_compute_mask(idx, 0), 2);
      mm128_xor_mask_region(&cval[2], Ablock[1].w64, mm128_compute_mask(idx, 1), 2);
    }
  }
  mm128_store(&cblock->w64[0], mm128_xor(cval[0], cval[2]));
  mm128_store(&cblock->w64[2], mm128_xor(cval[1], cval[3]));
}

#if defined(WITH_LOWMC_128_128_20)
ATTR_TARGET_S128
void mzd_mul_v_s128_128_640(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  const word* vptr      = CONST_BLOCK(v, 0)->w64;
  const block_t* Ablock = CONST_BLOCK(A, 0);

  word128 cval[5] = {mm128_zero, mm128_zero, mm128_zero, mm128_zero, mm128_zero};
  for (unsigned int w = 2; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int i = sizeof(word) * 8; i; i -= 1, idx >>= 1, Ablock += 3) {
      const word128 mask = mm128_compute_mask(idx, 0);
      mm128_xor_mask_region(&cval[0], Ablock[0].w64, mask, 2);
      mm128_xor_mask_region(&cval[2], Ablock[1].w64, mask, 2);
      cval[4] = mm128_xor_mask(cval[4], mm128_load(Ablock[2].w64), mask);
    }
  }

  block_t* cblock1 = BLOCK(c, 0);
  block_t* cblock2 = BLOCK(c, 1);
  block_t* cblock3 = BLOCK(c, 2);
  mm128_store(&cblock1->w64[0], cval[0]);
  mm128_store(&cblock1->w64[2], cval[1]);
  mm128_store(&cblock2->w64[0], cval[2]);
  mm128_store(&cblock2->w64[2], cval[3]);
  mm128_store(&cblock3->w64[0], cval[4]);
}
#endif

#if defined(WITH_LOWMC_192_192_30)
ATTR_TARGET_S128
void mzd_mul_v_s128_192_1024(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  const word* vptr      = CONST_BLOCK(v, 0)->w64;
  const block_t* Ablock = CONST_BLOCK(A, 0);

  word128 cval[8] = {mm128_zero, mm128_zero, mm128_zero, mm128_zero,
                     mm128_zero, mm128_zero, mm128_zero, mm128_zero};
  for (unsigned int w = 3; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int i = sizeof(word) * 8; i; i -= 1, idx >>= 1, Ablock += 4) {
      const word128 mask = mm128_compute_mask(idx, 0);
      mm128_xor_mask_region(&cval[0], Ablock[0].w64, mask, 2);
      mm128_xor_mask_region(&cval[2], Ablock[1].w64, mask, 2);
      mm128_xor_mask_region(&cval[4], Ablock[2].w64, mask, 2);
      mm128_xor_mask_region(&cval[6], Ablock[3].w64, mask, 2);
    }
  }

  block_t* cblock1 = BLOCK(c, 0);
  block_t* cblock2 = BLOCK(c, 1);
  block_t* cblock3 = BLOCK(c, 2);
  block_t* cblock4 = BLOCK(c, 3);
  mm128_store(&cblock1->w64[0], cval[0]);
  mm128_store(&cblock1->w64[2], cval[1]);
  mm128_store(&cblock2->w64[0], cval[2]);
  mm128_store(&cblock2->w64[2], cval[3]);
  mm128_store(&cblock3->w64[0], cval[4]);
  mm128_store(&cblock3->w64[2], cval[5]);
  mm128_store(&cblock4->w64[0], cval[6]);
  mm128_store(&cblock4->w64[2], cval[7]);
}
#endif

#if defined(WITH_LOWMC_256_256_38)
ATTR_TARGET_S128
void mzd_mul_v_s128_256_1280(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  const word* vptr      = CONST_BLOCK(v, 0)->w64;
  const block_t* Ablock = CONST_BLOCK(A, 0);

  word128 cval[10] = {mm128_zero, mm128_zero, mm128_zero, mm128_zero, mm128_zero,
                      mm128_zero, mm128_zero, mm128_zero, mm128_zero, mm128_zero};
  for (unsigned int w = 4; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int i = sizeof(word) * 8; i; i -= 1, idx >>= 1, Ablock += 5) {
      const word128 mask = mm128_compute_mask(idx, 0);
      mm128_xor_mask_region(&cval[0], Ablock[0].w64, mask, 2);
      mm128_xor_mask_region(&cval[2], Ablock[1].w64, mask, 2);
      mm128_xor_mask_region(&cval[4], Ablock[2].w64, mask, 2);
      mm128_xor_mask_region(&cval[6], Ablock[3].w64, mask, 2);
      mm128_xor_mask_region(&cval[8], Ablock[4].w64, mask, 2);
    }
  }

  block_t* cblock1 = BLOCK(c, 0);
  block_t* cblock2 = BLOCK(c, 1);
  block_t* cblock3 = BLOCK(c, 2);
  block_t* cblock4 = BLOCK(c, 3);
  block_t* cblock5 = BLOCK(c, 4);
  mm128_store(&cblock1->w64[0], cval[0]);
  mm128_store(&cblock1->w64[2], cval[1]);
  mm128_store(&cblock2->w64[0], cval[2]);
  mm128_store(&cblock2->w64[2], cval[3]);
  mm128_store(&cblock3->w64[0], cval[4]);
  mm128_store(&cblock3->w64[2], cval[5]);
  mm128_store(&cblock4->w64[0], cval[6]);
  mm128_store(&cblock4->w64[2], cval[7]);
  mm128_store(&cblock5->w64[0], cval[8]);
  mm128_store(&cblock5->w64[2], cval[9]);
}
#endif
#endif

#if defined(WITH_AVX2)
ATTR_TARGET_AVX2 ATTR_ARTIFICIAL ATTR_CONST static inline word256
mm256_compute_mask(const word idx, const size_t bit) {
  return _mm256_set1_epi64x(-((idx >> bit) & 1));
}

ATTR_TARGET_AVX2 ATTR_ARTIFICIAL ATTR_CONST static inline word256
mm256_compute_mask_2(const word idx, const size_t bit) {
  const word m1 = -((idx >> bit) & 1);
  const word m2 = -((idx >> (bit + 1)) & 1);
  return _mm256_set_epi64x(m2, m2, m1, m1);
}

ATTR_TARGET_AVX2
void mzd_addmul_v_s256_128(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  block_t* cblock       = BLOCK(c, 0);
  const word* vptr      = CONST_BLOCK(v, 0)->w64;
  const block_t* Ablock = CONST_BLOCK(A, 0);

  word256 cval[2] = {_mm256_setr_m128i(mm128_load(cblock->w64), mm128_zero), mm256_zero};
  for (unsigned int w = 2; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int i = sizeof(word) * 8; i; i -= 8, idx >>= 8, Ablock += 4) {
      cval[0] = mm256_xor_mask(cval[0], mm256_load(Ablock[0].w64), mm256_compute_mask_2(idx, 0));
      cval[1] = mm256_xor_mask(cval[1], mm256_load(Ablock[1].w64), mm256_compute_mask_2(idx, 2));
      cval[0] = mm256_xor_mask(cval[0], mm256_load(Ablock[2].w64), mm256_compute_mask_2(idx, 4));
      cval[1] = mm256_xor_mask(cval[1], mm256_load(Ablock[3].w64), mm256_compute_mask_2(idx, 6));
    }
  }
  cval[0] = mm256_xor(cval[0], cval[1]);
  mm128_store(
      cblock->w64,
      _mm256_extracti128_si256(
          mm256_xor(cval[0], _mm256_permute4x64_epi64(cval[0], _MM_SHUFFLE(3, 2, 3, 2))), 0));
}

ATTR_TARGET_AVX2
void mzd_mul_v_s256_128(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  block_t* cblock       = BLOCK(c, 0);
  const word* vptr      = CONST_BLOCK(v, 0)->w64;
  const block_t* Ablock = CONST_BLOCK(A, 0);

  word256 cval[2] = {mm256_zero, mm256_zero};
  for (unsigned int w = 2; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int i = sizeof(word) * 8; i; i -= 8, idx >>= 8, Ablock += 4) {
      cval[0] = mm256_xor_mask(cval[0], mm256_load(Ablock[0].w64), mm256_compute_mask_2(idx, 0));
      cval[1] = mm256_xor_mask(cval[1], mm256_load(Ablock[1].w64), mm256_compute_mask_2(idx, 2));
      cval[0] = mm256_xor_mask(cval[0], mm256_load(Ablock[2].w64), mm256_compute_mask_2(idx, 4));
      cval[1] = mm256_xor_mask(cval[1], mm256_load(Ablock[3].w64), mm256_compute_mask_2(idx, 6));
    }
  }
  cval[0] = mm256_xor(cval[0], cval[1]);
  mm128_store(
      cblock->w64,
      _mm256_extracti128_si256(
          mm256_xor(cval[0], _mm256_permute4x64_epi64(cval[0], _MM_SHUFFLE(3, 2, 3, 2))), 0));
}

ATTR_TARGET_AVX2
void mzd_addmul_v_s256_129(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  block_t* cblock       = BLOCK(c, 0);
  const word* vptr      = CONST_BLOCK(v, 0)->w64;
  const block_t* Ablock = CONST_BLOCK(A, 0);

  word256 cval[2] = {mm256_load(cblock->w64), mm256_zero};
  {
    cval[0] =
        mm256_xor_mask(cval[0], mm256_load(Ablock[63].w64), mm256_compute_mask((*vptr) >> 63, 0));
    vptr++;
    Ablock += 64;
  }
  for (unsigned int w = 2; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int i = sizeof(word) * 8; i; i -= 4, idx >>= 4, Ablock += 4) {
      cval[0] = mm256_xor_mask(cval[0], mm256_load(Ablock[0].w64), mm256_compute_mask(idx, 0));
      cval[1] = mm256_xor_mask(cval[1], mm256_load(Ablock[1].w64), mm256_compute_mask(idx, 1));
      cval[0] = mm256_xor_mask(cval[0], mm256_load(Ablock[2].w64), mm256_compute_mask(idx, 2));
      cval[1] = mm256_xor_mask(cval[1], mm256_load(Ablock[3].w64), mm256_compute_mask(idx, 3));
    }
  }
  mm256_store(cblock->w64, mm256_xor(cval[0], cval[1]));
}

ATTR_TARGET_AVX2
void mzd_mul_v_s256_129(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  block_t* cblock       = BLOCK(c, 0);
  const word* vptr      = CONST_BLOCK(v, 0)->w64;
  const block_t* Ablock = CONST_BLOCK(A, 0);

  word256 cval[2] = {mm256_zero, mm256_zero};
  {
    cval[0] =
        mm256_xor_mask(cval[0], mm256_load(Ablock[63].w64), mm256_compute_mask((*vptr) >> 63, 0));
    vptr++;
    Ablock += 64;
  }
  for (unsigned int w = 2; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int i = sizeof(word) * 8; i; i -= 4, idx >>= 4, Ablock += 4) {
      cval[0] = mm256_xor_mask(cval[0], mm256_load(Ablock[0].w64), mm256_compute_mask(idx, 0));
      cval[1] = mm256_xor_mask(cval[1], mm256_load(Ablock[1].w64), mm256_compute_mask(idx, 1));
      cval[0] = mm256_xor_mask(cval[0], mm256_load(Ablock[2].w64), mm256_compute_mask(idx, 2));
      cval[1] = mm256_xor_mask(cval[1], mm256_load(Ablock[3].w64), mm256_compute_mask(idx, 3));
    }
  }
  mm256_store(cblock->w64, mm256_xor(cval[0], cval[1]));
}

ATTR_TARGET_AVX2
void mzd_addmul_v_s256_192(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  block_t* cblock       = BLOCK(c, 0);
  const word* vptr      = CONST_BLOCK(v, 0)->w64;
  const block_t* Ablock = CONST_BLOCK(A, 0);

  word256 cval[2] = {mm256_load(cblock->w64), mm256_zero};
  for (unsigned int w = 3; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int i = sizeof(word) * 8; i; i -= 4, idx >>= 4, Ablock += 4) {
      cval[0] = mm256_xor_mask(cval[0], mm256_load(Ablock[0].w64), mm256_compute_mask(idx, 0));
      cval[1] = mm256_xor_mask(cval[1], mm256_load(Ablock[1].w64), mm256_compute_mask(idx, 1));
      cval[0] = mm256_xor_mask(cval[0], mm256_load(Ablock[2].w64), mm256_compute_mask(idx, 2));
      cval[1] = mm256_xor_mask(cval[1], mm256_load(Ablock[3].w64), mm256_compute_mask(idx, 3));
    }
  }
  mm256_store(cblock->w64, mm256_xor(cval[0], cval[1]));
}

ATTR_TARGET_AVX2
void mzd_mul_v_s256_192(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  block_t* cblock       = BLOCK(c, 0);
  const word* vptr      = CONST_BLOCK(v, 0)->w64;
  const block_t* Ablock = CONST_BLOCK(A, 0);

  word256 cval[2] = {mm256_zero, mm256_zero};
  for (unsigned int w = 3; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int i = sizeof(word) * 8; i; i -= 4, idx >>= 4, Ablock += 4) {
      cval[0] = mm256_xor_mask(cval[0], mm256_load(Ablock[0].w64), mm256_compute_mask(idx, 0));
      cval[1] = mm256_xor_mask(cval[1], mm256_load(Ablock[1].w64), mm256_compute_mask(idx, 1));
      cval[0] = mm256_xor_mask(cval[0], mm256_load(Ablock[2].w64), mm256_compute_mask(idx, 2));
      cval[1] = mm256_xor_mask(cval[1], mm256_load(Ablock[3].w64), mm256_compute_mask(idx, 3));
    }
  }
  mm256_store(cblock->w64, mm256_xor(cval[0], cval[1]));
}

ATTR_TARGET_AVX2
void mzd_addmul_v_s256_256(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  block_t* cblock       = BLOCK(c, 0);
  const word* vptr      = CONST_BLOCK(v, 0)->w64;
  const block_t* Ablock = CONST_BLOCK(A, 0);

  word256 cval[2] = {mm256_load(cblock->w64), mm256_zero};
  for (unsigned int w = 4; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int i = sizeof(word) * 8; i; i -= 4, idx >>= 4, Ablock += 4) {
      cval[0] = mm256_xor_mask(cval[0], mm256_load(Ablock[0].w64), mm256_compute_mask(idx, 0));
      cval[1] = mm256_xor_mask(cval[1], mm256_load(Ablock[1].w64), mm256_compute_mask(idx, 1));
      cval[0] = mm256_xor_mask(cval[0], mm256_load(Ablock[2].w64), mm256_compute_mask(idx, 2));
      cval[1] = mm256_xor_mask(cval[1], mm256_load(Ablock[3].w64), mm256_compute_mask(idx, 3));
    }
  }
  mm256_store(cblock->w64, mm256_xor(cval[0], cval[1]));
}

ATTR_TARGET_AVX2
void mzd_mul_v_s256_256(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  block_t* cblock       = BLOCK(c, 0);
  const word* vptr      = CONST_BLOCK(v, 0)->w64;
  const block_t* Ablock = CONST_BLOCK(A, 0);

  word256 cval[2] = {mm256_zero, mm256_zero};
  for (unsigned int w = 4; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int i = sizeof(word) * 8; i; i -= 4, idx >>= 4, Ablock += 4) {
      cval[0] = mm256_xor_mask(cval[0], mm256_load(Ablock[0].w64), mm256_compute_mask(idx, 0));
      cval[1] = mm256_xor_mask(cval[1], mm256_load(Ablock[1].w64), mm256_compute_mask(idx, 1));
      cval[0] = mm256_xor_mask(cval[0], mm256_load(Ablock[2].w64), mm256_compute_mask(idx, 2));
      cval[1] = mm256_xor_mask(cval[1], mm256_load(Ablock[3].w64), mm256_compute_mask(idx, 3));
    }
  }
  mm256_store(cblock->w64, mm256_xor(cval[0], cval[1]));
}

#if defined(WITH_LOWMC_128_128_20)
ATTR_TARGET_AVX2
void mzd_mul_v_s256_128_768(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  const word* vptr      = CONST_BLOCK(v, 0)->w64;
  const block_t* Ablock = CONST_BLOCK(A, 0);

  word256 cval[3] = {mm256_zero, mm256_zero, mm256_zero};
  for (unsigned int w = 2; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int i = sizeof(word) * 8; i; i -= 1, idx >>= 1) {
      const word256 mask = mm256_compute_mask(idx, 0);
      for (unsigned int j = 0; j < 3; ++j, ++Ablock) {
        cval[j] = mm256_xor_mask(cval[j], mm256_load(Ablock->w64), mask);
      }
    }
  }

  for (unsigned int j = 0; j < 3; ++j) {
    mm256_store(BLOCK(c, j)->w64, cval[j]);
  }
}
#endif

#if defined(WITH_LOWMC_192_192_30)
ATTR_TARGET_AVX2
void mzd_mul_v_s256_192_1024(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  const word* vptr      = CONST_BLOCK(v, 0)->w64;
  const block_t* Ablock = CONST_BLOCK(A, 0);

  word256 cval[4] = {mm256_zero, mm256_zero, mm256_zero, mm256_zero};
  for (unsigned int w = 3; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int i = sizeof(word) * 8; i; i -= 1, idx >>= 1) {
      const word256 mask = mm256_compute_mask(idx, 0);
      for (unsigned int j = 0; j < 4; ++j, ++Ablock) {
        cval[j] = mm256_xor_mask(cval[j], mm256_load(Ablock->w64), mask);
      }
    }
  }

  for (unsigned int j = 0; j < 4; ++j) {
    mm256_store(BLOCK(c, j)->w64, cval[j]);
  }
}
#endif

#if defined(WITH_LOWMC_256_256_38)
ATTR_TARGET_AVX2
void mzd_mul_v_s256_256_1280(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  const word* vptr      = CONST_BLOCK(v, 0)->w64;
  const block_t* Ablock = CONST_BLOCK(A, 0);

  word256 cval[5] = {mm256_zero, mm256_zero, mm256_zero, mm256_zero, mm256_zero};
  for (unsigned int w = 4; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int i = sizeof(word) * 8; i; i -= 1, idx >>= 1) {
      const word256 mask = mm256_compute_mask(idx, 0);
      for (unsigned int j = 0; j < 5; ++j, ++Ablock) {
        cval[j] = mm256_xor_mask(cval[j], mm256_load(Ablock->w64), mask);
      }
    }
  }

  for (unsigned int j = 0; j < 5; ++j) {
    mm256_store(BLOCK(c, j)->w64, cval[j]);
  }
}
#endif
#endif
#endif

static inline void clear_uint64_block(block_t* block, const unsigned int idx) {
  for (unsigned int i = 0; i < idx; ++i) {
    block->w64[i] = 0;
  }
}

// clang-format off
#if defined(WITH_LOWMC_128_128_20) || defined(WITH_LOWMC_192_192_30) || defined(WITH_LOWMC_256_256_38)
// clang-format on
static inline void clear_uint64_blocks(block_t* block, unsigned int len) {
  for (; len; --len, ++block) {
    clear_uint64_block(block, 4);
  }
}
#endif

static inline void mzd_xor_mask_uint64_block(block_t* rblock, const block_t* fblock,
                                             const word mask, const unsigned int idx) {
  for (unsigned int i = 0; i < idx; ++i) {
    rblock->w64[i] ^= fblock->w64[i] & mask;
  }
}

void mzd_addmul_v_uint64_128(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  block_t* cblock       = BLOCK(c, 0);
  const word* vptr      = CONST_BLOCK(v, 0)->w64;
  const block_t* Ablock = CONST_BLOCK(A, 0);

  for (unsigned int w = 2; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int i = sizeof(word) * 8; i; i -= 2, idx >>= 2, Ablock += 1) {
      const word mask1 = -(idx & 1);
      const word mask2 = -((idx >> 1) & 1);
      cblock->w64[0] ^= (Ablock->w64[0] & mask1) ^ (Ablock->w64[2] & mask2);
      cblock->w64[1] ^= (Ablock->w64[1] & mask1) ^ (Ablock->w64[3] & mask2);
    }
  }
}

void mzd_mul_v_uint64_128(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  clear_uint64_block(BLOCK(c, 0), 2);
  mzd_addmul_v_uint64_128(c, v, A);
}

void mzd_addmul_v_uint64_129(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  block_t* cblock       = BLOCK(c, 0);
  const word* vptr      = CONST_BLOCK(v, 0)->w64;
  const block_t* Ablock = CONST_BLOCK(A, 0);

  Ablock += 63;
  {
    word idx        = (*vptr) >> 63;
    const word mask = -(idx & 1);
    mzd_xor_mask_uint64_block(cblock, Ablock, mask, 3);
    Ablock++;
    vptr++;
  }

  for (unsigned int w = 2; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int i = sizeof(word) * 8; i; --i, idx >>= 1, ++Ablock) {
      const word mask = -(idx & 1);
      mzd_xor_mask_uint64_block(cblock, Ablock, mask, 3);
    }
  }
}

void mzd_mul_v_uint64_129(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  clear_uint64_block(BLOCK(c, 0), 3);
  mzd_addmul_v_uint64_129(c, v, A);
}

void mzd_addmul_v_uint64_192(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  block_t* cblock       = BLOCK(c, 0);
  const word* vptr      = CONST_BLOCK(v, 0)->w64;
  const block_t* Ablock = CONST_BLOCK(A, 0);

  for (unsigned int w = 3; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int i = sizeof(word) * 8; i; --i, idx >>= 1, ++Ablock) {
      const word mask = -(idx & 1);
      mzd_xor_mask_uint64_block(cblock, Ablock, mask, 3);
    }
  }
}

void mzd_mul_v_uint64_192(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  clear_uint64_block(BLOCK(c, 0), 3);
  mzd_addmul_v_uint64_192(c, v, A);
}

void mzd_addmul_v_uint64_256(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  block_t* cblock       = BLOCK(c, 0);
  const word* vptr      = CONST_BLOCK(v, 0)->w64;
  const block_t* Ablock = CONST_BLOCK(A, 0);

  for (unsigned int w = 4; w; --w, ++vptr) {
    word idx = *vptr;

    for (unsigned int i = sizeof(word) * 8; i; --i, idx >>= 1, ++Ablock) {
      const word mask = -(idx & 1);
      mzd_xor_mask_uint64_block(cblock, Ablock, mask, 4);
    }
  }
}

void mzd_mul_v_uint64_256(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  clear_uint64_block(BLOCK(c, 0), 4);
  mzd_addmul_v_uint64_256(c, v, A);
}

#if defined(WITH_LOWMC_128_128_20)
void mzd_mul_v_uint64_128_640(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  const word* vptr      = CONST_BLOCK(v, 0)->w64;
  const block_t* Ablock = CONST_BLOCK(A, 0);

  clear_uint64_blocks(BLOCK(c, 0), 2);
  clear_uint64_block(BLOCK(c, 2), 2);

  for (unsigned int w = 2; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int i = sizeof(word) * 8; i; --i, idx >>= 1, ++Ablock) {
      const word mask = -(idx & 1);
      for (unsigned int j = 0; j < 2; ++j, ++Ablock) {
        mzd_xor_mask_uint64_block(BLOCK(c, j), Ablock, mask, 4);
      }
      mzd_xor_mask_uint64_block(BLOCK(c, 2), Ablock, mask, 2);
    }
  }
}
#endif

#if defined(WITH_LOWMC_192_192_30)
void mzd_mul_v_uint64_192_960(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  const word* vptr      = CONST_BLOCK(v, 0)->w64;
  const block_t* Ablock = CONST_BLOCK(A, 0);

  clear_uint64_blocks(BLOCK(c, 0), 3);
  clear_uint64_block(BLOCK(c, 3), 3);

  for (unsigned int w = 3; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int i = sizeof(word) * 8; i; --i, idx >>= 1, ++Ablock) {
      const word mask = -(idx & 1);
      for (unsigned int j = 0; j < 3; ++j, ++Ablock) {
        mzd_xor_mask_uint64_block(BLOCK(c, j), Ablock, mask, 4);
      }
      mzd_xor_mask_uint64_block(BLOCK(c, 3), Ablock, mask, 3);
    }
  }
}
#endif

#if defined(WITH_LOWMC_256_256_38)
void mzd_mul_v_uint64_256_1216(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  const word* vptr      = CONST_BLOCK(v, 0)->w64;
  const block_t* Ablock = CONST_BLOCK(A, 0);

  clear_uint64_blocks(BLOCK(c, 0), 4);
  clear_uint64_block(BLOCK(c, 4), 3);

  for (unsigned int w = 4; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int i = sizeof(word) * 8; i; --i, idx >>= 1, ++Ablock) {
      const word mask = -(idx & 1);
      for (unsigned int j = 0; j < 4; ++j, ++Ablock) {
        mzd_xor_mask_uint64_block(BLOCK(c, j), Ablock, mask, 4);
      }
      mzd_xor_mask_uint64_block(BLOCK(c, 4), Ablock, mask, 3);
    }
  }
}
#endif

// clang-format off
#if defined(WITH_LOWMC_128_128_20) || defined(WITH_LOWMC_192_192_30) || defined(WITH_LOWMC_256_256_38)
// clang-format on
// specific instances
// bit extract, non-constant time for mask, but mask is public in our calls
static word extract_bits(word in, word mask) {
  word res = 0;
  for (word bb = 1; mask != 0; bb <<= 1, mask &= (mask - 1)) {
    res |= bb & (-((word) !!(in & mask & -mask)));
  }
  return res;
}

static inline void mzd_shuffle_30_idx(mzd_local_t* x, const word mask, unsigned int idx) {
  const word w          = CONST_BLOCK(x, 0)->w64[idx];
  const word a          = extract_bits(w, mask) << 34;
  BLOCK(x, 0)->w64[idx] = a | extract_bits(w, ~mask);
}
#endif

#if defined(WITH_LOWMC_128_128_20)
void mzd_shuffle_128_30(mzd_local_t* x, const word mask) {
  mzd_shuffle_30_idx(x, mask, 1);
}
#endif

#if defined(WITH_LOWMC_192_192_30)
void mzd_shuffle_192_30(mzd_local_t* x, const word mask) {
  mzd_shuffle_30_idx(x, mask, 2);
}
#endif

#if defined(WITH_LOWMC_256_256_38)
void mzd_shuffle_256_30(mzd_local_t* x, const word mask) {
  mzd_shuffle_30_idx(x, mask, 3);
}
#endif

// no SIMD
#if defined(WITH_LOWMC_128_128_20)
void mzd_addmul_v_uint64_30_128(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  block_t* cblock       = BLOCK(c, 0);
  const block_t* Ablock = CONST_BLOCK(A, 0);

  word idx = CONST_BLOCK(v, 0)->w64[1] >> 34;
  for (unsigned int i = 15; i; --i, idx >>= 2, ++Ablock) {
    const word mask1 = -(idx & 1);
    const word mask2 = -((idx >> 1) & 1);
    for (unsigned int j = 0; j < 2; ++j) {
      cblock->w64[j] ^= (Ablock->w64[j] & mask1) ^ (Ablock->w64[j + 2] & mask2);
    }
  }
}
#endif

#if defined(WITH_LOWMC_192_192_30)
void mzd_addmul_v_uint64_30_192(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  block_t* cblock       = BLOCK(c, 0);
  const block_t* Ablock = CONST_BLOCK(A, 0);

  word idx = CONST_BLOCK(v, 0)->w64[2] >> 34;
  for (unsigned int i = 30; i; --i, idx >>= 1, ++Ablock) {
    const word mask = -(idx & 1);
    mzd_xor_mask_uint64_block(cblock, Ablock, mask, 3);
  }
}
#endif

#if defined(WITH_LOWMC_256_256_38)
void mzd_addmul_v_uint64_30_256(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  block_t* cblock       = BLOCK(c, 0);
  const block_t* Ablock = CONST_BLOCK(A, 0);

  word idx = CONST_BLOCK(v, 0)->w64[3] >> 34;
  for (unsigned int i = 30; i; --i, idx >>= 1, ++Ablock) {
    const word mask = -(idx & 1);
    mzd_xor_mask_uint64_block(cblock, Ablock, mask, 4);
  }
}
#endif

#if defined(WITH_SSE2) || defined(WITH_NEON)
#if defined(WITH_LOWMC_128_128_20)
ATTR_TARGET_S128
void mzd_addmul_v_s128_30_128(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  block_t* cblock       = BLOCK(c, 0);
  const block_t* Ablock = CONST_BLOCK(A, 0);
  word idx              = CONST_BLOCK(v, 0)->w64[1] >> 34;

  word128 cval[2];
  cval[0] = mm128_xor_mask(mm128_load(cblock->w64), mm128_load(&Ablock[0].w64[0]),
                           mm128_compute_mask(idx, 0));
  cval[1] = mm128_and(mm128_load(&Ablock[0].w64[2]), mm128_compute_mask(idx, 1));
  idx >>= 2;
  Ablock += 1;

  for (unsigned int i = 28; i; i -= 4, idx >>= 4, Ablock += 2) {
    cval[0] = mm128_xor_mask(cval[0], mm128_load(&Ablock[0].w64[0]), mm128_compute_mask(idx, 0));
    cval[1] = mm128_xor_mask(cval[1], mm128_load(&Ablock[0].w64[2]), mm128_compute_mask(idx, 1));
    cval[0] = mm128_xor_mask(cval[0], mm128_load(&Ablock[1].w64[0]), mm128_compute_mask(idx, 2));
    cval[1] = mm128_xor_mask(cval[1], mm128_load(&Ablock[1].w64[2]), mm128_compute_mask(idx, 3));
  }
  mm128_store(cblock->w64, mm128_xor(cval[0], cval[1]));
}
#endif

#if defined(WITH_LOWMC_192_192_30) || defined(WITH_LOWMC_256_256_38)
ATTR_TARGET_S128
static inline void mzd_addmul_v_s128_30_256_idx(mzd_local_t* c, mzd_local_t const* A, word idx) {
  block_t* cblock       = BLOCK(c, 0);
  const block_t* Ablock = CONST_BLOCK(A, 0);

  word128 cval[4] = {mm128_load(&cblock->w64[0]), mm128_load(&cblock->w64[2]), mm128_zero,
                     mm128_zero};
  for (unsigned int i = 30; i; i -= 2, idx >>= 2, Ablock += 2) {
    mm128_xor_mask_region(&cval[0], Ablock[0].w64, mm128_compute_mask(idx, 0), 2);
    mm128_xor_mask_region(&cval[2], Ablock[1].w64, mm128_compute_mask(idx, 1), 2);
  }
  mm128_store(&cblock->w64[0], mm128_xor(cval[0], cval[2]));
  mm128_store(&cblock->w64[2], mm128_xor(cval[1], cval[3]));
}
#endif

#if defined(WITH_LOWMC_192_192_30)
ATTR_TARGET_S128
void mzd_addmul_v_s128_30_192(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  mzd_addmul_v_s128_30_256_idx(c, A, CONST_BLOCK(v, 0)->w64[2] >> 34);
}
#endif

#if defined(WITH_LOWMC_256_256_38)
ATTR_TARGET_S128
void mzd_addmul_v_s128_30_256(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  mzd_addmul_v_s128_30_256_idx(c, A, CONST_BLOCK(v, 0)->w64[3] >> 34);
}
#endif
#endif

#if defined(WITH_AVX2)
#if defined(WITH_LOWMC_128_128_20)
ATTR_TARGET_AVX2
void mzd_addmul_v_s256_30_128(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  block_t* cblock       = BLOCK(c, 0);
  const block_t* vblock = CONST_BLOCK(v, 0);
  const block_t* Ablock = CONST_BLOCK(A, 0);
  word idx              = vblock->w64[1] >> 34;

  word256 cval[2];
  cval[0] = mm256_xor_mask(_mm256_setr_m128i(mm128_load(cblock->w64), mm128_zero),
                           mm256_load(Ablock[0].w64), mm256_compute_mask_2(idx, 0));
  cval[1] = mm256_and(mm256_load(Ablock[1].w64), mm256_compute_mask_2(idx, 2));
  idx >>= 4;
  Ablock += 2;

  for (unsigned int i = 24; i; i -= 8, idx >>= 8, Ablock += 4) {
    cval[0] = mm256_xor_mask(cval[0], mm256_load(Ablock[0].w64), mm256_compute_mask_2(idx, 0));
    cval[1] = mm256_xor_mask(cval[1], mm256_load(Ablock[1].w64), mm256_compute_mask_2(idx, 2));
    cval[0] = mm256_xor_mask(cval[0], mm256_load(Ablock[2].w64), mm256_compute_mask_2(idx, 4));
    cval[1] = mm256_xor_mask(cval[1], mm256_load(Ablock[3].w64), mm256_compute_mask_2(idx, 6));
  }
  cval[0] = mm256_xor_mask(cval[0], mm256_load(Ablock[0].w64), mm256_compute_mask_2(idx, 0));

  cval[0] = mm256_xor(cval[0], cval[1]);
  mm128_store(cblock->w64, mm128_xor(_mm256_extractf128_si256(cval[0], 0),
                                     _mm256_extractf128_si256(cval[0], 1)));
}
#endif

#if defined(WITH_LOWMC_192_192_30) || defined(WITH_LOWMC_256_256_38)
ATTR_TARGET_AVX2
static inline void mzd_addmul_v_s256_30_256_idx(mzd_local_t* c, mzd_local_t const* A, word idx) {
  block_t* cblock       = BLOCK(c, 0);
  const block_t* Ablock = CONST_BLOCK(A, 0);

  word256 cval[2];
  cval[0] = mm256_xor_mask(mm256_load(cblock->w64), mm256_load(Ablock[0].w64),
                           mm256_compute_mask(idx, 0));
  cval[1] = mm256_and(mm256_load(Ablock[1].w64), mm256_compute_mask(idx, 1));
  idx >>= 2;
  Ablock += 2;

  for (unsigned int i = 28; i; i -= 4, idx >>= 4, Ablock += 4) {
    cval[0] = mm256_xor_mask(cval[0], mm256_load(Ablock[0].w64), mm256_compute_mask(idx, 0));
    cval[1] = mm256_xor_mask(cval[1], mm256_load(Ablock[1].w64), mm256_compute_mask(idx, 1));
    cval[0] = mm256_xor_mask(cval[0], mm256_load(Ablock[2].w64), mm256_compute_mask(idx, 2));
    cval[1] = mm256_xor_mask(cval[1], mm256_load(Ablock[3].w64), mm256_compute_mask(idx, 3));
  }
  mm256_store(cblock->w64, mm256_xor(cval[0], cval[1]));
}
#endif

#if defined(WITH_LOWMC_192_192_30)
ATTR_TARGET_AVX2
void mzd_addmul_v_s256_30_192(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  mzd_addmul_v_s256_30_256_idx(c, A, CONST_BLOCK(v, 0)->w64[2] >> 34);
}
#endif

#if defined(WITH_LOWMC_256_256_38)
ATTR_TARGET_AVX2
void mzd_addmul_v_s256_30_256(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  mzd_addmul_v_s256_30_256_idx(c, A, CONST_BLOCK(v, 0)->w64[3] >> 34);
}
#endif

// clang-format off
#if defined(WITH_LOWMC_128_128_20) || defined(WITH_LOWMC_192_192_30) || defined(WITH_LOWMC_256_256_38)
// clang-format on
#if !defined(__x86_64__) && !defined(_M_X64)
ATTR_TARGET_AVX2 ATTR_CONST static uint8_t popcount_32(uint32_t value) {
  uint64_t result =
      ((value & 0xfff) * UINT64_C(0x1001001001001) & UINT64_C(0x84210842108421)) % 0x1f;
  result +=
      (((value & 0xfff000) >> 12) * UINT64_C(0x1001001001001) & UINT64_C(0x84210842108421)) % 0x1f;
  result += ((value >> 24) * UINT64_C(0x1001001001001) & UINT64_C(0x84210842108421)) % 0x1f;
  return result;
}

ATTR_TARGET_AVX2 ATTR_CONST static uint64_t _pext_u64(uint64_t a, uint64_t mask) {
  const uint32_t low  = _pext_u32(a, mask);
  const uint32_t high = _pext_u32(a >> 32, mask >> 32);

  return (((uint64_t)high) << popcount_32(mask)) | low;
}
#endif

ATTR_TARGET_AVX2
static inline void mzd_shuffle_pext_30_idx(mzd_local_t* x, const word mask, unsigned int idx) {
  const word w          = CONST_BLOCK(x, 0)->w64[idx];
  const word a          = _pext_u64(w, mask) << 34;
  BLOCK(x, 0)->w64[idx] = a | _pext_u64(w, ~mask);
}
#endif

#if defined(WITH_LOWMC_128_128_20)
ATTR_TARGET_AVX2
void mzd_shuffle_pext_128_30(mzd_local_t* x, const word mask) {
  mzd_shuffle_pext_30_idx(x, mask, 1);
}
#endif

#if defined(WITH_LOWMC_192_192_30)
ATTR_TARGET_AVX2
void mzd_shuffle_pext_192_30(mzd_local_t* x, const word mask) {
  mzd_shuffle_pext_30_idx(x, mask, 2);
}
#endif

#if defined(WITH_LOWMC_256_256_38)
ATTR_TARGET_AVX2
void mzd_shuffle_pext_256_30(mzd_local_t* x, const word mask) {
  mzd_shuffle_pext_30_idx(x, mask, 3);
}

/**
 * Compute the parity of three 256-bit words and store them in the lowest bit (v1), the second
 * lowest bit (v2), and the third lowest bit (v3) of each quad word
 *
 * Based on Wojciech Mula, Nathan Kurz and Daniel Lemire: Faster Population Counts Using AVX2
 * Instructions
 */
ATTR_TARGET_AVX2
static inline word256 mm256_parity_3(word256 v1, word256 v2, word256 v3) {
  const word256 lookup   = _mm256_setr_epi8(0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 0, 1, 1,
                                            2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4);
  const word256 low_mask = _mm256_set1_epi8(0x0f);
  const word256 all1     = _mm256_set1_epi64x(0x1);

  const word256 lo1     = mm256_and(v1, low_mask);
  const word256 hi1     = mm256_and(_mm256_srli_epi32(v1, 4), low_mask);
  const word256 popcnt1 = _mm256_shuffle_epi8(lookup, lo1);
  const word256 popcnt2 = _mm256_shuffle_epi8(lookup, hi1);

  const word256 lo2     = mm256_and(v2, low_mask);
  const word256 hi2     = mm256_and(_mm256_srli_epi32(v2, 4), low_mask);
  const word256 popcnt3 = _mm256_shuffle_epi8(lookup, lo2);
  const word256 popcnt4 = _mm256_shuffle_epi8(lookup, hi2);

  const word256 lo3     = mm256_and(v3, low_mask);
  const word256 hi3     = mm256_and(_mm256_srli_epi32(v3, 4), low_mask);
  const word256 popcnt5 = _mm256_shuffle_epi8(lookup, lo3);
  const word256 popcnt6 = _mm256_shuffle_epi8(lookup, hi3);

  const word256 total1 = mm256_and(_mm256_sad_epu8(popcnt1, popcnt2), all1);
  const word256 total2 = mm256_and(_mm256_sad_epu8(popcnt3, popcnt4), all1);
  const word256 total3 = mm256_and(_mm256_sad_epu8(popcnt5, popcnt6), all1);

  return mm256_xor(mm256_xor(total1, _mm256_slli_epi64(total2, 1)), _mm256_slli_epi64(total3, 2));
}

ATTR_TARGET_AVX2
void mzd_mul_v_parity_s256_256_30(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* At) {
  const word256 vblock = mm256_load(v->w64);
  word256 res          = mm256_zero;
  // process 3 rows at a time
  for (unsigned int i = 30; i; i -= 3) {
    const word256 Ablock1 = mm256_and(vblock, mm256_load(CONST_BLOCK(At, 30 - i)->w64));
    const word256 Ablock2 = mm256_and(vblock, mm256_load(CONST_BLOCK(At, 30 - i + 1)->w64));
    const word256 Ablock3 = mm256_and(vblock, mm256_load(CONST_BLOCK(At, 30 - i + 2)->w64));

    res = mm256_xor(
        res, _mm256_sll_epi64(mm256_parity_3(Ablock1, Ablock2, Ablock3), _mm_set1_epi64x(64 - i)));
  }

  for (unsigned int j = 0; j < 3; j++) {
    c->w64[j] = 0;
  }

  // combine the parities of each quad word
#if defined(__x86_64__) || defined(_M_X64)
  c->w64[3] = _mm256_extract_epi64(res, 0) ^ _mm256_extract_epi64(res, 1) ^
              _mm256_extract_epi64(res, 2) ^ _mm256_extract_epi64(res, 3);
#else
  // workaround for missing _mm256_extract_epi64 on x86-32
  uint64_t tmp[4] ATTR_ALIGNED(sizeof(word256));
  mm256_store(tmp, res);
  c->w64[3] = tmp[0] ^ tmp[1] ^ tmp[2] ^ tmp[3];
#endif
}
#endif
#endif
