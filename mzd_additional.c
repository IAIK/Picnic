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

#include "compat.h"
#include "mzd_additional.h"

#if !defined(_MSC_VER)
#include <stdalign.h>
#endif
#include <assert.h>
#include <stdlib.h>
#include <string.h>

#if !defined(_MSC_VER) && !defined(static_assert)
#define static_assert _Static_assert
#endif

static const size_t mzd_local_t_size = (sizeof(mzd_local_t) + 0x1f) & ~0x1f;
static_assert(((sizeof(mzd_local_t) + 0x1f) & ~0x1f) == 32, "sizeof mzd_local_t not supported");


#if defined(WITH_OPT)
#include "simd.h"

#if defined(WITH_SSE2) || defined(WITH_AVX2) || defined(WITH_NEON)
static const unsigned int word_size_bits = 8 * sizeof(word);
#endif
#endif
static const unsigned int align_bound = 128 / (8 * sizeof(word));

static uint32_t calculate_rowstride(uint32_t width) {
  // As soon as we hit the AVX bound, use 32 byte alignment. Otherwise use 16
  // byte alignment for SSE2 and 128 bit vectors.
  if (width > align_bound) {
    return ((width * sizeof(word) + 31) & ~31) / sizeof(word);
  } else {
    return ((width * sizeof(word) + 15) & ~15) / sizeof(word);
  }
}

static uint32_t calculate_width(uint32_t c) {
  return (c + sizeof(word) * 8 - 1) / (sizeof(word) * 8);
}

// Notes on the memory layout: mzd_init allocates multiple memory blocks (one
// for mzd_local_t, one for rows and multiple for the buffers). We use one memory
// block for mzd_local_t, rows and the buffer. This improves memory locality and
// requires less calls to malloc.
//
// In mzd_local_init_multiple we do the same, but store n mzd_local_t instances in one
// memory block.

mzd_local_t* mzd_local_init_ex(uint32_t r, uint32_t c, bool clear) {
  const uint32_t width     = calculate_width(c);
  const uint32_t rowstride = calculate_rowstride(width);

  const size_t buffer_size = r * rowstride * sizeof(word);

  /* We always align mzd_local_ts to 32 bytes. Thus the first row is always
   * aligned to 32 bytes as well. For 128 bit and SSE all other rows are then
   * aligned to 16 bytes. */
  unsigned char* buffer = aligned_alloc(32, (mzd_local_t_size + buffer_size + 31) & ~31);

  mzd_local_t* A = (mzd_local_t*)buffer;
  buffer += mzd_local_t_size;

  if (clear) {
    memset(buffer, 0, buffer_size);
  }

  // assign in order
  A->nrows     = r;
  A->ncols     = c;
  A->width     = width;
  A->rowstride = rowstride;

  return A;
}

void mzd_local_free(mzd_local_t* v) {
  aligned_free(v);
}

void mzd_local_init_multiple_ex(mzd_local_t** dst, size_t n, uint32_t r, uint32_t c, bool clear) {
  const uint32_t width     = calculate_width(c);
  const uint32_t rowstride = calculate_rowstride(width);

  const size_t buffer_size   = r * rowstride * sizeof(word);
  const size_t size_per_elem = (mzd_local_t_size + buffer_size + 31) & ~31;

  unsigned char* full_buffer = aligned_alloc(32, size_per_elem * n);

  for (size_t s = 0; s < n; ++s, full_buffer += size_per_elem) {
    unsigned char* buffer = full_buffer;
    mzd_local_t* A        = (mzd_local_t*)buffer;
    dst[s]                = A;

    buffer += mzd_local_t_size;

    if (clear) {
      memset(buffer, 0, buffer_size);
    }

    // assign in order
    A->nrows     = r;
    A->ncols     = c;
    A->width     = width;
    A->rowstride = rowstride;
  }
}

void mzd_local_free_multiple(mzd_local_t** vs) {
  if (vs) {
    aligned_free(vs[0]);
  }
}

mzd_local_t* mzd_local_copy(mzd_local_t* dst, mzd_local_t const* src) {
  if (dst == src) {
    return dst;
  }

  if (!dst) {
    dst = mzd_local_init(src->nrows, src->ncols);
  }

  memcpy(ASSUME_ALIGNED(FIRST_ROW(dst), 32), ASSUME_ALIGNED(CONST_FIRST_ROW(src), 32),
         src->nrows * sizeof(word) * src->rowstride);
  return dst;
}

void mzd_local_clear(mzd_local_t* c) {
  memset(ASSUME_ALIGNED(FIRST_ROW(c), 32), 0, c->nrows * sizeof(word) * c->rowstride);
}

void mzd_shift_right(mzd_local_t* res, mzd_local_t const* val, unsigned count) {
  if (!count) {
    mzd_local_copy(res, val);
    return;
  }

  const unsigned int nwords     = val->width;
  const unsigned int left_count = 8 * sizeof(word) - count;

  word* resptr       = ASSUME_ALIGNED(FIRST_ROW(res), 32);
  word const* valptr = ASSUME_ALIGNED(CONST_FIRST_ROW(val), 32);

  for (unsigned int i = nwords - 1; i; --i, ++resptr) {
    const word tmp = *valptr >> count;
    *resptr        = tmp | (*++valptr << left_count);
  }
  *resptr = *valptr >> count;
}

void mzd_shift_left(mzd_local_t* res, mzd_local_t const* val, unsigned count) {
  if (!count) {
    mzd_local_copy(res, val);
    return;
  }

  const unsigned int nwords      = val->width;
  const unsigned int right_count = 8 * sizeof(word) - count;

  word* resptr       = FIRST_ROW(res) + nwords - 1;
  word const* valptr = CONST_FIRST_ROW(val) + nwords - 1;

  for (unsigned int i = nwords - 1; i; --i, --resptr) {
    const word tmp = *valptr << count;
    *resptr        = tmp | (*--valptr >> right_count);
  }
  *resptr = *valptr << count;
}

#if defined(WITH_OPT)
#if defined(WITH_SSE2)
ATTR_TARGET("sse2")
static inline void mzd_and_sse(mzd_local_t* res, mzd_local_t const* first,
                               mzd_local_t const* second) {
  unsigned int width        = first->rowstride;
  __m128i* mresptr          = ASSUME_ALIGNED(FIRST_ROW(res), alignof(__m128i));
  __m128i const* mfirstptr  = ASSUME_ALIGNED(CONST_FIRST_ROW(first), alignof(__m128i));
  __m128i const* msecondptr = ASSUME_ALIGNED(CONST_FIRST_ROW(second), alignof(__m128i));

  do {
    *mresptr++ = _mm_and_si128(*mfirstptr++, *msecondptr++);
    width -= sizeof(__m128i) / sizeof(word);
  } while (width);
}
#endif

#if defined(WITH_AVX2)
ATTR_TARGET("avx2")
static inline void mzd_and_avx(mzd_local_t* res, mzd_local_t const* first,
                               mzd_local_t const* second) {
  unsigned int width        = first->rowstride;
  __m256i* mresptr          = ASSUME_ALIGNED(FIRST_ROW(res), alignof(__m256i));
  __m256i const* mfirstptr  = ASSUME_ALIGNED(CONST_FIRST_ROW(first), alignof(__m256i));
  __m256i const* msecondptr = ASSUME_ALIGNED(CONST_FIRST_ROW(second), alignof(__m256i));

  do {
    *mresptr++ = _mm256_and_si256(*mfirstptr++, *msecondptr++);
    width -= sizeof(__m256i) / sizeof(word);
  } while (width);
}
#endif

#if defined(WITH_NEON)
static inline void mzd_and_neon(mzd_local_t* res, mzd_local_t const* first,
                                mzd_local_t const* second) {
  unsigned int width           = first->rowstride;
  uint32x4_t* mresptr          = ASSUME_ALIGNED(FIRST_ROW(res), alignof(uint32x4_t));
  uint32x4_t const* mfirstptr  = ASSUME_ALIGNED(CONST_FIRST_ROW(first), alignof(uint32x4_t));
  uint32x4_t const* msecondptr = ASSUME_ALIGNED(CONST_FIRST_ROW(second), alignof(uint32x4_t));

  do {
    *mresptr++ = vandq_u32(*mfirstptr++, *msecondptr++);
    width -= sizeof(uint32x4_t) / sizeof(word);
  } while (width);
}
#endif
#endif

static inline void mzd_and_uint64(mzd_local_t* res, mzd_local_t const* first,
                                  mzd_local_t const* second) {
  unsigned int width    = first->width;
  word* resptr          = ASSUME_ALIGNED(FIRST_ROW(res), 32);
  word const* firstptr  = ASSUME_ALIGNED(CONST_FIRST_ROW(first), 32);
  word const* secondptr = ASSUME_ALIGNED(CONST_FIRST_ROW(second), 32);

  while (width--) {
    *resptr++ = *firstptr++ & *secondptr++;
  }
}

void mzd_and(mzd_local_t* res, mzd_local_t const* first, mzd_local_t const* second) {
#if defined(WITH_OPT)
#if defined(WITH_AVX2)
  if (CPU_SUPPORTS_AVX2 && first->ncols >= 256 && ((first->ncols & (word_size_bits - 1)) == 0)) {
    mzd_and_avx(res, first, second);
    return;
  }
#endif
#if defined(WITH_SSE2)
  if (CPU_SUPPORTS_SSE2 && ((first->ncols & (word_size_bits - 1)) == 0)) {
    mzd_and_sse(res, first, second);
    return;
  }
#endif
#if defined(WITH_NEON)
  if (CPU_SUPPORTS_NEON && first->ncols % ((first->ncols & (word_size_bits - 1)) == 0)) {
    mzd_and_neon(res, first, second);
    return;
  }
#endif
#endif
  mzd_and_uint64(res, first, second);
}

#if defined(WITH_OPT)
#if defined(WITH_SSE2)
ATTR_TARGET("sse2")
void mzd_xor_sse(mzd_local_t* res, mzd_local_t const* first, mzd_local_t const* second) {
  unsigned int width        = first->rowstride;
  __m128i* mresptr          = ASSUME_ALIGNED(FIRST_ROW(res), alignof(__m128i));
  __m128i const* mfirstptr  = ASSUME_ALIGNED(CONST_FIRST_ROW(first), alignof(__m128i));
  __m128i const* msecondptr = ASSUME_ALIGNED(CONST_FIRST_ROW(second), alignof(__m128i));

  do {
    *mresptr++ = _mm_xor_si128(*mfirstptr++, *msecondptr++);
    width -= sizeof(__m128i) / sizeof(word);
  } while (width);
}

ATTR_TARGET("sse2")
void mzd_xor_sse_128(mzd_local_t* res, mzd_local_t const* first, mzd_local_t const* second) {
  __m128i* mresptr          = ASSUME_ALIGNED(FIRST_ROW(res), alignof(__m128i));
  __m128i const* mfirstptr  = ASSUME_ALIGNED(CONST_FIRST_ROW(first), alignof(__m128i));
  __m128i const* msecondptr = ASSUME_ALIGNED(CONST_FIRST_ROW(second), alignof(__m128i));

  *mresptr = _mm_xor_si128(*mfirstptr, *msecondptr);
}

ATTR_TARGET("sse2")
void mzd_xor_sse_256(mzd_local_t* res, mzd_local_t const* first, mzd_local_t const* second) {
  __m128i* mresptr          = ASSUME_ALIGNED(FIRST_ROW(res), alignof(__m128i));
  __m128i const* mfirstptr  = ASSUME_ALIGNED(CONST_FIRST_ROW(first), alignof(__m128i));
  __m128i const* msecondptr = ASSUME_ALIGNED(CONST_FIRST_ROW(second), alignof(__m128i));

  mresptr[0] = _mm_xor_si128(mfirstptr[0], msecondptr[0]);
  mresptr[1] = _mm_xor_si128(mfirstptr[1], msecondptr[1]);
}
#endif

#if defined(WITH_AVX2)
ATTR_TARGET("avx2")
void mzd_xor_avx(mzd_local_t* res, mzd_local_t const* first, mzd_local_t const* second) {
  unsigned int width        = first->rowstride;
  __m256i* mresptr          = ASSUME_ALIGNED(FIRST_ROW(res), alignof(__m256i));
  __m256i const* mfirstptr  = ASSUME_ALIGNED(CONST_FIRST_ROW(first), alignof(__m256i));
  __m256i const* msecondptr = ASSUME_ALIGNED(CONST_FIRST_ROW(second), alignof(__m256i));
  do {
    *mresptr++ = _mm256_xor_si256(*mfirstptr++, *msecondptr++);
    width -= sizeof(__m256i) / sizeof(word);
  } while (width);
}

ATTR_TARGET("avx2")
void mzd_xor_avx_256(mzd_local_t* res, mzd_local_t const* first, mzd_local_t const* second) {
  __m256i* mresptr          = ASSUME_ALIGNED(FIRST_ROW(res), alignof(__m256i));
  __m256i const* mfirstptr  = ASSUME_ALIGNED(CONST_FIRST_ROW(first), alignof(__m256i));
  __m256i const* msecondptr = ASSUME_ALIGNED(CONST_FIRST_ROW(second), alignof(__m256i));

  *mresptr = _mm256_xor_si256(*mfirstptr, *msecondptr);
}
#endif

#if defined(WITH_NEON)
void mzd_xor_neon(mzd_local_t* res, mzd_local_t const* first, mzd_local_t const* second) {
  unsigned int width           = first->rowstride;
  uint32x4_t* mresptr          = ASSUME_ALIGNED(FIRST_ROW(res), alignof(uint32x4_t));
  uint32x4_t const* mfirstptr  = ASSUME_ALIGNED(CONST_FIRST_ROW(first), alignof(uint32x4_t));
  uint32x4_t const* msecondptr = ASSUME_ALIGNED(CONST_FIRST_ROW(second), alignof(uint32x4_t));

  do {
    *mresptr++ = veorq_u32(*mfirstptr++, *msecondptr++);
    width -= sizeof(uint32x4_t) / sizeof(word);
  } while (width);
}

void mzd_xor_neon_128(mzd_local_t* res, mzd_local_t const* first, mzd_local_t const* second) {
  uint32x4_t* mresptr          = ASSUME_ALIGNED(FIRST_ROW(res), alignof(uint32x4_t));
  uint32x4_t const* mfirstptr  = ASSUME_ALIGNED(CONST_FIRST_ROW(first), alignof(uint32x4_t));
  uint32x4_t const* msecondptr = ASSUME_ALIGNED(CONST_FIRST_ROW(second), alignof(uint32x4_t));

  *mresptr = veorq_u32(*mfirstptr, *msecondptr);
}

void mzd_xor_neon_256(mzd_local_t* res, mzd_local_t const* first, mzd_local_t const* second) {
  uint32x4_t* mresptr          = ASSUME_ALIGNED(FIRST_ROW(res), alignof(uint32x4_t));
  uint32x4_t const* mfirstptr  = ASSUME_ALIGNED(CONST_FIRST_ROW(first), alignof(uint32x4_t));
  uint32x4_t const* msecondptr = ASSUME_ALIGNED(CONST_FIRST_ROW(second), alignof(uint32x4_t));

  mresptr[0] = veorq_u32(mfirstptr[0], msecondptr[0]);
  mresptr[1] = veorq_u32(mfirstptr[1], msecondptr[1]);
}
#endif
#endif

void mzd_xor(mzd_local_t* res, mzd_local_t const* first, mzd_local_t const* second) {
#if defined(WITH_OPT)
#if defined(WITH_AVX2)
  if (CPU_SUPPORTS_AVX2 && first->ncols >= 256 && ((first->ncols & (word_size_bits - 1)) == 0)) {
    mzd_xor_avx(res, first, second);
    return;
  }
#endif
#if defined(WITH_SSE2)
  if (CPU_SUPPORTS_SSE2 && ((first->ncols & (word_size_bits - 1)) == 0)) {
    mzd_xor_sse(res, first, second);
    return;
  }
#endif
#if defined(WITH_NEON)
  if (CPU_SUPPORTS_NEON && ((first->ncols & (word_size_bits - 1)) == 0)) {
    mzd_xor_neon(res, first, second);
    return;
  }
#endif
#endif
  mzd_xor_uint64(res, first, second);
}

void mzd_xor_uint64(mzd_local_t* res, mzd_local_t const* first, mzd_local_t const* second) {
  unsigned int width    = first->width;
  word* resptr          = ASSUME_ALIGNED(FIRST_ROW(res), 32);
  word const* firstptr  = ASSUME_ALIGNED(CONST_FIRST_ROW(first), 32);
  word const* secondptr = ASSUME_ALIGNED(CONST_FIRST_ROW(second), 32);

  while (width--) {
    *resptr++ = *firstptr++ ^ *secondptr++;
  }
}

void mzd_mul_v(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* At) {
  mzd_local_clear(c);
  mzd_addmul_v(c, v, At);
}

void mzd_mul_v_uint64(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* At) {
  mzd_local_clear(c);
  mzd_addmul_v_uint64(c, v, At);
}

void mzd_mul_v_popcnt(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* At) {
//    if (At->ncols != v->ncols) {
//        // number of columns does not match
//        // since we use popcnt we want cols to be equal, A is transposed
//        exit(1);
//        return;
//    }

  mzd_local_clear(c);
  const unsigned int rowstride = At->rowstride;
  word* cptr                   = ASSUME_ALIGNED(FIRST_ROW(c), 32);
  word const* vptr             = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  const unsigned int width     = v->width;
  word const* Aptr             = ASSUME_ALIGNED(CONST_FIRST_ROW(At), 32);

  assert(c->width == 1); // only use for 32 bit or less
  for(unsigned i = c->ncols; i; --i) {
    word popcnt = 0;
    for (unsigned int w = 0; w < width; ++w) {
      word idx = vptr[w] & Aptr[w+(i-1)*rowstride];
      popcnt += __builtin_popcountll(idx);
    }
    *cptr <<= 1;
    *cptr |= (popcnt & WORD_C(0x1));
  }
}

#if defined(WITH_OPT)
#if defined(WITH_SSE2)
ATTR_TARGET("sse2") ATTR_CONST
static inline __m128i mm128_compute_mask(const word idx, const size_t bit) {
  return _mm_set1_epi64x(-((idx >> bit) & 1));
}

ATTR_TARGET("sse2")
void mzd_addmul_v_sse(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  word const* vptr              = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  const unsigned int width      = v->width;
  const unsigned int rowstride  = A->rowstride;
  const unsigned int mrowstride = rowstride * sizeof(word) / sizeof(__m128i);
  const unsigned int len        = mrowstride;

  __m128i* mcptr = ASSUME_ALIGNED(FIRST_ROW(c), alignof(__m128i));

  for (unsigned int w = 0; w < width; ++w, ++vptr) {
    word idx             = *vptr;
    __m128i const* mAptr = ASSUME_ALIGNED(CONST_ROW(A, w * sizeof(word) * 8), alignof(__m128i));

    for (unsigned int i = sizeof(word) * 8; i; --i, idx >>= 1, mAptr += mrowstride) {
      mm128_xor_mask_region(mcptr, mAptr, mm128_compute_mask(idx, 0), len);
    }
  }
}

ATTR_TARGET("sse2")
void mzd_mul_v_sse(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  mzd_local_clear(c);
  mzd_addmul_v_sse(c, v, A);
}

ATTR_TARGET("sse2")
void mzd_mul_v_sse_128(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  word const* vptr     = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  __m128i* mcptr       = ASSUME_ALIGNED(FIRST_ROW(c), alignof(__m128i));
  __m128i const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(__m128i));

  __m128i cval[2] ATTR_ALIGNED(alignof(__m128i)) = {_mm_setzero_si128(), _mm_setzero_si128()};
  for (unsigned int w = 2; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int i = sizeof(word) * 8; i; i -= 4, idx >>= 4, mAptr += 4) {
      cval[0] = mm128_xor_mask(cval[0], mAptr[0], mm128_compute_mask(idx, 0));
      cval[1] = mm128_xor_mask(cval[1], mAptr[1], mm128_compute_mask(idx, 1));
      cval[0] = mm128_xor_mask(cval[0], mAptr[2], mm128_compute_mask(idx, 2));
      cval[1] = mm128_xor_mask(cval[1], mAptr[3], mm128_compute_mask(idx, 3));
    }
  }
  *mcptr = _mm_xor_si128(cval[0], cval[1]);
}

ATTR_TARGET("sse2")
void mzd_addmul_v_sse_128(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  word const* vptr     = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  __m128i* mcptr       = ASSUME_ALIGNED(FIRST_ROW(c), alignof(__m128i));
  __m128i const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(__m128i));

  __m128i cval[2] ATTR_ALIGNED(alignof(__m128i)) = {*mcptr, _mm_setzero_si128()};
  for (unsigned int w = 2; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int i = sizeof(word) * 8; i; i -= 4, idx >>= 4, mAptr += 4) {
      cval[0] = mm128_xor_mask(cval[0], mAptr[0], mm128_compute_mask(idx, 0));
      cval[1] = mm128_xor_mask(cval[1], mAptr[1], mm128_compute_mask(idx, 1));
      cval[0] = mm128_xor_mask(cval[0], mAptr[2], mm128_compute_mask(idx, 2));
      cval[1] = mm128_xor_mask(cval[1], mAptr[3], mm128_compute_mask(idx, 3));
    }
  }
  *mcptr = _mm_xor_si128(cval[0], cval[1]);
}

ATTR_TARGET("sse2")
void mzd_mul_v_sse_192(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  word const* vptr     = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  __m128i* mcptr       = ASSUME_ALIGNED(FIRST_ROW(c), alignof(__m128i));
  __m128i const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(__m128i));

  __m128i cval[4] ATTR_ALIGNED(alignof(__m128i)) = {_mm_setzero_si128(), _mm_setzero_si128(),
                                                    _mm_setzero_si128(), _mm_setzero_si128()};
  for (unsigned int w = 3; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int i = sizeof(word) * 8; i; i -= 2, idx >>= 2, mAptr += 4) {
      mm128_xor_mask_region(&cval[0], mAptr + 0, mm128_compute_mask(idx, 0), 2);
      mm128_xor_mask_region(&cval[2], mAptr + 2, mm128_compute_mask(idx, 1), 2);
    }
  }
  mcptr[0] = _mm_xor_si128(cval[0], cval[2]);
  mcptr[1] = _mm_xor_si128(cval[1], cval[3]);
}

ATTR_TARGET("sse2")
void mzd_addmul_v_sse_192(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  word const* vptr     = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  __m128i* mcptr       = ASSUME_ALIGNED(FIRST_ROW(c), alignof(__m128i));
  __m128i const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(__m128i));

  __m128i cval[4] ATTR_ALIGNED(alignof(__m128i)) = {mcptr[0], mcptr[1], _mm_setzero_si128(),
                                                    _mm_setzero_si128()};
  for (unsigned int w = 3; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int i = sizeof(word) * 8; i; i -= 2, idx >>= 2, mAptr += 4) {
      mm128_xor_mask_region(&cval[0], mAptr + 0, mm128_compute_mask(idx, 0), 2);
      mm128_xor_mask_region(&cval[2], mAptr + 2, mm128_compute_mask(idx, 1), 2);
    }
  }
  mcptr[0] = _mm_xor_si128(cval[0], cval[2]);
  mcptr[1] = _mm_xor_si128(cval[1], cval[3]);
}

ATTR_TARGET("sse2")
void mzd_mul_v_sse_256(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  word const* vptr     = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  __m128i* mcptr       = ASSUME_ALIGNED(FIRST_ROW(c), alignof(__m128i));
  __m128i const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(__m128i));

  __m128i cval[4] ATTR_ALIGNED(alignof(__m128i)) = {_mm_setzero_si128(), _mm_setzero_si128(),
                                                    _mm_setzero_si128(), _mm_setzero_si128()};
  for (unsigned int w = 4; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int i = sizeof(word) * 8; i; i -= 2, idx >>= 2, mAptr += 4) {
      mm128_xor_mask_region(&cval[0], mAptr + 0, mm128_compute_mask(idx, 0), 2);
      mm128_xor_mask_region(&cval[2], mAptr + 2, mm128_compute_mask(idx, 1), 2);
    }
  }
  mcptr[0] = _mm_xor_si128(cval[0], cval[2]);
  mcptr[1] = _mm_xor_si128(cval[1], cval[3]);
}

ATTR_TARGET("sse2")
void mzd_addmul_v_sse_256(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  word const* vptr     = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  __m128i* mcptr       = ASSUME_ALIGNED(FIRST_ROW(c), alignof(__m128i));
  __m128i const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(__m128i));

  __m128i cval[4] ATTR_ALIGNED(alignof(__m128i)) = {mcptr[0], mcptr[1], _mm_setzero_si128(),
                                                    _mm_setzero_si128()};
  for (unsigned int w = 4; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int i = sizeof(word) * 8; i; i -= 2, idx >>= 2, mAptr += 4) {
      mm128_xor_mask_region(&cval[0], mAptr + 0, mm128_compute_mask(idx, 0), 2);
      mm128_xor_mask_region(&cval[2], mAptr + 2, mm128_compute_mask(idx, 1), 2);
    }
  }
  mcptr[0] = _mm_xor_si128(cval[0], cval[2]);
  mcptr[1] = _mm_xor_si128(cval[1], cval[3]);
}
#endif

#if defined(WITH_AVX2)
ATTR_TARGET("avx2")
ATTR_CONST static inline __m256i mm256_compute_mask(const word idx, const size_t bit) {
  return _mm256_set1_epi64x(-((idx >> bit) & 1));
}

ATTR_TARGET("avx2")
ATTR_CONST static inline __m256i mm256_compute_mask_2(const word idx, const size_t bit) {
  const uint64_t m1 = -((idx >> bit) & 1);
  const uint64_t m2 = -((idx >> (bit + 1)) & 1);
  return _mm256_set_epi64x(m2, m2, m1, m1);
}

ATTR_TARGET("avx2")
void mzd_mul_v_avx(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  mzd_local_clear(c);
  mzd_addmul_v_avx(c, v, A);
}

ATTR_TARGET("avx2")
void mzd_addmul_v_avx(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  word const* vptr              = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  const unsigned int width      = v->width;
  const unsigned int rowstride  = A->rowstride;
  const unsigned int mrowstride = rowstride * sizeof(word) / sizeof(__m256i);
  const unsigned int len        = mrowstride;

  __m256i* mcptr = ASSUME_ALIGNED(FIRST_ROW(c), alignof(__m256i));

  for (unsigned int w = 0; w < width; ++w, ++vptr) {
    word idx             = *vptr;
    __m256i const* mAptr = ASSUME_ALIGNED(CONST_ROW(A, w * sizeof(word) * 8), alignof(__m256i));

    for (unsigned int i = sizeof(word) * 8; i; --i, idx >>= 1, mAptr += mrowstride) {
      mm256_xor_mask_region(mcptr, mAptr, mm256_compute_mask(idx, 0), len);
    }
  }
}

ATTR_TARGET("avx2")
void mzd_addmul_v_avx_128(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  word const* vptr     = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  __m128i* mcptr       = ASSUME_ALIGNED(FIRST_ROW(c), alignof(__m128i));
  __m256i const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(__m256i));

  __m256i cval[2] ATTR_ALIGNED(alignof(__m256i)) = {_mm256_castsi128_si256(*mcptr),
                                                    _mm256_setzero_si256()};
  for (unsigned int w = 2; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int i = sizeof(word) * 8; i; i -= 8, idx >>= 8, mAptr += 4) {
      cval[0] = mm256_xor_mask(cval[0], mAptr[0], mm256_compute_mask_2(idx, 0));
      cval[1] = mm256_xor_mask(cval[1], mAptr[1], mm256_compute_mask_2(idx, 2));
      cval[0] = mm256_xor_mask(cval[0], mAptr[2], mm256_compute_mask_2(idx, 4));
      cval[1] = mm256_xor_mask(cval[1], mAptr[3], mm256_compute_mask_2(idx, 6));
    }
  }
  cval[0] = _mm256_xor_si256(cval[0], cval[1]);
  *mcptr =
      _mm_xor_si128(_mm256_extractf128_si256(cval[0], 0), _mm256_extractf128_si256(cval[0], 1));
}

ATTR_TARGET("avx2")
void mzd_mul_v_avx_128(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  word const* vptr     = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  __m128i* mcptr       = ASSUME_ALIGNED(FIRST_ROW(c), alignof(__m128i));
  __m256i const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(__m256i));

  __m256i cval[2] ATTR_ALIGNED(alignof(__m256i)) = {_mm256_setzero_si256(), _mm256_setzero_si256()};
  for (unsigned int w = 2; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int i = sizeof(word) * 8; i; i -= 8, idx >>= 8, mAptr += 4) {
      cval[0] = mm256_xor_mask(cval[0], mAptr[0], mm256_compute_mask_2(idx, 0));
      cval[1] = mm256_xor_mask(cval[1], mAptr[1], mm256_compute_mask_2(idx, 2));
      cval[0] = mm256_xor_mask(cval[0], mAptr[2], mm256_compute_mask_2(idx, 4));
      cval[1] = mm256_xor_mask(cval[1], mAptr[3], mm256_compute_mask_2(idx, 6));
    }
  }
  cval[0] = _mm256_xor_si256(cval[0], cval[1]);
  *mcptr =
      _mm_xor_si128(_mm256_extractf128_si256(cval[0], 0), _mm256_extractf128_si256(cval[0], 1));
}

ATTR_TARGET("avx2")
void mzd_addmul_v_avx_192(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  word const* vptr     = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  __m256i* mcptr       = ASSUME_ALIGNED(FIRST_ROW(c), alignof(__m256i));
  __m256i const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(__m256i));

  __m256i cval[2] ATTR_ALIGNED(alignof(__m256i)) = {*mcptr, _mm256_setzero_si256()};
  for (unsigned int w = 3; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int i = sizeof(word) * 8; i; i -= 4, idx >>= 4, mAptr += 4) {
      cval[0] = mm256_xor_mask(cval[0], mAptr[0], mm256_compute_mask(idx, 0));
      cval[1] = mm256_xor_mask(cval[1], mAptr[1], mm256_compute_mask(idx, 1));
      cval[0] = mm256_xor_mask(cval[0], mAptr[2], mm256_compute_mask(idx, 2));
      cval[1] = mm256_xor_mask(cval[1], mAptr[3], mm256_compute_mask(idx, 3));
    }
  }
  *mcptr = _mm256_xor_si256(cval[0], cval[1]);
}

ATTR_TARGET("avx2")
void mzd_mul_v_avx_192(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  word const* vptr     = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  __m256i* mcptr       = ASSUME_ALIGNED(FIRST_ROW(c), alignof(__m256i));
  __m256i const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(__m256i));

  __m256i cval[2] ATTR_ALIGNED(alignof(__m256i)) = {_mm256_setzero_si256(), _mm256_setzero_si256()};
  for (unsigned int w = 3; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int i = sizeof(word) * 8; i; i -= 4, idx >>= 4, mAptr += 4) {
      cval[0] = mm256_xor_mask(cval[0], mAptr[0], mm256_compute_mask(idx, 0));
      cval[1] = mm256_xor_mask(cval[1], mAptr[1], mm256_compute_mask(idx, 1));
      cval[0] = mm256_xor_mask(cval[0], mAptr[2], mm256_compute_mask(idx, 2));
      cval[1] = mm256_xor_mask(cval[1], mAptr[3], mm256_compute_mask(idx, 3));
    }
  }
  *mcptr = _mm256_xor_si256(cval[0], cval[1]);
}

ATTR_TARGET("avx2")
void mzd_addmul_v_avx_256(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  word const* vptr     = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  __m256i* mcptr       = ASSUME_ALIGNED(FIRST_ROW(c), alignof(__m256i));
  __m256i const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(__m256i));

  __m256i cval[2] ATTR_ALIGNED(alignof(__m256i)) = {*mcptr, _mm256_setzero_si256()};
  for (unsigned int w = 4; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int i = sizeof(word) * 8; i; i -= 4, idx >>= 4, mAptr += 4) {
      cval[0] = mm256_xor_mask(cval[0], mAptr[0], mm256_compute_mask(idx, 0));
      cval[1] = mm256_xor_mask(cval[1], mAptr[1], mm256_compute_mask(idx, 1));
      cval[0] = mm256_xor_mask(cval[0], mAptr[2], mm256_compute_mask(idx, 2));
      cval[1] = mm256_xor_mask(cval[1], mAptr[3], mm256_compute_mask(idx, 3));
    }
  }
  *mcptr = _mm256_xor_si256(cval[0], cval[1]);
}

ATTR_TARGET("avx2")
void mzd_mul_v_avx_256(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  word const* vptr     = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  __m256i* mcptr       = ASSUME_ALIGNED(FIRST_ROW(c), alignof(__m256i));
  __m256i const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(__m256i));

  __m256i cval[2] ATTR_ALIGNED(alignof(__m256i)) = {_mm256_setzero_si256(), _mm256_setzero_si256()};
  for (unsigned int w = 4; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int i = sizeof(word) * 8; i; i -= 4, idx >>= 4, mAptr += 4) {
      cval[0] = mm256_xor_mask(cval[0], mAptr[0], mm256_compute_mask(idx, 0));
      cval[1] = mm256_xor_mask(cval[1], mAptr[1], mm256_compute_mask(idx, 1));
      cval[0] = mm256_xor_mask(cval[0], mAptr[2], mm256_compute_mask(idx, 2));
      cval[1] = mm256_xor_mask(cval[1], mAptr[3], mm256_compute_mask(idx, 3));
    }
  }
  *mcptr = _mm256_xor_si256(cval[0], cval[1]);
}

#endif

#if defined(WITH_NEON)
ATTR_CONST
static inline uint32x4_t mm128_compute_mask(const word idx, size_t bit) {
  return vreinterpretq_u32_u64(vdupq_n_u64(-((idx >> bit) & 1)));
}

void mzd_mul_v_neon(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  mzd_local_clear(c);
  mzd_addmul_v_neon(c, v, A);
}

void mzd_addmul_v_neon(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  const unsigned int width      = v->width;
  const unsigned int rowstride  = A->rowstride;
  const unsigned int mrowstride = rowstride * sizeof(word) / sizeof(uint32x4_t);
  const unsigned int len        = mrowstride;

  word const* vptr  = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  uint32x4_t* mcptr = ASSUME_ALIGNED(FIRST_ROW(c), alignof(uint32x4_t));

  for (unsigned int w = 0; w < width; ++w, ++vptr) {
    word idx = *vptr;
    uint32x4_t const* mAptr =
        ASSUME_ALIGNED(CONST_ROW(A, w * sizeof(word) * 8), alignof(uint32x4_t));

    for (unsigned int i = sizeof(word) * 8; i; --i, idx >>= 1, mAptr += mrowstride) {
      mm128_xor_mask_region(mcptr, mAptr, mm128_compute_mask(idx, 0), len);
    }
  }
}

void mzd_mul_v_neon_128(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  word const* vptr = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);

  uint32x4_t* mcptr       = ASSUME_ALIGNED(FIRST_ROW(c), alignof(uint32x4_t));
  uint32x4_t const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(uint32x4_t));

  uint32x4_t cval[2] ATTR_ALIGNED(alignof(uint32x4_t)) = {vmovq_n_u32(0), vmovq_n_u32(0)};
  for (unsigned int w = 2; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int i = sizeof(word) * 8; i; i -= 4, idx >>= 4, mAptr += 4) {
      cval[0] = mm128_xor_mask(cval[0], mAptr[0], mm128_compute_mask(idx, 0));
      cval[1] = mm128_xor_mask(cval[1], mAptr[1], mm128_compute_mask(idx, 1));
      cval[0] = mm128_xor_mask(cval[0], mAptr[2], mm128_compute_mask(idx, 2));
      cval[1] = mm128_xor_mask(cval[1], mAptr[3], mm128_compute_mask(idx, 3));
    }
  }
  *mcptr = veorq_u32(cval[0], cval[1]);
}

void mzd_addmul_v_neon_128(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  word const* vptr = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);

  uint32x4_t* mcptr       = ASSUME_ALIGNED(FIRST_ROW(c), alignof(uint32x4_t));
  uint32x4_t const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(uint32x4_t));

  uint32x4_t cval[2] ATTR_ALIGNED(alignof(uint32x4_t)) = {*mcptr, vmovq_n_u32(0)};
  for (unsigned int w = 2; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int i = sizeof(word) * 8; i; i -= 4, idx >>= 4, mAptr += 4) {
      cval[0] = mm128_xor_mask(cval[0], mAptr[0], mm128_compute_mask(idx, 0));
      cval[1] = mm128_xor_mask(cval[1], mAptr[1], mm128_compute_mask(idx, 1));
      cval[0] = mm128_xor_mask(cval[0], mAptr[2], mm128_compute_mask(idx, 2));
      cval[1] = mm128_xor_mask(cval[1], mAptr[3], mm128_compute_mask(idx, 3));
    }
  }
  *mcptr = veorq_u32(cval[0], cval[1]);
}

void mzd_mul_v_neon_192(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  word const* vptr = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);

  uint32x4_t* mcptr       = ASSUME_ALIGNED(FIRST_ROW(c), alignof(uint32x4_t));
  uint32x4_t const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(uint32x4_t));

  uint32x4_t cval[4] ATTR_ALIGNED(alignof(uint32x4_t)) = {vmovq_n_u32(0), vmovq_n_u32(0),
                                                          vmovq_n_u32(0), vmovq_n_u32(0)};
  for (unsigned int w = 3; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int i = sizeof(word) * 8; i; i -= 2, idx >>= 2, mAptr += 4) {
      mm128_xor_mask_region(&cval[0], mAptr + 0, mm128_compute_mask(idx, 0), 2);
      mm128_xor_mask_region(&cval[2], mAptr + 2, mm128_compute_mask(idx, 1), 2);
    }
  }
  mcptr[0] = veorq_u32(cval[0], cval[2]);
  mcptr[1] = veorq_u32(cval[1], cval[3]);
}

void mzd_addmul_v_neon_192(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  word const* vptr = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);

  uint32x4_t* mcptr       = ASSUME_ALIGNED(FIRST_ROW(c), alignof(uint32x4_t));
  uint32x4_t const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(uint32x4_t));

  uint32x4_t cval[4] ATTR_ALIGNED(alignof(uint32x4_t)) = {mcptr[0], mcptr[1], vmovq_n_u32(0),
                                                          vmovq_n_u32(0)};
  for (unsigned int w = 3; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int i = sizeof(word) * 8; i; i -= 2, idx >>= 2, mAptr += 4) {
      mm128_xor_mask_region(&cval[0], mAptr + 0, mm128_compute_mask(idx, 0), 2);
      mm128_xor_mask_region(&cval[2], mAptr + 2, mm128_compute_mask(idx, 1), 2);
    }
  }
  mcptr[0] = veorq_u32(cval[0], cval[2]);
  mcptr[1] = veorq_u32(cval[1], cval[3]);
}

void mzd_mul_v_neon_256(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  word const* vptr = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);

  uint32x4_t* mcptr       = ASSUME_ALIGNED(FIRST_ROW(c), alignof(uint32x4_t));
  uint32x4_t const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(uint32x4_t));

  uint32x4_t cval[4] ATTR_ALIGNED(alignof(uint32x4_t)) = {vmovq_n_u32(0), vmovq_n_u32(0),
                                                          vmovq_n_u32(0), vmovq_n_u32(0)};
  for (unsigned int w = 4; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int i = sizeof(word) * 8; i; i -= 2, idx >>= 2, mAptr += 4) {
      mm128_xor_mask_region(&cval[0], mAptr + 0, mm128_compute_mask(idx, 0), 2);
      mm128_xor_mask_region(&cval[2], mAptr + 2, mm128_compute_mask(idx, 1), 2);
    }
  }
  mcptr[0] = veorq_u32(cval[0], cval[2]);
  mcptr[1] = veorq_u32(cval[1], cval[3]);
}

void mzd_addmul_v_neon_256(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  word const* vptr = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);

  uint32x4_t* mcptr       = ASSUME_ALIGNED(FIRST_ROW(c), alignof(uint32x4_t));
  uint32x4_t const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(uint32x4_t));

  uint32x4_t ATTR_ALIGNED(alignof(uint32x4_t))
      cval[4] = {mcptr[0], mcptr[1], vmovq_n_u32(0), vmovq_n_u32(0)};
  for (unsigned int w = 4; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int i = sizeof(word) * 8; i; i -= 2, idx >>= 2, mAptr += 4) {
      mm128_xor_mask_region(&cval[0], mAptr + 0, mm128_compute_mask(idx, 0), 2);
      mm128_xor_mask_region(&cval[2], mAptr + 2, mm128_compute_mask(idx, 1), 2);
    }
  }
  mcptr[0] = veorq_u32(cval[0], cval[2]);
  mcptr[1] = veorq_u32(cval[1], cval[3]);
}
#endif
#endif



void mzd_addmul_v(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
#if defined(WITH_OPT)
  if (A->nrows % (sizeof(word) * 8) == 0) {
#if defined(WITH_AVX2)
    if (CPU_SUPPORTS_AVX2 && (A->ncols & 0xff) == 0) {
      mzd_addmul_v_avx(c, v, A);
      return;
    }
#endif
#if defined(WITH_SSE2)
    if (CPU_SUPPORTS_SSE2 && (A->ncols & 0x7f) == 0) {
      mzd_addmul_v_sse(c, v, A);
      return;
    }
#endif
#if defined(WITH_NEON)
    if (CPU_SUPPORTS_NEON && (A->ncols & 0x7f) == 0) {
      mzd_addmul_v_neon(c, v, A);
      return;
    }
#endif
  }
#endif

  mzd_addmul_v_uint64(c, v, A);
}

void mzd_addmul_v_uint64(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  const unsigned int len       = A->width;
  const unsigned int rowstride = A->rowstride;
  word* cptr                   = ASSUME_ALIGNED(FIRST_ROW(c), 32);
  word const* vptr             = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  const unsigned int width     = v->width;
  word const* Aptr             = ASSUME_ALIGNED(CONST_FIRST_ROW(A), 32);

  for (unsigned int w = width; w; --w, ++vptr) {
    word idx = *vptr;

    for (unsigned int i = sizeof(word) * 8; i; --i, idx >>= 1, Aptr += rowstride) {
      const uint64_t mask = -(idx & 1);
      for (unsigned int j = 0; j < len; ++j) {
        cptr[j] ^= (Aptr[j] & mask);
      }
    }
  }
}

bool mzd_local_equal(mzd_local_t const* first, mzd_local_t const* second) {
  if (first == second) {
    return true;
  }
  if (first->ncols != second->ncols || first->nrows != second->nrows) {
    return false;
  }

  const unsigned int rows  = first->nrows;
  const unsigned int width = first->width;

  for (unsigned int r = 0; r < rows; ++r) {
    if (memcmp(ASSUME_ALIGNED(CONST_ROW(first, r), 32), ASSUME_ALIGNED(CONST_ROW(second, r), 32),
               sizeof(word) * width) != 0) {
      return false;
    }
  }

  return true;
}

static void xor_comb(const unsigned int len, word* Brow, mzd_local_t const* A,
                     unsigned int r_offset, unsigned comb) {
  while (comb) {
    const word* Arow = CONST_ROW(A, r_offset);
    if (comb & 0x1) {
      for (unsigned int i = 0; i < len; ++i) {
        Brow[i] ^= Arow[i];
      }
    }

    comb >>= 1;
    ++r_offset;
  }
}

/**
 * Pre-compute matrices for faster mzd_addmul_v computions.
 */
mzd_local_t* mzd_precompute_matrix_lookup(mzd_local_t const* A) {
  mzd_local_t* B = mzd_local_init_ex(32 * A->nrows, A->ncols, true);

  const unsigned int len = A->width;
  for (unsigned int r = 0; r < B->nrows; ++r) {
    const unsigned int comb     = r & 0xff;
    const unsigned int r_offset = (r >> 8) << 3;
    if (!comb) {
      continue;
    }

    xor_comb(len, ROW(B, r), A, r_offset, comb);
  }

  return B;
}

#if defined(WITH_OPT)
#if defined(WITH_SSE2)
ATTR_TARGET("sse2")
void mzd_mul_vl_sse(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  mzd_local_clear(c);
  mzd_addmul_vl_sse(c, v, A);
}

ATTR_TARGET("sse2")
void mzd_addmul_vl_sse(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  word const* vptr              = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  const unsigned int width      = v->width;
  const unsigned int rowstride  = A->rowstride;
  const unsigned int mrowstride = rowstride * sizeof(word) / sizeof(__m128i);
  const unsigned int len        = mrowstride;
  const unsigned int moff2      = 256 * mrowstride;

  __m128i* mcptr       = ASSUME_ALIGNED(FIRST_ROW(c), alignof(__m128i));
  __m128i const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(__m128i));

  for (unsigned int w = width; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int s = sizeof(word); s; --s, idx >>= 8, mAptr += moff2) {
      const word comb = idx & 0xff;
      mm128_xor_region(mcptr, mAptr + comb * mrowstride, len);
    }
  }
}

ATTR_TARGET("sse2")
void mzd_addmul_vl_sse_128(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  static const unsigned int moff2 = 256;

  word const* vptr     = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  __m128i* mcptr       = ASSUME_ALIGNED(FIRST_ROW(c), alignof(__m128i));
  __m128i const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(__m128i));

  __m128i cval[2] ATTR_ALIGNED(alignof(__m128i)) = {*mcptr, _mm_setzero_si128()};
  for (unsigned int w = 2; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int s = sizeof(word); s; s -= 2, idx >>= 16) {
      mm128_xor_region(&cval[0], mAptr + ((idx >> 0) & 0xff), 1);
      mAptr += moff2;
      mm128_xor_region(&cval[1], mAptr + ((idx >> 8) & 0xff), 1);
      mAptr += moff2;
    }
  }
  *mcptr = _mm_xor_si128(cval[0], cval[1]);
}

ATTR_TARGET("sse2")
void mzd_mul_vl_sse_128(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  static const unsigned int moff2 = 256;

  word const* vptr     = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  __m128i* mcptr       = ASSUME_ALIGNED(FIRST_ROW(c), alignof(__m128i));
  __m128i const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(__m128i));

  __m128i cval[2] ATTR_ALIGNED(alignof(__m128i)) = {_mm_setzero_si128(), _mm_setzero_si128()};
  for (unsigned int w = 2; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int s = sizeof(word); s; s -= 2, idx >>= 16) {
      mm128_xor_region(&cval[0], mAptr + ((idx >> 0) & 0xff), 1);
      mAptr += moff2;
      mm128_xor_region(&cval[1], mAptr + ((idx >> 8) & 0xff), 1);
      mAptr += moff2;
    }
  }
  *mcptr = _mm_xor_si128(cval[0], cval[1]);
}

ATTR_TARGET("sse2")
void mzd_addmul_vl_sse_192(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  static const unsigned int moff2 = 512;

  word const* vptr     = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  __m128i* mcptr       = ASSUME_ALIGNED(FIRST_ROW(c), alignof(__m128i));
  __m128i const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(__m128i));

  __m128i cval[4] ATTR_ALIGNED(alignof(__m128i)) = {mcptr[0], mcptr[1], _mm_setzero_si128(),
                                                    _mm_setzero_si128()};
  for (unsigned int w = 3; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int s = sizeof(word); s; s -= 2, idx >>= 16) {
      mm128_xor_region(&cval[0], mAptr + 2 * ((idx >> 0) & 0xff), 2);
      mAptr += moff2;
      mm128_xor_region(&cval[2], mAptr + 2 * ((idx >> 8) & 0xff), 2);
      mAptr += moff2;
    }
  }
  mcptr[0] = _mm_xor_si128(cval[0], cval[2]);
  mcptr[1] = _mm_xor_si128(cval[1], cval[3]);
}

ATTR_TARGET("sse2")
void mzd_mul_vl_sse_192(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  static const unsigned int moff2 = 512;

  word const* vptr     = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  __m128i* mcptr       = ASSUME_ALIGNED(FIRST_ROW(c), alignof(__m128i));
  __m128i const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(__m128i));

  __m128i cval[4] ATTR_ALIGNED(alignof(__m128i)) = {_mm_setzero_si128(), _mm_setzero_si128(),
                                                    _mm_setzero_si128(), _mm_setzero_si128()};
  for (unsigned int w = 3; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int s = sizeof(word); s; s -= 2, idx >>= 16) {
      mm128_xor_region(&cval[0], mAptr + 2 * ((idx >> 0) & 0xff), 2);
      mAptr += moff2;
      mm128_xor_region(&cval[2], mAptr + 2 * ((idx >> 8) & 0xff), 2);
      mAptr += moff2;
    }
  }
  mcptr[0] = _mm_xor_si128(cval[0], cval[2]);
  mcptr[1] = _mm_xor_si128(cval[1], cval[3]);
}

ATTR_TARGET("sse2")
void mzd_addmul_vl_sse_256(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  static const unsigned int moff2 = 512;

  word const* vptr     = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  __m128i* mcptr       = ASSUME_ALIGNED(FIRST_ROW(c), alignof(__m128i));
  __m128i const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(__m128i));

  __m128i cval[4] ATTR_ALIGNED(alignof(__m128i)) = {mcptr[0], mcptr[1], _mm_setzero_si128(),
                                                    _mm_setzero_si128()};
  for (unsigned int w = 4; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int s = sizeof(word); s; s -= 2, idx >>= 16) {
      mm128_xor_region(&cval[0], mAptr + 2 * ((idx >> 0) & 0xff), 2);
      mAptr += moff2;
      mm128_xor_region(&cval[2], mAptr + 2 * ((idx >> 8) & 0xff), 2);
      mAptr += moff2;
    }
  }
  mcptr[0] = _mm_xor_si128(cval[0], cval[2]);
  mcptr[1] = _mm_xor_si128(cval[1], cval[3]);
}

ATTR_TARGET("sse2")
void mzd_mul_vl_sse_256(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  static const unsigned int moff2 = 512;

  word const* vptr     = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  __m128i* mcptr       = ASSUME_ALIGNED(FIRST_ROW(c), alignof(__m128i));
  __m128i const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(__m128i));

  __m128i cval[4] ATTR_ALIGNED(alignof(__m128i)) = {_mm_setzero_si128(), _mm_setzero_si128(),
                                                    _mm_setzero_si128(), _mm_setzero_si128()};
  for (unsigned int w = 4; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int s = sizeof(word); s; s -= 2, idx >>= 16) {
      mm128_xor_region(&cval[0], mAptr + 2 * ((idx >> 0) & 0xff), 2);
      mAptr += moff2;
      mm128_xor_region(&cval[2], mAptr + 2 * ((idx >> 8) & 0xff), 2);
      mAptr += moff2;
    }
  }
  mcptr[0] = _mm_xor_si128(cval[0], cval[2]);
  mcptr[1] = _mm_xor_si128(cval[1], cval[3]);
}
#endif

#if defined(WITH_AVX2)
ATTR_TARGET("avx2")
void mzd_mul_vl_avx_256(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  word const* vptr                = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  static const unsigned int moff2 = 256;

  __m256i* mcptr       = ASSUME_ALIGNED(FIRST_ROW(c), alignof(__m256i));
  __m256i const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(__m256i));

  __m256i cval[2] ATTR_ALIGNED(alignof(__m256i)) = {_mm256_setzero_si256(), _mm256_setzero_si256()};
  for (unsigned int w = 4; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int s = sizeof(word); s; s -= 2, idx >>= 16) {
      mm256_xor_region(&cval[0], mAptr + ((idx >> 0) & 0xff), 1);
      mAptr += moff2;
      mm256_xor_region(&cval[1], mAptr + ((idx >> 8) & 0xff), 1);
      mAptr += moff2;
    }
  }
  *mcptr = _mm256_xor_si256(cval[0], cval[1]);
}

ATTR_TARGET("avx2")
void mzd_addmul_vl_avx_256(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  word const* vptr                = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  static const unsigned int moff2 = 256;

  __m256i* mcptr       = ASSUME_ALIGNED(FIRST_ROW(c), alignof(__m256i));
  __m256i const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(__m256i));

  __m256i cval[2] ATTR_ALIGNED(alignof(__m256i)) = {*mcptr, _mm256_setzero_si256()};
  for (unsigned int w = 4; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int s = sizeof(word); s; s -= 2, idx >>= 16) {
      mm256_xor_region(&cval[0], mAptr + ((idx >> 0) & 0xff), 1);
      mAptr += moff2;
      mm256_xor_region(&cval[1], mAptr + ((idx >> 8) & 0xff), 1);
      mAptr += moff2;
    }
  }
  *mcptr = _mm256_xor_si256(cval[0], cval[1]);
}

ATTR_TARGET("avx2")
void mzd_mul_vl_avx_192(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  word const* vptr                = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  static const unsigned int moff2 = 256;

  __m256i* mcptr       = ASSUME_ALIGNED(FIRST_ROW(c), alignof(__m256i));
  __m256i const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(__m256i));

  __m256i cval[2] ATTR_ALIGNED(alignof(__m256i)) = {_mm256_setzero_si256(), _mm256_setzero_si256()};
  for (unsigned int w = 3; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int s = sizeof(word); s; s -= 2, idx >>= 16) {
      mm256_xor_region(&cval[0], mAptr + ((idx >> 0) & 0xff), 1);
      mAptr += moff2;
      mm256_xor_region(&cval[1], mAptr + ((idx >> 8) & 0xff), 1);
      mAptr += moff2;
    }
  }
  *mcptr = _mm256_xor_si256(cval[0], cval[1]);
}

ATTR_TARGET("avx2")
void mzd_addmul_vl_avx_192(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  word const* vptr                = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  static const unsigned int moff2 = 256;

  __m256i* mcptr       = ASSUME_ALIGNED(FIRST_ROW(c), alignof(__m256i));
  __m256i const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(__m256i));

  __m256i cval[2] ATTR_ALIGNED(alignof(__m256i)) = {*mcptr, _mm256_setzero_si256()};
  for (unsigned int w = 3; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int s = sizeof(word); s; s -= 2, idx >>= 16) {
      mm256_xor_region(&cval[0], mAptr + ((idx >> 0) & 0xff), 1);
      mAptr += moff2;
      mm256_xor_region(&cval[1], mAptr + ((idx >> 8) & 0xff), 1);
      mAptr += moff2;
    }
  }
  *mcptr = _mm256_xor_si256(cval[0], cval[1]);
}

ATTR_TARGET("avx2")
void mzd_mul_vl_avx_128(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  word const* vptr                = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  static const unsigned int moff2 = 256;

  __m128i* mcptr       = ASSUME_ALIGNED(FIRST_ROW(c), alignof(__m128i));
  __m128i const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(__m128i));

  __m256i cval[2] ATTR_ALIGNED(alignof(__m256i)) = {_mm256_setzero_si256(), _mm256_setzero_si256()};
  for (unsigned int w = 2; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int s = sizeof(word); s; s -= 4, idx >>= 32) {
      const __m256i t1 =
          _mm256_set_m128i(mAptr[(idx >> 0) & 0xff], mAptr[((idx >> 8) & 0xff) + moff2]);
      mm256_xor_region(&cval[0], &t1, 1);
      mAptr += 2 * moff2;

      const __m256i t2 =
          _mm256_set_m128i(mAptr[(idx >> 16) & 0xff], mAptr[((idx >> 24) & 0xff) + moff2]);
      mm256_xor_region(&cval[1], &t2, 1);
      mAptr += 2 * moff2;
    }
  }
  cval[0] = _mm256_xor_si256(cval[0], cval[1]);
  *mcptr =
      _mm_xor_si128(_mm256_extractf128_si256(cval[0], 0), _mm256_extractf128_si256(cval[0], 1));
}

ATTR_TARGET("avx2")
void mzd_addmul_vl_avx_128(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  word const* vptr                = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  static const unsigned int moff2 = 256;

  __m128i* mcptr       = ASSUME_ALIGNED(FIRST_ROW(c), alignof(__m128i));
  __m128i const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(__m128i));

  __m256i cval[2] ATTR_ALIGNED(alignof(__m256i)) = {_mm256_castsi128_si256(*mcptr),
                                                    _mm256_setzero_si256()};
  for (unsigned int w = 2; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int s = sizeof(word); s; s -= 4, idx >>= 32) {
      const __m256i t1 =
          _mm256_set_m128i(mAptr[(idx >> 0) & 0xff], mAptr[((idx >> 8) & 0xff) + moff2]);
      mm256_xor_region(&cval[0], &t1, 1);
      mAptr += 2 * moff2;

      const __m256i t2 =
          _mm256_set_m128i(mAptr[(idx >> 16) & 0xff], mAptr[((idx >> 24) & 0xff) + moff2]);
      mm256_xor_region(&cval[1], &t2, 1);
      mAptr += 2 * moff2;
    }
  }
  cval[0] = _mm256_xor_si256(cval[0], cval[1]);
  *mcptr =
      _mm_xor_si128(_mm256_extractf128_si256(cval[0], 0), _mm256_extractf128_si256(cval[0], 1));
}

ATTR_TARGET("avx2")
void mzd_mul_vl_avx(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  mzd_local_clear(c);
  mzd_addmul_vl_avx(c, v, A);
}

ATTR_TARGET("avx2")
void mzd_addmul_vl_avx(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  const unsigned int width      = v->width;
  const unsigned int rowstride  = A->rowstride;
  const unsigned int mrowstride = rowstride * sizeof(word) / sizeof(__m256i);
  const unsigned int moff2      = 256 * mrowstride;
  const unsigned int len        = mrowstride;

  word const* vptr     = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  __m256i* mcptr       = ASSUME_ALIGNED(FIRST_ROW(c), alignof(__m256i));
  __m256i const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(__m256i));

  for (unsigned int w = width; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int s = sizeof(word); s; --s, idx >>= 8, mAptr += moff2) {
      const word comb = idx & 0xff;
      mm256_xor_region(mcptr, mAptr + comb * mrowstride, len);
    }
  }
}
#endif

#if defined(WITH_NEON)
void mzd_mul_vl_neon_128(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  static const unsigned int moff2 = 256;

  word const* vptr        = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  uint32x4_t* mcptr       = ASSUME_ALIGNED(FIRST_ROW(c), alignof(uint32x4_t));
  uint32x4_t const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(uint32x4_t));

  uint32x4_t mc ATTR_ALIGNED(alignof(uint32x4_t)) = vmovq_n_u32(0);
  for (unsigned int w = 2; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int s = sizeof(word); s; --s, idx >>= 8, mAptr += moff2) {
      const word comb = idx & 0xff;
      mc              = veorq_u32(mc, mAptr[comb]);
    }
  }

  *mcptr = mc;
}

void mzd_addmul_vl_neon_128(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  static const unsigned int moff2 = 256;

  word const* vptr        = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  uint32x4_t* mcptr       = ASSUME_ALIGNED(FIRST_ROW(c), alignof(uint32x4_t));
  uint32x4_t const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(uint32x4_t));

  uint32x4_t mc ATTR_ALIGNED(alignof(uint32x4_t)) = *mcptr;
  for (unsigned int w = 2; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int s = sizeof(word); s; --s, idx >>= 8, mAptr += moff2) {
      const word comb = idx & 0xff;
      mc              = veorq_u32(mc, mAptr[comb]);
    }
  }
  *mcptr = mc;
}

void mzd_mul_vl_neon(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  mzd_local_clear(c);
  mzd_addmul_vl_neon(c, v, A);
}

void mzd_addmul_vl_neon(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  word const* vptr              = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  const unsigned int width      = v->width;
  const unsigned int rowstride  = A->rowstride;
  const unsigned int mrowstride = rowstride * sizeof(word) / sizeof(uint32x4_t);
  const unsigned int len        = mrowstride;
  const unsigned int moff2      = 256 * mrowstride;

  uint32x4_t* mcptr       = ASSUME_ALIGNED(FIRST_ROW(c), alignof(uint32x4_t));
  uint32x4_t const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(uint32x4_t));

  for (unsigned int w = width; w; --w, ++vptr) {
    word idx = *vptr;
    for (unsigned int s = sizeof(word); s; --s, idx >>= 8, mAptr += moff2) {
      const word comb = idx & 0xff;
      mm128_xor_region(mcptr, mAptr + comb * mrowstride, len);
    }
  }
}
#endif
#endif

void mzd_mul_vl(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
#if defined(WITH_OPT)
  if (A->nrows % (sizeof(word) * 8) == 0) {
#if defined(WITH_AVX2)
    if (CPU_SUPPORTS_AVX2) {
      if (A->ncols == 256 && v->ncols == 256) {
        mzd_mul_vl_avx_256(c, v, A);
        return;
      }
    }
#endif
#if defined(WITH_SSE2)
    if (CPU_SUPPORTS_SSE2) {
      if (A->ncols == 128 && v->ncols == 128) {
        mzd_mul_vl_sse_128(c, v, A);
        return;
      }
    }
#endif
#if defined(WITH_NEON)
    if (CPU_SUPPORTS_NEON) {
      if (A->ncols == 128 && v->ncols == 128) {
        mzd_mul_vl_neon_128(c, v, A);
        return;
      }
    }
#endif
  }
#endif
  mzd_local_clear(c);
  mzd_addmul_vl(c, v, A);
}

void mzd_mul_vl_uint64(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  mzd_local_clear(c);
  mzd_addmul_vl_uint64(c, v, A);
}

void mzd_addmul_vl(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
#if defined(WITH_OPT)
  if (A->nrows % (sizeof(word) * 8) == 0) {
#if defined(WITH_AVX2)
    if (CPU_SUPPORTS_AVX2) {
      if (A->ncols == 256 && v->ncols == 256) {
        mzd_addmul_vl_avx_256(c, v, A);
        return;
      }
      if ((A->ncols & 0xff) == 0) {
        mzd_addmul_vl_avx(c, v, A);
        return;
      }
    }
#endif
#if defined(WITH_SSE2)
    if (CPU_SUPPORTS_SSE2) {
      if (A->ncols == 128 && v->ncols == 128) {
        mzd_addmul_vl_sse_128(c, v, A);
        return;
      }
      if ((A->ncols & 0x7f) == 0) {
        mzd_addmul_vl_sse(c, v, A);
        return;
      }
    }
#endif
#if defined(WITH_NEON)
    if (CPU_SUPPORTS_NEON) {
      if (A->ncols == 128 && v->ncols == 128) {
        mzd_addmul_vl_neon_128(c, v, A);
        return;
      }
      if ((A->ncols & 0x7f) == 0) {
        mzd_addmul_vl_neon(c, v, A);
        return;
      }
    }
#endif
  }
#endif
  mzd_addmul_vl_uint64(c, v, A);
  return;
}

void mzd_addmul_vl_uint64(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  const unsigned int len   = A->width;
  word* cptr               = ASSUME_ALIGNED(FIRST_ROW(c), 32);
  word const* vptr         = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  const unsigned int width = v->width;

  for (unsigned int w = 0; w < width; ++w, ++vptr) {
    word idx         = *vptr;
    unsigned int add = 0;

    while (idx) {
      const word comb = idx & 0xff;

      word const* Aptr = CONST_ROW(A, w * sizeof(word) * 8 * 32 + add + comb);
      for (unsigned int i = 0; i < len; ++i) {
        cptr[i] ^= Aptr[i];
      }

      idx >>= 8;
      add += 256;
    }
  }
}

// specific instances
#if defined(REDUCED_LINEAR_LAYER_NEXT)
//no simd
void mzd_mul_v_uint64_30(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
    const unsigned int rowstride = A->rowstride;
    word* cptr                   = ASSUME_ALIGNED(FIRST_ROW(c), 32);
    word const* vptr             = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
    word const* Aptr             = ASSUME_ALIGNED(CONST_FIRST_ROW(A), 32);
    const unsigned int width     = v->width;

    word idx = vptr[width-1] >> 34;
    for (unsigned int i = 30; i; --i, idx >>= 1, Aptr += rowstride) {
      const uint64_t mask = -(idx & 1);
      for(unsigned int j = 0; j < width; j++) {
        cptr[j] ^= (Aptr[j] & mask);
      }
    }
}

void mzd_mul_v_uint64_3(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  const unsigned int rowstride = A->rowstride;
  word* cptr                   = ASSUME_ALIGNED(FIRST_ROW(c), 32);
  word const* vptr             = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  word const* Aptr             = ASSUME_ALIGNED(CONST_FIRST_ROW(A), 32);
  const unsigned int width     = v->width;

  word idx = vptr[width-1] >> 61;
  for (unsigned int i = 3; i; --i, idx >>= 1, Aptr += rowstride) {
    const uint64_t mask = -(idx & 1);
    for(unsigned int j = 0; j < width; j++) {
      cptr[j] ^= (Aptr[j] & mask);
    }
  }
}

//popcnt
void mzd_mul_v_30_popcnt(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* At) {
  const unsigned int rowstride = At->rowstride;
  word* cptr                   = ASSUME_ALIGNED(FIRST_ROW(c), 32);
  word const* vptr             = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  word const* Aptr             = ASSUME_ALIGNED(CONST_FIRST_ROW(At), 32);
  const unsigned int width     = v->width;

  for(unsigned int j = 0; j < width; j++) {
    cptr[j] = 0;
  }
  for(unsigned i = 30; i; --i) {
    word const* A = Aptr + (30-i)*rowstride;
    word popcnt = 0;
    for(unsigned int j = 0; j < width; j++) {
      popcnt += __builtin_popcountll(vptr[j] & A[j]);
    }
    cptr[width-1] |= (popcnt & WORD_C(0x1)) << (64-i);
  }
}

void mzd_mul_v_98_30_popcnt(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* At) {
  const unsigned int rowstride = At->rowstride;
  word* cptr                   = ASSUME_ALIGNED(FIRST_ROW(c), 32);
  word const* vptr             = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  word const* Aptr             = ASSUME_ALIGNED(CONST_FIRST_ROW(At), 32);

  for(unsigned int j = 0; j < 2; j++) {
    cptr[j] = 0;
  }
  for(unsigned i = 30; i; --i) {
    word const* A = Aptr + (30-i)*rowstride;
    word popcnt = __builtin_popcountll(vptr[0] & A[0])
                + __builtin_popcountll(vptr[1] & A[1]);
    cptr[1] |= (popcnt & WORD_C(0x1)) << (64-i);
  }
}
void mzd_mul_v_162_30_popcnt(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* At) {
  const unsigned int rowstride = At->rowstride;
  word* cptr                   = ASSUME_ALIGNED(FIRST_ROW(c), 32);
  word const* vptr             = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  word const* Aptr             = ASSUME_ALIGNED(CONST_FIRST_ROW(At), 32);

  for(unsigned int j = 0; j < 3; j++) {
    cptr[j] = 0;
  }
  for(unsigned i = 30; i; --i) {
    word const* A = Aptr + (30-i)*rowstride;
    word popcnt = __builtin_popcountll(vptr[0] & A[0])
                + __builtin_popcountll(vptr[1] & A[1])
                + __builtin_popcountll(vptr[2] & A[2]);
    cptr[2] |= (popcnt & WORD_C(0x1)) << (64-i);
  }
}
void mzd_mul_v_226_30_popcnt(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* At) {
  const unsigned int rowstride = At->rowstride;
  word* cptr                   = ASSUME_ALIGNED(FIRST_ROW(c), 32);
  word const* vptr             = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  word const* Aptr             = ASSUME_ALIGNED(CONST_FIRST_ROW(At), 32);

  for(unsigned int j = 0; j < 4; j++) {
    cptr[j] = 0;
  }
  for(unsigned i = 30; i; --i) {
    word const* A = Aptr + (30-i)*rowstride;
    word popcnt = __builtin_popcountll(vptr[0] & A[0])
                + __builtin_popcountll(vptr[1] & A[1])
                + __builtin_popcountll(vptr[2] & A[2])
                + __builtin_popcountll(vptr[3] & A[3]);
    cptr[3] |= (popcnt & WORD_C(0x1)) << (64-i);
  }
}

void mzd_mul_v_3_popcnt(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* At) {
  const unsigned int rowstride = At->rowstride;
  word* cptr                   = ASSUME_ALIGNED(FIRST_ROW(c), 32);
  word const* vptr             = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  word const* Aptr             = ASSUME_ALIGNED(CONST_FIRST_ROW(At), 32);
  const unsigned int width     = v->width;

  for(unsigned int j = 0; j < width; j++) {
    cptr[j] = 0;
  }
  for(unsigned i = 3; i; --i) {
    word const* A = Aptr + (3-i)*rowstride;
    word popcnt = 0;
    for(unsigned int j = 0; j < width; j++) {
      popcnt += __builtin_popcountll(vptr[j] & A[j]);
    }
    cptr[width-1] |= (popcnt & WORD_C(0x1)) << (64-i);
  }
}

void mzd_mul_v_125_3_popcnt(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* At) {
  const unsigned int rowstride = At->rowstride;
  word* cptr                   = ASSUME_ALIGNED(FIRST_ROW(c), 32);
  word const* vptr             = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  word const* Aptr             = ASSUME_ALIGNED(CONST_FIRST_ROW(At), 32);

  for(unsigned int j = 0; j < 2; j++) {
    cptr[j] = 0;
  }
  for(unsigned i = 3; i; --i) {
    word const* A = Aptr + (3-i)*rowstride;
    word popcnt = __builtin_popcountll(vptr[0] & A[0])
                  + __builtin_popcountll(vptr[1] & A[1]);
   cptr[1] |= (popcnt & WORD_C(0x1)) << (64-i);
  }
}

void mzd_mul_v_189_3_popcnt(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* At) {
  const unsigned int rowstride = At->rowstride;
  word* cptr                   = ASSUME_ALIGNED(FIRST_ROW(c), 32);
  word const* vptr             = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  word const* Aptr             = ASSUME_ALIGNED(CONST_FIRST_ROW(At), 32);

  for(unsigned int j = 0; j < 4; j++) {
    cptr[j] = 0;
  }
  for(unsigned i = 3; i; --i) {
    word const* A = Aptr + (3-i)*rowstride;
    word popcnt = __builtin_popcountll(vptr[0] & A[0])
                  + __builtin_popcountll(vptr[1] & A[1])
                  + __builtin_popcountll(vptr[2] & A[2]);
    cptr[2] |= (popcnt & WORD_C(0x1)) << (64-i);
  }
}

void mzd_mul_v_253_3_popcnt(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* At) {
  const unsigned int rowstride = At->rowstride;
  word* cptr                   = ASSUME_ALIGNED(FIRST_ROW(c), 32);
  word const* vptr             = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  word const* Aptr             = ASSUME_ALIGNED(CONST_FIRST_ROW(At), 32);

  for(unsigned int j = 0; j < 4; j++) {
    cptr[j] = 0;
  }
  for(unsigned i = 3; i; --i) {
    word const* A = Aptr + (3-i)*rowstride;
    word popcnt = __builtin_popcountll(vptr[0] & A[0])
                  + __builtin_popcountll(vptr[1] & A[1])
                  + __builtin_popcountll(vptr[2] & A[2])
                  + __builtin_popcountll(vptr[3] & A[3]);
    cptr[3] |= (popcnt & WORD_C(0x1)) << (64-i);
  }
}

#if defined(WITH_SSE2)
ATTR_TARGET("sse2")
void mzd_mul_v_sse_30_128(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  word const* vptr     = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  __m128i* mcptr       = ASSUME_ALIGNED(FIRST_ROW(c), alignof(__m128i));
  __m128i const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(__m128i));

  __m128i cval[2] ATTR_ALIGNED(alignof(__m128i)) = {_mm_setzero_si128(), _mm_setzero_si128()};
  word idx = vptr[1] >> 34;
  for (unsigned int i = 28; i; i -= 4, idx >>= 4, mAptr += 4) {
    mm128_xor_mask_region(&cval[0], mAptr + 0, _mm_set1_epi64x(-(idx & 1)), 1);
    mm128_xor_mask_region(&cval[1], mAptr + 1, _mm_set1_epi64x(-((idx >> 1) & 1)), 1);
    mm128_xor_mask_region(&cval[0], mAptr + 2, _mm_set1_epi64x(-((idx >> 2) & 1)), 1);
    mm128_xor_mask_region(&cval[1], mAptr + 3, _mm_set1_epi64x(-((idx >> 3) & 1)), 1);
  }
  mm128_xor_mask_region(&cval[0], mAptr + 0, _mm_set1_epi64x(-(idx & 1)), 1);
  mm128_xor_mask_region(&cval[1], mAptr + 1, _mm_set1_epi64x(-((idx >> 1) & 1)), 1);
  *mcptr ^= _mm_xor_si128(cval[0], cval[1]);
}

ATTR_TARGET("sse2")
void mzd_mul_v_sse_30_192(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  word const* vptr     = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  __m128i* mcptr       = ASSUME_ALIGNED(FIRST_ROW(c), alignof(__m128i));
  __m128i const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(__m128i));

  __m128i cval[4] ATTR_ALIGNED(alignof(__m128i)) = {_mm_setzero_si128(), _mm_setzero_si128(),
                                                    _mm_setzero_si128(), _mm_setzero_si128()};
  word idx = vptr[2] >> 34;
  for (unsigned int i = 30; i; i -= 2, idx >>= 2, mAptr += 4) {
    mm128_xor_mask_region(&cval[0], mAptr + 0, _mm_set1_epi64x(-(idx & 1)), 2);
    mm128_xor_mask_region(&cval[2], mAptr + 2, _mm_set1_epi64x(-((idx >> 1) & 1)), 2);
  }
  mcptr[0] ^= _mm_xor_si128(cval[0], cval[2]);
  mcptr[1] ^= _mm_xor_si128(cval[1], cval[3]);
}

ATTR_TARGET("sse2")
void mzd_mul_v_sse_30_256(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  word const* vptr     = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  __m128i* mcptr       = ASSUME_ALIGNED(FIRST_ROW(c), alignof(__m128i));
  __m128i const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(__m128i));

  __m128i cval[4] ATTR_ALIGNED(alignof(__m128i)) = {_mm_setzero_si128(), _mm_setzero_si128(),
                                                    _mm_setzero_si128(), _mm_setzero_si128()};
  word idx = vptr[3] >> 34;
  for (unsigned int i = 30; i; i -= 2, idx >>= 2, mAptr += 4) {
    mm128_xor_mask_region(&cval[0], mAptr + 0, _mm_set1_epi64x(-(idx & 1)), 2);
    mm128_xor_mask_region(&cval[2], mAptr + 2, _mm_set1_epi64x(-((idx >> 1) & 1)), 2);
  }
  mcptr[0] ^= _mm_xor_si128(cval[0], cval[2]);
  mcptr[1] ^= _mm_xor_si128(cval[1], cval[3]);
}

ATTR_TARGET("sse2")
void mzd_mul_v_sse_3_128(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  word const* vptr     = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  __m128i* mcptr       = ASSUME_ALIGNED(FIRST_ROW(c), alignof(__m128i));
  __m128i const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(__m128i));

  __m128i cval[2] ATTR_ALIGNED(alignof(__m128i)) = {_mm_setzero_si128(), _mm_setzero_si128()};
  word idx = vptr[1] >> 61;
  mm128_xor_mask_region(&cval[0], mAptr + 0, _mm_set1_epi64x(-(idx & 1)), 1);
  mm128_xor_mask_region(&cval[1], mAptr + 1, _mm_set1_epi64x(-((idx >> 1) & 1)), 1);
  mm128_xor_mask_region(&cval[0], mAptr + 2, _mm_set1_epi64x(-((idx >> 2) & 1)), 1);
  *mcptr ^= _mm_xor_si128(cval[0], cval[1]);
}

ATTR_TARGET("sse2")
void mzd_mul_v_sse_3_192(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  word const* vptr     = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  __m128i* mcptr       = ASSUME_ALIGNED(FIRST_ROW(c), alignof(__m128i));
  __m128i const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(__m128i));

  __m128i cval[4] ATTR_ALIGNED(alignof(__m128i)) = {_mm_setzero_si128(), _mm_setzero_si128(),
                                                    _mm_setzero_si128(), _mm_setzero_si128()};
  word idx = vptr[2] >> 61;
  mm128_xor_mask_region(&cval[0], mAptr + 0, _mm_set1_epi64x(-(idx & 1)), 2);
  mm128_xor_mask_region(&cval[2], mAptr + 2, _mm_set1_epi64x(-((idx >> 1) & 1)), 2);
  mm128_xor_mask_region(&cval[0], mAptr + 4, _mm_set1_epi64x(-((idx >> 2) & 1)), 2);
  mcptr[0] ^= _mm_xor_si128(cval[0], cval[2]);
  mcptr[1] ^= _mm_xor_si128(cval[1], cval[3]);
}

ATTR_TARGET("sse2")
void mzd_mul_v_sse_3_256(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  word const* vptr     = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  __m128i* mcptr       = ASSUME_ALIGNED(FIRST_ROW(c), alignof(__m128i));
  __m128i const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(__m128i));

  __m128i cval[4] ATTR_ALIGNED(alignof(__m128i)) = {_mm_setzero_si128(), _mm_setzero_si128(),
                                                    _mm_setzero_si128(), _mm_setzero_si128()};
  word idx = vptr[3] >> 61;
  mm128_xor_mask_region(&cval[0], mAptr + 0, _mm_set1_epi64x(-(idx & 1)), 2);
  mm128_xor_mask_region(&cval[2], mAptr + 2, _mm_set1_epi64x(-((idx >> 1) & 1)), 2);
  mm128_xor_mask_region(&cval[0], mAptr + 4, _mm_set1_epi64x(-((idx >> 2) & 1)), 2);

  mcptr[0] ^= _mm_xor_si128(cval[0], cval[2]);
  mcptr[1] ^= _mm_xor_si128(cval[1], cval[3]);
}
#endif

#if defined(WITH_AVX2)
ATTR_TARGET("avx2")
void mzd_mul_v_avx_30_128(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  word const* vptr     = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  __m128i* mcptr       = ASSUME_ALIGNED(FIRST_ROW(c), alignof(__m128i));
  __m256i const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(__m256i));

  __m256i cval[2] ATTR_ALIGNED(alignof(__m256i)) = {_mm256_setzero_si256(), _mm256_setzero_si256()};
  word idx = vptr[1] >> 34;
  for (unsigned int i = 24; i; i -= 8, idx >>= 8, mAptr += 4) {
    const int64_t m1 = -(idx & 1);
    const int64_t m2 = -((idx >> 1) & 1);
    mm256_xor_mask_region(&cval[0], mAptr + 0, _mm256_set_epi64x(m2, m2, m1, m1), 1);
    const int64_t m3 = -((idx >> 2) & 1);
    const int64_t m4 = -((idx >> 3) & 1);
    mm256_xor_mask_region(&cval[1], mAptr + 1, _mm256_set_epi64x(m4, m4, m3, m3), 1);
    const int64_t m5 = -((idx >> 4) & 1);
    const int64_t m6 = -((idx >> 5) & 1);
    mm256_xor_mask_region(&cval[0], mAptr + 2, _mm256_set_epi64x(m6, m6, m5, m5), 1);
    const int64_t m7 = -((idx >> 6) & 1);
    const int64_t m8 = -((idx >> 7) & 1);
    mm256_xor_mask_region(&cval[1], mAptr + 3, _mm256_set_epi64x(m8, m8, m7, m7), 1);
  }
  const int64_t m1 = -(idx & 1);
  const int64_t m2 = -((idx >> 1) & 1);
  mm256_xor_mask_region(&cval[0], mAptr + 0, _mm256_set_epi64x(m2, m2, m1, m1), 1);
  const int64_t m3 = -((idx >> 2) & 1);
  const int64_t m4 = -((idx >> 3) & 1);
  mm256_xor_mask_region(&cval[1], mAptr + 1, _mm256_set_epi64x(m4, m4, m3, m3), 1);
  const int64_t m5 = -((idx >> 4) & 1);
  const int64_t m6 = -((idx >> 5) & 1);
  mm256_xor_mask_region(&cval[0], mAptr + 2, _mm256_set_epi64x(m6, m6, m5, m5), 1);

  cval[0] = _mm256_xor_si256(cval[0], cval[1]);
  *mcptr ^=
          _mm_xor_si128(_mm256_extractf128_si256(cval[0], 0), _mm256_extractf128_si256(cval[0], 1));
}

ATTR_TARGET("avx2")
void mzd_mul_v_avx_30_192(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  word const* vptr     = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  __m256i* mcptr       = ASSUME_ALIGNED(FIRST_ROW(c), alignof(__m256i));
  __m256i const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(__m256i));

  __m256i cval[2] ATTR_ALIGNED(alignof(__m256i)) = {_mm256_setzero_si256(), _mm256_setzero_si256()};
  word idx = vptr[2] >> 34;
  // do 7x4 and then 2 extra to get 30
  for (unsigned int i = 28; i; i -= 4, idx >>= 4, mAptr += 4) {
    mm256_xor_mask_region(&cval[0], mAptr + 0, _mm256_set1_epi64x(-(idx & 1)), 1);
    mm256_xor_mask_region(&cval[1], mAptr + 1, _mm256_set1_epi64x(-((idx >> 1) & 1)), 1);
    mm256_xor_mask_region(&cval[0], mAptr + 2, _mm256_set1_epi64x(-((idx >> 2) & 1)), 1);
    mm256_xor_mask_region(&cval[1], mAptr + 3, _mm256_set1_epi64x(-((idx >> 3) & 1)), 1);
  }
  mm256_xor_mask_region(&cval[0], mAptr + 0, _mm256_set1_epi64x(-(idx & 1)), 1);
  mm256_xor_mask_region(&cval[1], mAptr + 1, _mm256_set1_epi64x(-((idx >> 1) & 1)), 1);
  *mcptr ^= _mm256_xor_si256(cval[0], cval[1]);
}

ATTR_TARGET("avx2")
void mzd_mul_v_avx_30_256(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  word const* vptr     = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  __m256i* mcptr       = ASSUME_ALIGNED(FIRST_ROW(c), alignof(__m256i));
  __m256i const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(__m256i));

  __m256i cval[2] ATTR_ALIGNED(alignof(__m256i)) = {_mm256_setzero_si256(), _mm256_setzero_si256()};
  word idx = vptr[3] >> 34;
  // do 7x4 and then 2 extra to get 30
  for (unsigned int i = 28; i; i -= 4, idx >>= 4, mAptr += 4) {
    mm256_xor_mask_region(&cval[0], mAptr + 0, _mm256_set1_epi64x(-(idx & 1)), 1);
    mm256_xor_mask_region(&cval[1], mAptr + 1, _mm256_set1_epi64x(-((idx >> 1) & 1)), 1);
    mm256_xor_mask_region(&cval[0], mAptr + 2, _mm256_set1_epi64x(-((idx >> 2) & 1)), 1);
    mm256_xor_mask_region(&cval[1], mAptr + 3, _mm256_set1_epi64x(-((idx >> 3) & 1)), 1);
  }
  mm256_xor_mask_region(&cval[0], mAptr + 0, _mm256_set1_epi64x(-(idx & 1)), 1);
  mm256_xor_mask_region(&cval[1], mAptr + 1, _mm256_set1_epi64x(-((idx >> 1) & 1)), 1);
  *mcptr ^= _mm256_xor_si256(cval[0], cval[1]);
}

ATTR_TARGET("avx2")
void mzd_mul_v_avx_3_128(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
    (void)c, (void)v, (void)A;
    assert(false && "3_128 not faster for AVX, use sse_3_128");
    exit(-1);
}

ATTR_TARGET("avx2")
void mzd_mul_v_avx_3_192(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  word const* vptr     = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  __m256i* mcptr       = ASSUME_ALIGNED(FIRST_ROW(c), alignof(__m256i));
  __m256i const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(__m256i));

  __m256i cval[2] ATTR_ALIGNED(alignof(__m256i)) = {_mm256_setzero_si256(), _mm256_setzero_si256()};
  word idx = vptr[2] >> 61;
  mm256_xor_mask_region(&cval[0], mAptr + 0, _mm256_set1_epi64x(-(idx & 1)), 1);
  mm256_xor_mask_region(&cval[1], mAptr + 1, _mm256_set1_epi64x(-((idx >> 1) & 1)), 1);
  mm256_xor_mask_region(&cval[0], mAptr + 2, _mm256_set1_epi64x(-((idx >> 2) & 1)), 1);
  *mcptr ^= _mm256_xor_si256(cval[0], cval[1]);
}

ATTR_TARGET("avx2")
void mzd_mul_v_avx_3_256(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  word const* vptr     = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  __m256i* mcptr       = ASSUME_ALIGNED(FIRST_ROW(c), alignof(__m256i));
  __m256i const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(__m256i));

  __m256i cval[2] ATTR_ALIGNED(alignof(__m256i)) = {_mm256_setzero_si256(), _mm256_setzero_si256()};
  word idx = vptr[3] >> 61;
  mm256_xor_mask_region(&cval[0], mAptr + 0, _mm256_set1_epi64x(-(idx & 1)), 1);
  mm256_xor_mask_region(&cval[1], mAptr + 1, _mm256_set1_epi64x(-((idx >> 1) & 1)), 1);
  mm256_xor_mask_region(&cval[0], mAptr + 2, _mm256_set1_epi64x(-((idx >> 2) & 1)), 1);
  *mcptr ^= _mm256_xor_si256(cval[0], cval[1]);
}

// Standard multiplication using AVX, slower than 226_30_popcnt without AVX
//ATTR_TARGET("avx2")
//void mzd_mul_v_avx_226_30(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
//  word const* vptr     = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
//  word* cptr                   = ASSUME_ALIGNED(FIRST_ROW(c), 32);
//  __m256i const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(__m256i));
//
//  __m256i cval[2] ATTR_ALIGNED(alignof(__m256i)) = {_mm256_setzero_si256(), _mm256_setzero_si256()};
//  // do 3x2x4x8 and then 30 extra to get 226
//  for (unsigned int w = 3; w; --w, ++vptr) {
//    word idx = *vptr;
//    for (unsigned int i = sizeof(word)*8; i; i -= 32, idx >>= 32, mAptr += 4) {
//      mm256_xor_mask_region(&cval[0], mAptr + 0, _mm256_set_epi32(-((idx >>  0) & 1), -((idx >>  1) & 1), -((idx >>  2) & 1), -((idx >>  3) & 1), -((idx >>  4) & 1), -((idx >>  5) & 1), -((idx >>  6) & 1), -((idx >>  7) & 1)), 1);
//      mm256_xor_mask_region(&cval[1], mAptr + 1, _mm256_set_epi32(-((idx >>  8) & 1), -((idx >>  9) & 1), -((idx >> 10) & 1), -((idx >> 11) & 1), -((idx >> 12) & 1), -((idx >> 13) & 1), -((idx >> 14) & 1), -((idx >> 15) & 1)), 1);
//      mm256_xor_mask_region(&cval[0], mAptr + 2, _mm256_set_epi32(-((idx >> 16) & 1), -((idx >> 17) & 1), -((idx >> 18) & 1), -((idx >> 19) & 1), -((idx >> 20) & 1), -((idx >> 21) & 1), -((idx >> 22) & 1), -((idx >> 23) & 1)), 1);
//      mm256_xor_mask_region(&cval[1], mAptr + 3, _mm256_set_epi32(-((idx >> 24) & 1), -((idx >> 25) & 1), -((idx >> 26) & 1), -((idx >> 27) & 1), -((idx >> 28) & 1), -((idx >> 29) & 1), -((idx >> 30) & 1), -((idx >> 31) & 1)), 1);
//    }
//  }
//  word idx = vptr[3];
//  mm256_xor_mask_region(&cval[0], mAptr + 0, _mm256_set_epi32(-((idx >>  0) & 1), -((idx >>  1) & 1), -((idx >>  2) & 1), -((idx >>  3) & 1), -((idx >>  4) & 1), -((idx >>  5) & 1), -((idx >>  6) & 1), -((idx >>  7) & 1)), 1);
//  mm256_xor_mask_region(&cval[1], mAptr + 1, _mm256_set_epi32(-((idx >>  8) & 1), -((idx >>  9) & 1), -((idx >> 10) & 1), -((idx >> 11) & 1), -((idx >> 12) & 1), -((idx >> 13) & 1), -((idx >> 14) & 1), -((idx >> 15) & 1)), 1);
//  mm256_xor_mask_region(&cval[0], mAptr + 2, _mm256_set_epi32(-((idx >> 16) & 1), -((idx >> 17) & 1), -((idx >> 18) & 1), -((idx >> 19) & 1), -((idx >> 20) & 1), -((idx >> 21) & 1), -((idx >> 22) & 1), -((idx >> 23) & 1)), 1);
//  mm256_xor_mask_region(&cval[1], mAptr + 3, _mm256_set_epi32(-((idx >> 24) & 1), -((idx >> 25) & 1), -((idx >> 26) & 1), -((idx >> 27) & 1), -((idx >> 28) & 1), -((idx >> 29) & 1), 0, 0), 1);
//  cval[0] = _mm256_xor_si256(cval[0], cval[1]);
//  word result =   _mm256_extract_epi32(cval[0], 0) ^ _mm256_extract_epi32(cval[0], 1) ^
//                  _mm256_extract_epi32(cval[0], 2) ^ _mm256_extract_epi32(cval[0], 3) ^
//                  _mm256_extract_epi32(cval[0], 4) ^ _mm256_extract_epi32(cval[0], 5) ^
//                  _mm256_extract_epi32(cval[0], 6) ^ _mm256_extract_epi32(cval[0], 7);
////  printf("0x%016lX\n", result);
//  cptr[3] &= WORD_C(0x00000003FFFFFFFF); //clear nl part
//  cptr[3] |= result << 32;
//}

// Multiplication using AVX & popcnt, slower than 226_30_popcnt without AVX
ATTR_TARGET("avx2")
void mzd_mul_v_avx_226_30_popcnt(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* At) {
  __m256i const* vptr          = ASSUME_ALIGNED(CONST_FIRST_ROW(v), alignof(__m256i));
  __m256i const* mAptr         = ASSUME_ALIGNED(CONST_FIRST_ROW(At), alignof(__m256i));
  word* cptr                   = ASSUME_ALIGNED(FIRST_ROW(c), 32);
  cptr[3] &= WORD_C(0x00000003FFFFFFFF); //clear nl part

  for(unsigned i = 30; i; --i, mAptr++) {
    __m256i cnt = _mm256_and_si256(*vptr, *mAptr);
    word popcnt = __builtin_popcountll(_mm256_extract_epi64(cnt, 0)) +
                  __builtin_popcountll(_mm256_extract_epi64(cnt, 1)) +
                  __builtin_popcountll(_mm256_extract_epi64(cnt, 2)) +
                  __builtin_popcountll(_mm256_extract_epi64(cnt, 3));
    cptr[3] |= (popcnt & WORD_C(0x1)) << (64-i);
  }
}

ATTR_TARGET("bmi2")
void mzd_shuffle_pext_30(mzd_local_t* x, const word mask)  {
  word a = _pext_u64(CONST_FIRST_ROW(x)[x->width - 1], mask) << (34);
  FIRST_ROW(x)[x->width - 1] = a | _pext_u64(CONST_FIRST_ROW(x)[x->width - 1], ~(mask));
}
ATTR_TARGET("bmi2")
void mzd_shuffle_pext_3(mzd_local_t* x, const word mask)  {
  word a = _pext_u64(CONST_FIRST_ROW(x)[x->width - 1], mask) << (61);
  FIRST_ROW(x)[x->width - 1] = a | _pext_u64(CONST_FIRST_ROW(x)[x->width - 1], ~(mask));
}
#endif

#if defined(WITH_NEON)
void mzd_mul_v_neon_30_128(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  word const* vptr     = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  uint32x4_t* mcptr       = ASSUME_ALIGNED(FIRST_ROW(c), alignof(uint32x4_t));
  uint32x4_t const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(uint32x4_t));

  uint32x4_t cval[2] ATTR_ALIGNED(alignof(uint32x4_t)) = {vmovq_n_u32(0), vmovq_n_u32(0)};
  word idx = vptr[1] >> 34;
  for (unsigned int i = 28; i; i -= 4, idx >>= 4, mAptr += 4) {
    mm128_xor_mask_region(&cval[0], mAptr + 0, vreinterpretq_u32_u64(vdupq_n_u64(-(idx & 1))), 1);
    mm128_xor_mask_region(&cval[1], mAptr + 1, vreinterpretq_u32_u64(vdupq_n_u64(-((idx >> 1) & 1))), 1);
    mm128_xor_mask_region(&cval[0], mAptr + 2, vreinterpretq_u32_u64(vdupq_n_u64(-((idx >> 2) & 1))), 1);
    mm128_xor_mask_region(&cval[1], mAptr + 3, vreinterpretq_u32_u64(vdupq_n_u64(-((idx >> 3) & 1))), 1);
  }
  mm128_xor_mask_region(&cval[0], mAptr + 0, vreinterpretq_u32_u64(vdupq_n_u64(-(idx & 1))), 1);
  mm128_xor_mask_region(&cval[1], mAptr + 1, vreinterpretq_u32_u64(vdupq_n_u64(-((idx >> 1) & 1))), 1);
  *mcptr ^= veorq_u32(cval[0], cval[1]);
}

void mzd_mul_v_neon_30_192(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  word const* vptr     = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  uint32x4_t* mcptr       = ASSUME_ALIGNED(FIRST_ROW(c), alignof(uint32x4_t));
  uint32x4_t const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(uint32x4_t));

  uint32x4_t cval[4] ATTR_ALIGNED(alignof(uint32x4_t)) = {vmovq_n_u32(0), vmovq_n_u32(0),
                                                          vmovq_n_u32(0), vmovq_n_u32(0)};
  word idx = vptr[2] >> 34;
  for (unsigned int i = 30; i; i -= 2, idx >>= 2, mAptr += 4) {
    mm128_xor_mask_region(&cval[0], mAptr + 0, vreinterpretq_u32_u64(vdupq_n_u64(-(idx & 1))), 2);
    mm128_xor_mask_region(&cval[2], mAptr + 2, vreinterpretq_u32_u64(vdupq_n_u64(-((idx >> 1) & 1))), 2);
  }
  mcptr[0] ^= veorq_u32(cval[0], cval[2]);
  mcptr[1] ^= veorq_u32(cval[1], cval[3]);
}

void mzd_mul_v_neon_30_256(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  word const* vptr     = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  uint32x4_t* mcptr       = ASSUME_ALIGNED(FIRST_ROW(c), alignof(uint32x4_t));
  uint32x4_t const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(uint32x4_t));

  uint32x4_t cval[4] ATTR_ALIGNED(alignof(uint32x4_t)) = {vmovq_n_u32(0), vmovq_n_u32(0),
                                                          vmovq_n_u32(0), vmovq_n_u32(0)};
  word idx = vptr[3] >> 34;
  for (unsigned int i = 30; i; i -= 2, idx >>= 2, mAptr += 4) {
    mm128_xor_mask_region(&cval[0], mAptr + 0, vreinterpretq_u32_u64(vdupq_n_u64(-(idx & 1))), 2);
    mm128_xor_mask_region(&cval[2], mAptr + 2, vreinterpretq_u32_u64(vdupq_n_u64(-((idx >> 1) & 1))), 2);
  }
  mcptr[0] ^= veorq_u32(cval[0], cval[2]);
  mcptr[1] ^= veorq_u32(cval[1], cval[3]);
}

void mzd_mul_v_neon_3_128(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  word const* vptr     = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  uint32x4_t* mcptr       = ASSUME_ALIGNED(FIRST_ROW(c), alignof(uint32x4_t));
  uint32x4_t const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(uint32x4_t));

  uint32x4_t cval[2] ATTR_ALIGNED(alignof(uint32x4_t)) = {vmovq_n_u32(0), vmovq_n_u32(0)};
  word idx = vptr[1] >> 61;
  mm128_xor_mask_region(&cval[0], mAptr + 0, vreinterpretq_u32_u64(vdupq_n_u64(-(idx & 1))), 1);
  mm128_xor_mask_region(&cval[1], mAptr + 1, vreinterpretq_u32_u64(vdupq_n_u64(-((idx >> 1) & 1))), 1);
  mm128_xor_mask_region(&cval[0], mAptr + 2, vreinterpretq_u32_u64(vdupq_n_u64(-((idx >> 2) & 1))), 1);
  *mcptr ^= veorq_u32(cval[0], cval[1]);
}

void mzd_mul_v_neon_3_192(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  word const* vptr     = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  uint32x4_t* mcptr       = ASSUME_ALIGNED(FIRST_ROW(c), alignof(uint32x4_t));
  uint32x4_t const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(uint32x4_t));

  uint32x4_t cval[4] ATTR_ALIGNED(alignof(uint32x4_t)) = {vmovq_n_u32(0), vmovq_n_u32(0),
                                                          vmovq_n_u32(0), vmovq_n_u32(0)};
  word idx = vptr[2] >> 61;
  mm128_xor_mask_region(&cval[0], mAptr + 0, vreinterpretq_u32_u64(vdupq_n_u64(-(idx & 1))), 2);
  mm128_xor_mask_region(&cval[2], mAptr + 2, vreinterpretq_u32_u64(vdupq_n_u64(-((idx >> 1) & 1))), 2);
  mm128_xor_mask_region(&cval[0], mAptr + 4, vreinterpretq_u32_u64(vdupq_n_u64(-((idx >> 2) & 1))), 2);
  mcptr[0] ^= veorq_u32(cval[0], cval[2]);
  mcptr[1] ^= veorq_u32(cval[1], cval[3]);
}

void mzd_mul_v_neon_3_256(mzd_local_t* c, mzd_local_t const* v, mzd_local_t const* A) {
  word const* vptr     = ASSUME_ALIGNED(CONST_FIRST_ROW(v), 32);
  uint32x4_t* mcptr       = ASSUME_ALIGNED(FIRST_ROW(c), alignof(uint32x4_t));
  uint32x4_t const* mAptr = ASSUME_ALIGNED(CONST_FIRST_ROW(A), alignof(uint32x4_t));

  uint32x4_t cval[4] ATTR_ALIGNED(alignof(uint32x4_t)) = {vmovq_n_u32(0), vmovq_n_u32(0),
                                                          vmovq_n_u32(0), vmovq_n_u32(0)};
  word idx = vptr[3] >> 61;
  mm128_xor_mask_region(&cval[0], mAptr + 0, vreinterpretq_u32_u64(vdupq_n_u64(-(idx & 1))), 2);
  mm128_xor_mask_region(&cval[2], mAptr + 2, vreinterpretq_u32_u64(vdupq_n_u64(-((idx >> 1) & 1))), 2);
  mm128_xor_mask_region(&cval[0], mAptr + 4, vreinterpretq_u32_u64(vdupq_n_u64(-((idx >> 2) & 1))), 2);

  mcptr[0] ^= veorq_u32(cval[0], cval[2]);
  mcptr[1] ^= veorq_u32(cval[1], cval[3]);
}
#endif

#endif
