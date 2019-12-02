/*
 *  This file is part of the optimized implementation of the Picnic signature scheme.
 *  See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */

#include <assert.h>
#include <string.h>

#include "endian_compat.h"
#include "picnic2_simulate_mul.h"

#if defined(WITH_SSE2)
#define ATTR_TARGET_S128 ATTR_TARGET_SSE2
#else
#define ATTR_TARGET_S128
#endif

#define PACKING_FACTOR 4
#define PARTIES_LOG 4

static const block_t block_masks[] = {
    {{
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
    }},
    {{
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0xffffffffffffffff),
    }},
    {{
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0xffffffffffffffff),
        UINT64_C(0x0000000000000000),
    }},
    {{
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0xffffffffffffffff),
        UINT64_C(0xffffffffffffffff),
    }},
    {{
        UINT64_C(0x0000000000000000),
        UINT64_C(0xffffffffffffffff),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
    }},
    {{
        UINT64_C(0x0000000000000000),
        UINT64_C(0xffffffffffffffff),
        UINT64_C(0x0000000000000000),
        UINT64_C(0xffffffffffffffff),
    }},
    {{
        UINT64_C(0x0000000000000000),
        UINT64_C(0xffffffffffffffff),
        UINT64_C(0xffffffffffffffff),
        UINT64_C(0x0000000000000000),
    }},
    {{
        UINT64_C(0x0000000000000000),
        UINT64_C(0xffffffffffffffff),
        UINT64_C(0xffffffffffffffff),
        UINT64_C(0xffffffffffffffff),
    }},
    {{
        UINT64_C(0xffffffffffffffff),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
    }},
    {{
        UINT64_C(0xffffffffffffffff),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
        UINT64_C(0xffffffffffffffff),
    }},
    {{
        UINT64_C(0xffffffffffffffff),
        UINT64_C(0x0000000000000000),
        UINT64_C(0xffffffffffffffff),
        UINT64_C(0x0000000000000000),
    }},
    {{
        UINT64_C(0xffffffffffffffff),
        UINT64_C(0x0000000000000000),
        UINT64_C(0xffffffffffffffff),
        UINT64_C(0xffffffffffffffff),
    }},
    {{
        UINT64_C(0xffffffffffffffff),
        UINT64_C(0xffffffffffffffff),
        UINT64_C(0x0000000000000000),
        UINT64_C(0x0000000000000000),
    }},
    {{
        UINT64_C(0xffffffffffffffff),
        UINT64_C(0xffffffffffffffff),
        UINT64_C(0x0000000000000000),
        UINT64_C(0xffffffffffffffff),
    }},
    {{
        UINT64_C(0xffffffffffffffff),
        UINT64_C(0xffffffffffffffff),
        UINT64_C(0xffffffffffffffff),
        UINT64_C(0x0000000000000000),
    }},
    {{
        UINT64_C(0xffffffffffffffff),
        UINT64_C(0xffffffffffffffff),
        UINT64_C(0xffffffffffffffff),
        UINT64_C(0xffffffffffffffff),
    }},
};

#if defined(PICNIC_STATIC)
/* transpose a 64x64 bit matrix using Eklundh's algorithm
   this variant assumes that the bit with index 0 is the lsb of byte 0
   e.g., 76543210 fedcba98 ...
 */
void transpose_64_64_lsb(const uint64_t* in, uint64_t* out) {
  static const uint64_t TRANSPOSE_MASKS64[6] = {
      UINT64_C(0x00000000FFFFFFFF), UINT64_C(0x0000FFFF0000FFFF), UINT64_C(0x00FF00FF00FF00FF),
      UINT64_C(0x0F0F0F0F0F0F0F0F), UINT64_C(0x3333333333333333), UINT64_C(0x5555555555555555)};

  uint32_t width = 32, nswaps = 1;
  const uint32_t logn = 6;

  // copy in to out and transpose in-place
  memcpy(out, in, 64 * sizeof(uint64_t));

  for (uint32_t i = 0; i < logn; i++) {
    uint64_t mask     = TRANSPOSE_MASKS64[i];
    uint64_t inv_mask = ~mask;

    for (uint32_t j = 0; j < nswaps; j++) {
      for (uint32_t k = 0; k < width; k++) {
        uint32_t i1 = k + 2 * width * j;
        uint32_t i2 = k + width + 2 * width * j;

        uint64_t t1 = out[i1];
        uint64_t t2 = out[i2];

        out[i1] = (t1 & mask) ^ ((t2 & mask) << width);
        out[i2] = (t2 & inv_mask) ^ ((t1 & inv_mask) >> width);
      }
    }
    nswaps *= 2;
    width /= 2;
  }
}
#endif

/* transpose a 64x64 bit matrix using Eklundh's algorithm
   this variant assumes that the bit with index 0 is the msb of byte 0
   e.g., 01234567 89abcdef ...
 */
static void transpose_64_64_uint64(const uint64_t* in, uint64_t* out) {
  static const uint64_t TRANSPOSE_MASKS64[6] = {
      UINT64_C(0xFFFFFFFF00000000), UINT64_C(0xFFFF0000FFFF0000), UINT64_C(0xFF00FF00FF00FF00),
      UINT64_C(0xF0F0F0F0F0F0F0F0), UINT64_C(0xCCCCCCCCCCCCCCCC), UINT64_C(0xAAAAAAAAAAAAAAAA)};

  uint32_t width = 32, nswaps = 1;
  const uint32_t logn = 6;

  // copy in to out and transpose in-place
  for (uint32_t i = 0; i < 64; i++) {
    out[i] = bswap64(in[i]);
  }

  for (uint32_t i = 0; i < logn; i++) {
    uint64_t mask     = TRANSPOSE_MASKS64[i];
    uint64_t inv_mask = ~mask;

    for (uint32_t j = 0; j < nswaps; j++) {
      for (uint32_t k = 0; k < width; k++) {
        uint32_t i1 = k + 2 * width * j;
        uint32_t i2 = k + width + 2 * width * j;

        uint64_t t1 = out[i1];
        uint64_t t2 = out[i2];

        out[i1] = (t1 & mask) ^ ((t2 & mask) >> width);
        out[i2] = (t2 & inv_mask) ^ ((t1 & inv_mask) << width);
      }
    }
    nswaps *= 2;
    width /= 2;
  }
  for (uint32_t i = 0; i < 64; i++) {
    out[i] = bswap64(out[i]);
  }
}

#if defined(WITH_OPT)
#if defined(WITH_SSE2) || defined(WITH_NEON)
ATTR_TARGET_S128
static inline void memcpy_bswap64_64_s128(word128* out, const word128* in) {
#if defined(WITH_NEON)
  for (size_t i = 0; i < 32; ++i) {
    out[i] = vreinterpretq_u64_u8(vrev64q_u8(vreinterpretq_u8_u64(in[i])));
  }
#elif defined(WITH_SSE2)
  for (size_t i = 0; i < 32; ++i) {
    word128 tmp = _mm_or_si128(_mm_slli_epi16(in[i], 8), _mm_srli_epi16(in[i], 8));
    tmp         = _mm_shufflelo_epi16(tmp, _MM_SHUFFLE(0, 1, 2, 3));
    out[i]      = _mm_shufflehi_epi16(tmp, _MM_SHUFFLE(0, 1, 2, 3));
  }
#endif
}

/* transpose a 64x64 bit matrix using Eklundh's algorithm
   this variant assumes that the bit with index 0 is the msb of byte 0
   e.g., 01234567 89abcdef ...
 */
ATTR_TARGET_S128
static void transpose_64_64_s128(const uint64_t* in, uint64_t* out) {
  static const uint64_t TRANSPOSE_MASKS64[6] = {
      UINT64_C(0xFFFFFFFF00000000), UINT64_C(0xFFFF0000FFFF0000), UINT64_C(0xFF00FF00FF00FF00),
      UINT64_C(0xF0F0F0F0F0F0F0F0), UINT64_C(0xCCCCCCCCCCCCCCCC), UINT64_C(0xAAAAAAAAAAAAAAAA)};

  uint32_t width = 32, nswaps = 1;
  const uint32_t logn = 6;

  // copy in to out and transpose in-place
  word128* out128      = (word128*)out;
  const word128* in128 = (const word128*)in;
  memcpy_bswap64_64_s128(out128, in128);

  for (uint32_t i = 0; i < logn - 1; i++) {
    const word128 mask = mm128_broadcast_u64(TRANSPOSE_MASKS64[i]);
    for (uint32_t j = 0; j < nswaps; j++) {
      for (uint32_t k = 0; k < width; k += 2) {
        // uint32_t i1 = k/2 + width * j;
        // uint32_t i2 = i1 + width / 2;
        uint32_t i1 = k + 2 * width * j;
        uint32_t i2 = k + width + 2 * width * j;

        word128 t1 = out128[i1 / 2];
        word128 t2 = out128[i2 / 2];

        out128[i1 / 2] = mm128_xor(mm128_and(t1, mask), mm128_sr_u64(mm128_and(t2, mask), width));
        out128[i2 / 2] = mm128_xor(mm128_nand(mask, t2), mm128_sl_u64(mm128_nand(mask, t1), width));
      }
    }
    nswaps *= 2;
    width /= 2;
  }
  const uint64_t mask     = TRANSPOSE_MASKS64[5];
  const uint64_t inv_mask = ~mask;
  for (uint32_t j = 0; j < nswaps; j++) {
    for (uint32_t k = 0; k < width; k++) {
      uint32_t i1 = k + 2 * width * j;
      uint32_t i2 = k + width + 2 * width * j;

      uint64_t t1 = out[i1];
      uint64_t t2 = out[i2];

      out[i1] = (t1 & mask) ^ ((t2 & mask) >> width);
      out[i2] = (t2 & inv_mask) ^ ((t1 & inv_mask) << width);
    }
  }

  memcpy_bswap64_64_s128(out128, out128);
}
#endif

#if defined(WITH_AVX2)
ATTR_TARGET_AVX2
static inline void memcpy_bswap64_64_s256(word256* out, const word256* in) {
  const word256 shuffle = _mm256_set_epi8(
      // two bswap64s in first lane
      8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7,
      // two bswap64s in second lane
      8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7);

  for (size_t i = 0; i < 16; i++) {
    out[i] = _mm256_shuffle_epi8(in[i], shuffle);
  }
}

/* transpose a 64x64 bit matrix using Eklundh's algorithm
   this variant assumes that the bit with index 0 is the msb of byte 0
   e.g., 01234567 89abcdef ...
 */
ATTR_TARGET_AVX2
static void transpose_64_64_s256(const uint64_t* in, uint64_t* out) {
  static const uint64_t TRANSPOSE_MASKS64[6] = {
      UINT64_C(0xFFFFFFFF00000000), UINT64_C(0xFFFF0000FFFF0000), UINT64_C(0xFF00FF00FF00FF00),
      UINT64_C(0xF0F0F0F0F0F0F0F0), UINT64_C(0xCCCCCCCCCCCCCCCC), UINT64_C(0xAAAAAAAAAAAAAAAA)};

  uint32_t width = 32, nswaps = 1;
  static const uint32_t logn = 6;

  const word256* in256 = (const word256*)in;
  word256* out256      = (word256*)out;

  // copy in to out and swap bytes
  memcpy_bswap64_64_s256(out256, in256);

  for (uint32_t i = 0; i < logn - 2; i++) {
    const word256 mask = _mm256_set1_epi64x(TRANSPOSE_MASKS64[i]);
    for (uint32_t j = 0; j < nswaps; j++) {
      for (uint32_t k = 0; k < width; k += 4) {
        uint32_t i1 = k + 2 * width * j;
        uint32_t i2 = k + width + 2 * width * j;

        word256 t1 = out256[i1 / 4];
        word256 t2 = out256[i2 / 4];

        out256[i1 / 4] =
            mm256_xor(mm256_and(t1, mask), _mm256_srli_epi64(mm256_and(t2, mask), width));
        out256[i2 / 4] =
            mm256_xor(mm256_nand(mask, t2), _mm256_slli_epi64(mm256_nand(mask, t1), width));
      }
    }
    nswaps *= 2;
    width /= 2;
  }
  {
    word128* out128    = (word128*)out;
    const word128 mask = mm128_broadcast_u64(TRANSPOSE_MASKS64[4]);

    for (uint32_t j = 0; j < nswaps; j++) {
      for (uint32_t k = 0; k < width; k += 2) {
        uint32_t i1 = k + 2 * width * j;
        uint32_t i2 = k + width + 2 * width * j;

        word128 t1 = out128[i1 / 2];
        word128 t2 = out128[i2 / 2];

        out128[i1 / 2] = mm128_xor(mm128_and(t1, mask), mm128_sr_u64(mm128_and(t2, mask), width));
        out128[i2 / 2] = mm128_xor(mm128_nand(mask, t2), mm128_sl_u64(mm128_nand(mask, t1), width));
      }
    }
    nswaps *= 2;
    width /= 2;
  }

  const uint64_t mask     = TRANSPOSE_MASKS64[5];
  const uint64_t inv_mask = ~mask;
  for (uint32_t j = 0; j < nswaps; j++) {
    for (uint32_t k = 0; k < width; k++) {
      uint32_t i1 = k + 2 * width * j;
      uint32_t i2 = k + width + 2 * width * j;

      uint64_t t1 = out[i1];
      uint64_t t2 = out[i2];

      out[i1] = (t1 & mask) ^ ((t2 & mask) >> width);
      out[i2] = (t2 & inv_mask) ^ ((t1 & inv_mask) << width);
    }
  }

  memcpy_bswap64_64_s256(out256, out256);
}
#endif
#endif

void transpose_64_64(const uint64_t* in, uint64_t* out) {
#if defined(WITH_OPT)
#if defined(WITH_AVX2)
  if (CPU_SUPPORTS_AVX2) {
    transpose_64_64_s256(in, out);
    return;
  }
#endif
#if defined(WITH_SSE2) || defined(WITH_NEON)
  if (CPU_SUPPORTS_NEON || CPU_SUPPORTS_SSE2) {
    transpose_64_64_s128(in, out);
    return;
  }
#endif
#endif
  transpose_64_64_uint64(in, out);
}

uint64_t tapesToParityOfWord(randomTape_t* tapes, uint8_t without_last) {
  uint64_t shares;
  const size_t wordBits = sizeof(uint64_t) * 8;

  if (tapes->pos % wordBits == 0) {
    tapes->buffer[0] = 0;
    for (size_t i = 0; i < tapes->nTapes - 1; i++) {
      tapes->buffer[0] ^= ((uint64_t*)tapes->tape[i])[tapes->pos / wordBits];
    }
    tapes->buffer[1] = tapes->buffer[0];
    tapes->buffer[0] ^= ((uint64_t*)tapes->tape[tapes->nTapes - 1])[tapes->pos / wordBits];
  }

  shares = getBit((uint8_t*)&tapes->buffer[without_last ? 1 : 0], tapes->pos % wordBits);
  tapes->pos++;
  return shares;
}

uint64_t tapesToWord(randomTape_t* tapes) {
  uint64_t shares;
  const size_t wordBits = sizeof(uint64_t) * 8;

  if (tapes->pos % wordBits == 0) {
    for (size_t j = 0; j < PACKING_FACTOR; j++) {
      for (size_t i = 0; i < tapes->nTapes; i++) {
        tapes->buffer[j * tapes->nTapes + i] = ((uint64_t*)tapes[j].tape[i])[tapes->pos / wordBits];
      }
    }
    transpose_64_64(tapes->buffer, tapes->buffer);
  }

  shares = tapes->buffer[tapes->pos % wordBits];
  for (size_t j = 0; j < PACKING_FACTOR; j++) {
    tapes[j].pos++;
  }
  return shares;
}

void reconstructShares(uint32_t* output, shares_t* shares) {
  for (size_t i = 0; i < shares->numWords; i++) {
    setBit((uint8_t*)output, i, parity64_uint64(shares->shares[i]));
  }
}

void reconstructSharesN(uint32_t** output, shares_t* shares) {
  for (size_t i = 0; i < shares->numWords; i++) {
    uint64_t a = shares->shares[i];
    for (size_t k = 0; k < PARTIES_LOG; k++) {
      a ^= a >> (1 << (PARTIES_LOG - 1 - k));
    }
    for (size_t k = 0; k < PACKING_FACTOR; k++) {
      setBit((uint8_t*)output[k], i, a & 1);
      a >>= (64 / PACKING_FACTOR);
    }
  }
}

void xor_word_array(uint32_t* out, const uint32_t* in1, const uint32_t* in2, uint32_t length) {
  for (uint32_t i = 0; i < length; i++) {
    out[i] = in1[i] ^ in2[i];
  }
}

void copyShares(shares_t* dst, shares_t* src) {
  assert(dst->numWords == src->numWords);
  memcpy(dst->shares, src->shares, dst->numWords * sizeof(dst->shares[0]));
}

void mpc_matrix_mul_uint64_128(mzd_local_t** output, mzd_local_t** vec, const mzd_local_t* matrix,
                               shares_t* mask_shares) {
  const uint32_t rowstride = (128) / 8;
  shares_t* tmp_mask       = allocateShares(mask_shares->numWords);

  for (size_t i = 0; i < 128; i++) {
    const uint64_t mask_share = mask_shares->shares[128 - 1 - i];

    for (uint32_t j = 0; j < 128; j += 8) {
      uint8_t matrix_byte  = ((const uint8_t*)matrix->w64)[(i * rowstride) + (128 - 1 - j) / 8];
      const block_t* mask1 = &block_masks[(matrix_byte >> 4) & 0xF];
      const block_t* mask2 = &block_masks[(matrix_byte >> 0) & 0xF];

      tmp_mask->shares[j + 0] ^= mask_share & mask1->w64[0];
      tmp_mask->shares[j + 1] ^= mask_share & mask1->w64[1];
      tmp_mask->shares[j + 2] ^= mask_share & mask1->w64[2];
      tmp_mask->shares[j + 3] ^= mask_share & mask1->w64[3];
      tmp_mask->shares[j + 4] ^= mask_share & mask2->w64[0];
      tmp_mask->shares[j + 5] ^= mask_share & mask2->w64[1];
      tmp_mask->shares[j + 6] ^= mask_share & mask2->w64[2];
      tmp_mask->shares[j + 7] ^= mask_share & mask2->w64[3];
    }
  }
  for (uint32_t k = 0; k < PACKING_FACTOR; k++) {
    mzd_local_t tmp;
    mzd_copy_uint64_128(&tmp, vec[k]);
    mzd_mul_v_uint64_128(output[k], &tmp, matrix);
  }

  copyShares(mask_shares, tmp_mask);
  freeShares(tmp_mask);
}

void mpc_matrix_mul_uint64_192(mzd_local_t** output, mzd_local_t** vec, const mzd_local_t* matrix,
                               shares_t* mask_shares) {
  const uint32_t rowstride = (256) / 8;
  shares_t* tmp_mask       = allocateShares(mask_shares->numWords);

  for (size_t i = 0; i < 192; i++) {
    const uint64_t mask_share = mask_shares->shares[192 - 1 - i];

    for (uint32_t j = 0; j < 192; j += 8) {
      uint8_t matrix_byte = ((const uint8_t*)matrix->w64)[(i * rowstride) + (192 - 1 - j) / 8];

      const block_t* mask1 = &block_masks[(matrix_byte >> 4) & 0xF];
      const block_t* mask2 = &block_masks[(matrix_byte >> 0) & 0xF];

      tmp_mask->shares[j + 0] ^= mask_share & mask1->w64[0];
      tmp_mask->shares[j + 1] ^= mask_share & mask1->w64[1];
      tmp_mask->shares[j + 2] ^= mask_share & mask1->w64[2];
      tmp_mask->shares[j + 3] ^= mask_share & mask1->w64[3];
      tmp_mask->shares[j + 4] ^= mask_share & mask2->w64[0];
      tmp_mask->shares[j + 5] ^= mask_share & mask2->w64[1];
      tmp_mask->shares[j + 6] ^= mask_share & mask2->w64[2];
      tmp_mask->shares[j + 7] ^= mask_share & mask2->w64[3];
    }
  }
  for (uint32_t k = 0; k < PACKING_FACTOR; k++) {
    mzd_local_t tmp;
    mzd_copy_uint64_192(&tmp, vec[k]);
    mzd_mul_v_uint64_192(output[k], &tmp, matrix);
  }

  copyShares(mask_shares, tmp_mask);
  freeShares(tmp_mask);
}
void mpc_matrix_mul_uint64_256(mzd_local_t** output, mzd_local_t** vec, const mzd_local_t* matrix,
                               shares_t* mask_shares) {
  const uint32_t rowstride = (256) / 8;
  shares_t* tmp_mask       = allocateShares(mask_shares->numWords);

  for (size_t i = 0; i < 256; i++) {
    const uint64_t mask_share = mask_shares->shares[256 - 1 - i];

    for (uint32_t j = 0; j < 256; j += 8) {
      uint8_t matrix_byte = ((const uint8_t*)matrix->w64)[(i * rowstride) + (256 - 1 - j) / 8];

      const block_t* mask1 = &block_masks[(matrix_byte >> 4) & 0xF];
      const block_t* mask2 = &block_masks[(matrix_byte >> 0) & 0xF];

      tmp_mask->shares[j + 0] ^= mask_share & mask1->w64[0];
      tmp_mask->shares[j + 1] ^= mask_share & mask1->w64[1];
      tmp_mask->shares[j + 2] ^= mask_share & mask1->w64[2];
      tmp_mask->shares[j + 3] ^= mask_share & mask1->w64[3];
      tmp_mask->shares[j + 4] ^= mask_share & mask2->w64[0];
      tmp_mask->shares[j + 5] ^= mask_share & mask2->w64[1];
      tmp_mask->shares[j + 6] ^= mask_share & mask2->w64[2];
      tmp_mask->shares[j + 7] ^= mask_share & mask2->w64[3];
    }
  }
  for (uint32_t k = 0; k < PACKING_FACTOR; k++) {
    mzd_local_t tmp;
    mzd_copy_uint64_256(&tmp, vec[k]);
    mzd_mul_v_uint64_256(output[k], &tmp, matrix);
  }

  copyShares(mask_shares, tmp_mask);
  freeShares(tmp_mask);
}

#if defined(WITH_OPT)
#if defined(WITH_SSE2) || defined(WITH_NEON)
ATTR_TARGET_S128
void mpc_matrix_mul_s128_128(mzd_local_t** output, mzd_local_t** vec, const mzd_local_t* matrix,
                             shares_t* mask_shares) {
  const uint32_t rowstride = (128) / 8;
  shares_t* tmp_mask       = allocateShares(mask_shares->numWords);

  for (size_t i = 0; i < 128; i++) {
    const uint64_t mask_share = mask_shares->shares[128 - 1 - i];
    const block_t mask_share2 = {{mask_share, mask_share, mask_share, mask_share}};
    word128 mask1, mask2, mask3, mask4;

    word128* tmp_mask_block = (word128*)tmp_mask->shares;

    for (uint32_t j = 0; j < 128; j += 8, tmp_mask_block += 4) {
      uint8_t matrix_byte = ((const uint8_t*)matrix->w64)[(i * rowstride) + (128 - 1 - j) / 8];

      mask1 = block_masks[(matrix_byte >> 4) & 0xf].w128[0];
      mask2 = block_masks[(matrix_byte >> 4) & 0xf].w128[1];
      mask3 = block_masks[(matrix_byte >> 0) & 0xf].w128[0];
      mask4 = block_masks[(matrix_byte >> 0) & 0xf].w128[1];

      tmp_mask_block[0] = mm128_xor_mask(tmp_mask_block[0], mask_share2.w128[0], mask1);
      tmp_mask_block[1] = mm128_xor_mask(tmp_mask_block[1], mask_share2.w128[0], mask2);
      tmp_mask_block[2] = mm128_xor_mask(tmp_mask_block[2], mask_share2.w128[0], mask3);
      tmp_mask_block[3] = mm128_xor_mask(tmp_mask_block[3], mask_share2.w128[0], mask4);
    }
  }
  for (uint32_t k = 0; k < PACKING_FACTOR; k++) {
    mzd_mul_v_s128_128(output[k], vec[k], matrix);
  }

  copyShares(mask_shares, tmp_mask);
  freeShares(tmp_mask);
}

ATTR_TARGET_S128
void mpc_matrix_mul_s128_192(mzd_local_t** output, mzd_local_t** vec, const mzd_local_t* matrix,
                             shares_t* mask_shares) {
  const uint32_t rowstride = (256) / 8;
  shares_t* tmp_mask       = allocateShares(mask_shares->numWords);

  for (size_t i = 0; i < 192; i++) {
    const uint64_t mask_share = mask_shares->shares[192 - 1 - i];
    const block_t mask_share2 = {{mask_share, mask_share, mask_share, mask_share}};
    word128 mask1, mask2, mask3, mask4;

    word128* tmp_mask_block = (word128*)tmp_mask->shares;

    for (uint32_t j = 0; j < 192; j += 8, tmp_mask_block += 4) {
      uint8_t matrix_byte = ((const uint8_t*)matrix->w64)[(i * rowstride) + (192 - 1 - j) / 8];

      mask1 = block_masks[(matrix_byte >> 4) & 0xf].w128[0];
      mask2 = block_masks[(matrix_byte >> 4) & 0xf].w128[1];
      mask3 = block_masks[(matrix_byte >> 0) & 0xf].w128[0];
      mask4 = block_masks[(matrix_byte >> 0) & 0xf].w128[1];

      tmp_mask_block[0] = mm128_xor_mask(tmp_mask_block[0], mask_share2.w128[0], mask1);
      tmp_mask_block[1] = mm128_xor_mask(tmp_mask_block[1], mask_share2.w128[0], mask2);
      tmp_mask_block[2] = mm128_xor_mask(tmp_mask_block[2], mask_share2.w128[0], mask3);
      tmp_mask_block[3] = mm128_xor_mask(tmp_mask_block[3], mask_share2.w128[0], mask4);
    }
  }
  for (uint32_t k = 0; k < PACKING_FACTOR; k++) {
    mzd_mul_v_s128_192(output[k], vec[k], matrix);
  }

  copyShares(mask_shares, tmp_mask);
  freeShares(tmp_mask);
}

ATTR_TARGET_S128
void mpc_matrix_mul_s128_256(mzd_local_t** output, mzd_local_t** vec, const mzd_local_t* matrix,
                             shares_t* mask_shares) {
  const uint32_t rowstride = (256) / 8;
  shares_t* tmp_mask       = allocateShares(mask_shares->numWords);

  for (size_t i = 0; i < 256; i++) {
    const uint64_t mask_share = mask_shares->shares[256 - 1 - i];
    const block_t mask_share2 = {{mask_share, mask_share, mask_share, mask_share}};
    word128 mask1, mask2, mask3, mask4;

    word128* tmp_mask_block = (word128*)tmp_mask->shares;

    for (uint32_t j = 0; j < 256; j += 8, tmp_mask_block += 4) {
      uint8_t matrix_byte = ((const uint8_t*)matrix->w64)[(i * rowstride) + (256 - 1 - j) / 8];

      mask1 = block_masks[(matrix_byte >> 4) & 0xf].w128[0];
      mask2 = block_masks[(matrix_byte >> 4) & 0xf].w128[1];
      mask3 = block_masks[(matrix_byte >> 0) & 0xf].w128[0];
      mask4 = block_masks[(matrix_byte >> 0) & 0xf].w128[1];

      tmp_mask_block[0] = mm128_xor_mask(tmp_mask_block[0], mask_share2.w128[0], mask1);
      tmp_mask_block[1] = mm128_xor_mask(tmp_mask_block[1], mask_share2.w128[0], mask2);
      tmp_mask_block[2] = mm128_xor_mask(tmp_mask_block[2], mask_share2.w128[0], mask3);
      tmp_mask_block[3] = mm128_xor_mask(tmp_mask_block[3], mask_share2.w128[0], mask4);
    }
  }
  for (uint32_t k = 0; k < PACKING_FACTOR; k++) {
    mzd_mul_v_s128_256(output[k], vec[k], matrix);
  }

  copyShares(mask_shares, tmp_mask);
  freeShares(tmp_mask);
}
#endif

#if defined(WITH_AVX2)
ATTR_TARGET_AVX2
void mpc_matrix_mul_s256_128(mzd_local_t** output, mzd_local_t** vec, const mzd_local_t* matrix,
                             shares_t* mask_shares) {
  const uint32_t rowstride = (128) / 8;
  shares_t* tmp_mask       = allocateShares(mask_shares->numWords);

  for (size_t i = 0; i < 128; i++) {
    const uint64_t mask_share = mask_shares->shares[128 - 1 - i];
    const word256 mask_share2 = _mm256_set1_epi64x(mask_share);
    word256 mask1, mask2;

    word256* tmp_mask_block = (word256*)tmp_mask->shares;

    for (uint32_t j = 0; j < 128; j += 8, tmp_mask_block += 2) {
      uint8_t matrix_byte = ((const uint8_t*)matrix)[(i * rowstride) + (128 - 1 - j) / 8];
      mask1               = block_masks[(matrix_byte >> 4) & 0xf].w256;
      mask2               = block_masks[(matrix_byte >> 0) & 0xf].w256;

      tmp_mask_block[0] = mm256_xor_mask(tmp_mask_block[0], mask_share2, mask1);
      tmp_mask_block[1] = mm256_xor_mask(tmp_mask_block[1], mask_share2, mask2);
    }
  }
  for (uint32_t k = 0; k < PACKING_FACTOR; k++) {
    mzd_mul_v_s256_128(output[k], vec[k], matrix);
  }

  copyShares(mask_shares, tmp_mask);
  freeShares(tmp_mask);
}

ATTR_TARGET_AVX2
void mpc_matrix_mul_s256_192(mzd_local_t** output, mzd_local_t** vec, const mzd_local_t* matrix,
                             shares_t* mask_shares) {
  const uint32_t rowstride = (256) / 8;
  shares_t* tmp_mask       = allocateShares(mask_shares->numWords);

  for (size_t i = 0; i < 192; i++) {
    const uint64_t mask_share = mask_shares->shares[192 - 1 - i];
    const word256 mask_share2 = _mm256_set1_epi64x(mask_share);
    word256 mask1, mask2;

    word256* tmp_mask_block = (word256*)tmp_mask->shares;

    for (uint32_t j = 0; j < 192; j += 8, tmp_mask_block += 2) {
      uint8_t matrix_byte = ((const uint8_t*)matrix->w64)[(i * rowstride) + (192 - 1 - j) / 8];
      mask1               = block_masks[(matrix_byte >> 4) & 0xf].w256;
      mask2               = block_masks[(matrix_byte >> 0) & 0xf].w256;

      tmp_mask_block[0] = mm256_xor_mask(tmp_mask_block[0], mask_share2, mask1);
      tmp_mask_block[1] = mm256_xor_mask(tmp_mask_block[1], mask_share2, mask2);
    }
  }
  for (uint32_t k = 0; k < PACKING_FACTOR; k++) {
    mzd_mul_v_s256_192(output[k], vec[k], matrix);
  }

  copyShares(mask_shares, tmp_mask);
  freeShares(tmp_mask);
}

ATTR_TARGET_AVX2
void mpc_matrix_mul_s256_256(mzd_local_t** output, mzd_local_t** vec, const mzd_local_t* matrix,
                             shares_t* mask_shares) {
  const uint32_t rowstride = (256) / 8;
  shares_t* tmp_mask       = allocateShares(mask_shares->numWords);

  for (size_t i = 0; i < 256; i++) {
    const uint64_t mask_share = mask_shares->shares[256 - 1 - i];
    const word256 mask_share2 = _mm256_set1_epi64x(mask_share);
    word256 mask1, mask2;

    word256* tmp_mask_block = (word256*)tmp_mask->shares;

    for (uint32_t j = 0; j < 256; j += 8, tmp_mask_block += 2) {
      uint8_t matrix_byte = ((const uint8_t*)matrix->w64)[(i * rowstride) + (256 - 1 - j) / 8];
      mask1               = block_masks[(matrix_byte >> 4) & 0xf].w256;
      mask2               = block_masks[(matrix_byte >> 0) & 0xf].w256;

      tmp_mask_block[0] = mm256_xor_mask(tmp_mask_block[0], mask_share2, mask1);
      tmp_mask_block[1] = mm256_xor_mask(tmp_mask_block[1], mask_share2, mask2);
    }
  }
  for (uint32_t k = 0; k < PACKING_FACTOR; k++) {
    mzd_mul_v_s256_256(output[k], vec[k], matrix);
  }

  copyShares(mask_shares, tmp_mask);
  freeShares(tmp_mask);
}
#endif
#endif
