/*! @file picnic3_impl.c
 *  @brief This is the main file of the signature scheme for the Picnic3
 *  parameter sets.
 *
 *  This file is part of the reference implementation of the Picnic signature scheme.
 *  See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#if !defined(_MSC_VER)
#include <stdalign.h>
#endif

#include "io.h"
#include "picnic3_simulate.h"
#include "picnic3_types.h"

#define PACKING_FACTOR 4
#define PARTIES_LOG 4

/* transpose a 16x16 bit matrix using Eklundh's algorithm
   this variant assumes that the bit with index 0 is the msb of byte 0
   e.g., 01234567 89abcdef ...
 */
void transpose_16_16_uint64(const uint16_t* in, uint16_t* out) {
  static const uint16_t TRANSPOSE_MASKS16[4] = {UINT16_C(0xFF00), UINT16_C(0xF0F0),
                                                UINT16_C(0xCCCC), UINT16_C(0xAAAA)};

  uint32_t width = 8, nswaps = 1;
  const uint32_t logn = 4;

  // copy in to out and transpose in-place
  for (uint32_t i = 0; i < 16; i++) {
    out[i] = htobe16(in[i]);
  }

  for (uint32_t i = 0; i < logn; i++) {
    uint64_t mask     = TRANSPOSE_MASKS16[i];
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
  for (uint32_t i = 0; i < 16; i++) {
    out[i] = be16toh(out[i]);
  }
}

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
    out[i] = htobe64(in[i]);
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
    out[i] = be64toh(out[i]);
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

#if !defined(PICNIC_STATIC)
static
#endif
    void
    transpose_64_64(const uint64_t* in, uint64_t* out) {
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

static uint16_t tapesToWord(randomTape_t* tapes) {
  uint16_t shares;
  const size_t wordBits = sizeof(uint16_t) * 8;

  if (tapes->pos % wordBits == 0) {
    for (size_t i = 0; i < tapes->nTapes; i++) {
      tapes->buffer[i] = ((uint16_t*)tapes->tape[i])[tapes->pos / wordBits];
    }
    transpose_16_16_uint64(tapes->buffer, tapes->buffer);
  }

  shares = tapes->buffer[tapes->pos % wordBits];
  tapes->pos++;
  return shares;
}

static void wordToMsgsNoTranspose(uint16_t w, msgs_t* msgs) {
  ((uint16_t*)msgs->msgs[msgs->pos % 16])[msgs->pos / 16] = w;
  msgs->pos++;
}

static void msgsTranspose(msgs_t* msgs) {
  uint16_t buffer[16] ATTR_ALIGNED(32);
  size_t pos;
  for (pos = 0; pos < msgs->pos / 16; pos++) {
    for (size_t i = 0; i < 16; i++) {
      buffer[i] = ((uint16_t*)msgs->msgs[i])[pos];
    }
    transpose_16_16_uint64(buffer, buffer);
    for (size_t i = 0; i < 16; i++) {
      ((uint16_t*)msgs->msgs[i])[pos] = buffer[i];
    }
  }
  if (msgs->pos % 16) {
    memset(buffer, 0, 16 * sizeof(uint16_t));
    for (size_t i = 0; i < msgs->pos % 16; i++) {
      buffer[i] = ((uint16_t*)msgs->msgs[i])[pos];
    }
    transpose_16_16_uint64(buffer, buffer);
    for (size_t i = 0; i < 16; i++) {
      ((uint16_t*)msgs->msgs[i])[pos] = buffer[i];
    }
  }
}

static uint16_t mpc_AND(uint16_t a, uint16_t b, uint16_t mask_a, uint16_t mask_b,
                        randomTape_t* tapes, msgs_t* msgs, uint8_t* unopened_msg) {
  uint16_t and_helper =
      tapesToWord(tapes); // The special mask value setup during preprocessing for each AND gate
  uint16_t s_shares =
      ((a * 0xFFFF) & mask_b) ^ ((b * 0xFFFF) & mask_a) ^ and_helper; // TODO: nicer extend

  if (msgs->unopened != -1) {
    uint8_t unopenedPartyBit = getBit(unopened_msg, msgs->pos);
    setBit((uint8_t*)&s_shares, msgs->unopened, unopenedPartyBit);
  }

  // Broadcast each share of s
  wordToMsgsNoTranspose(s_shares, msgs);

  return parity64_uint16(s_shares) ^ (a & b);
}

static void mpc_sbox(mzd_local_t* statein, randomTape_t* tapes, msgs_t* msgs,
                     uint8_t* unopenened_msg, const picnic_instance_t* params) {
  uint8_t state[MAX_LOWMC_BLOCK_SIZE];
  uint16_t masks[MAX_LOWMC_BLOCK_SIZE * 8];
  mzd_to_char_array(state, statein, params->output_size);
  for (size_t i = 0; i < params->lowmc.n; i++) {
    masks[i] = tapesToWord(tapes);
  }
  for (size_t i = 0; i < params->lowmc.m * 3; i += 3) {
    uint16_t mask_a = masks[i + 2];
    uint16_t mask_b = masks[i + 1];
    uint16_t mask_c = masks[i];
    uint16_t a      = getBit(state, i + 2);
    uint16_t b      = getBit(state, i + 1);
    uint16_t c      = getBit(state, i);

    uint16_t ab = mpc_AND(a, b, mask_a, mask_b, tapes, msgs, unopenened_msg);
    uint16_t bc = mpc_AND(b, c, mask_b, mask_c, tapes, msgs, unopenened_msg);
    uint16_t ca = mpc_AND(c, a, mask_c, mask_a, tapes, msgs, unopenened_msg);

    uint16_t d = a ^ bc;
    uint16_t e = a ^ b ^ ca;
    uint16_t f = a ^ b ^ c ^ ab;

    setBit(state, i + 2, d);
    setBit(state, i + 1, e);
    setBit(state, i + 0, f);
  }
  mzd_from_char_array(statein, state, params->output_size);
}

#include "lowmc_192_192_4.h"
static void mpc_sbox_192_uint64(mzd_local_t* statein, randomTape_t* tapes, msgs_t* msgs,
                                uint8_t* unopenened_msg, const picnic_instance_t* params) {
  mzd_local_t a[1], b[1], c[1];
  // a
  mzd_and_uint64_192(a, mask_192_192_64_a, statein);
  // b
  mzd_and_uint64_192(b, mask_192_192_64_b, statein);
  // c
  mzd_and_uint64_192(c, mask_192_192_64_c, statein);

  mzd_shift_left_uint64_192(a, a, 2);
  mzd_shift_left_uint64_192(b, b, 1);

  mzd_local_t t0[1], t1[1], t2[1];

  mzd_local_t s_ab[1] = {{{0, 0, 0, 0}}}, s_bc[1] = {{{0, 0, 0, 0}}}, s_ca[1] = {{{0, 0, 0, 0}}};
  for (int i = 0; i < 16; i++) {
    mzd_local_t tmp[1];
    if (i == msgs->unopened) {
      // we are in verify, just grab the broadcast s from the msgs array
      mzd_from_char_array(tmp, &msgs->msgs[i][msgs->pos / 8], 24);
      // a
      mzd_and_uint64_192(t0, mask_192_192_64_a, tmp);
      // b
      mzd_and_uint64_192(t1, mask_192_192_64_b, tmp);
      // c
      mzd_and_uint64_192(t2, mask_192_192_64_c, tmp);
      mzd_shift_left_uint64_192(t0, t0, 2);
      mzd_shift_left_uint64_192(t1, t1, 1);
      mzd_xor_uint64_192(s_ab, t2, s_ab);
      mzd_xor_uint64_192(s_bc, t1, s_bc);
      mzd_xor_uint64_192(s_ca, t0, s_ca);

      continue;
    }
    // make a mzd_local from tape[i] for input_masks
    mzd_local_t mask_a[1], mask_b[1], mask_c[1];
    mzd_from_char_array(tmp, &tapes->tape[i][tapes->pos / 8], 24);
    // a
    mzd_and_uint64_192(mask_a, mask_192_192_64_a, tmp);
    // b
    mzd_and_uint64_192(mask_b, mask_192_192_64_b, tmp);
    // c
    mzd_and_uint64_192(mask_c, mask_192_192_64_c, tmp);
    mzd_shift_left_uint64_192(mask_a, mask_a, 2);
    mzd_shift_left_uint64_192(mask_b, mask_b, 1);

    // make a mzd_local from tape[i] for and_helper
    mzd_local_t and_helper_ab[1], and_helper_bc[1], and_helper_ca[1];
    mzd_from_char_array(tmp, &tapes->tape[i][tapes->pos / 8 + 24], 24);
    // a
    mzd_and_uint64_192(and_helper_ab, mask_192_192_64_c, tmp);
    // b
    mzd_and_uint64_192(and_helper_bc, mask_192_192_64_b, tmp);
    // c
    mzd_and_uint64_192(and_helper_ca, mask_192_192_64_a, tmp);
    mzd_shift_left_uint64_192(and_helper_ca, and_helper_ca, 2);
    mzd_shift_left_uint64_192(and_helper_bc, and_helper_bc, 1);

    // s_ab
    mzd_and_uint64_192(t0, a, mask_b);
    mzd_and_uint64_192(t1, b, mask_a);
    mzd_xor_uint64_192(t0, t0, t1);
    mzd_xor_uint64_192(tmp, t0, and_helper_ab);
    mzd_xor_uint64_192(s_ab, tmp, s_ab);
    // s_bc
    mzd_and_uint64_192(t0, b, mask_c);
    mzd_and_uint64_192(t1, c, mask_b);
    mzd_xor_uint64_192(t0, t0, t1);
    mzd_xor_uint64_192(t0, t0, and_helper_bc);
    mzd_xor_uint64_192(s_bc, t0, s_bc);

    mzd_shift_right_uint64_192(t0, t0, 1);
    mzd_xor_uint64_192(tmp, tmp, t0);
    // s_ca
    mzd_and_uint64_192(t0, c, mask_a);
    mzd_and_uint64_192(t1, a, mask_c);
    mzd_xor_uint64_192(t0, t0, t1);
    mzd_xor_uint64_192(t0, t0, and_helper_ca);
    mzd_xor_uint64_192(s_ca, t0, s_ca);

    mzd_shift_right_uint64_192(t0, t0, 2);
    mzd_xor_uint64_192(tmp, tmp, t0);
    mzd_to_char_array(&msgs->msgs[i][msgs->pos / 8], tmp, 24);
  }
  tapes->pos += 192;
  tapes->pos += 192;
  msgs->pos += 192;
  // b & c
  mzd_and_uint64_192(t0, b, c);
  // c & a
  mzd_and_uint64_192(t1, c, a);
  // a & b
  mzd_and_uint64_192(t2, a, b);
  // (b & c) ^ s_bc
  mzd_xor_uint64_192(t0, t0, s_bc);
  // (c & a) ^ s_ca
  mzd_xor_uint64_192(t1, t1, s_ca);
  // (a & b) ^ s_ab
  mzd_xor_uint64_192(t2, t2, s_ab);

  // (b & c) ^ a
  mzd_xor_uint64_192(t0, t0, a);

  // (c & a) ^ a ^ b
  mzd_xor_uint64_192(t1, t1, a);
  mzd_xor_uint64_192(t1, t1, b);

  // (a & b) ^ a ^ b ^c
  mzd_xor_uint64_192(t2, t2, a);
  mzd_xor_uint64_192(t2, t2, b);
  mzd_xor_uint64_192(t2, t2, c);

  mzd_shift_right_uint64_192(t0, t0, 2);
  mzd_shift_right_uint64_192(t1, t1, 1);

  mzd_xor_uint64_192(t2, t2, t1);
  mzd_xor_uint64_192(statein, t2, t0);
}
static void mpc_sbox_192_s256(mzd_local_t* statein, randomTape_t* tapes, msgs_t* msgs,
                              uint8_t* unopenened_msg, const picnic_instance_t* params) {
  mzd_local_t a[1], b[1], c[1];
  // a
  mzd_and_s256_256(a, mask_192_192_64_a, statein);
  // b
  mzd_and_s256_256(b, mask_192_192_64_b, statein);
  // c
  mzd_and_s256_256(c, mask_192_192_64_c, statein);

  mzd_shift_left_uint64_192(a, a, 2);
  mzd_shift_left_uint64_192(b, b, 1);

  mzd_local_t t0[1], t1[1], t2[1];

  mzd_local_t s_ab[1] = {{{0, 0, 0, 0}}}, s_bc[1] = {{{0, 0, 0, 0}}}, s_ca[1] = {{{0, 0, 0, 0}}};
  for (int i = 0; i < 16; i++) {
    mzd_local_t tmp[1];
    if (i == msgs->unopened) {
      // we are in verify, just grab the broadcast s from the msgs array
      mzd_from_char_array(tmp, &msgs->msgs[i][msgs->pos / 8], 24);
      // a
      mzd_and_s256_256(t0, mask_192_192_64_a, tmp);
      // b
      mzd_and_s256_256(t1, mask_192_192_64_b, tmp);
      // c
      mzd_and_s256_256(t2, mask_192_192_64_c, tmp);
      mzd_shift_left_uint64_192(t0, t0, 2);
      mzd_shift_left_uint64_192(t1, t1, 1);
      mzd_xor_s256_256(s_ab, t2, s_ab);
      mzd_xor_s256_256(s_bc, t1, s_bc);
      mzd_xor_s256_256(s_ca, t0, s_ca);

      continue;
    }
    // make a mzd_local from tape[i] for input_masks
    mzd_local_t mask_a[1], mask_b[1], mask_c[1];
    mzd_from_char_array(tmp, &tapes->tape[i][tapes->pos / 8], 24);
    // a
    mzd_and_s256_256(mask_a, mask_192_192_64_a, tmp);
    // b
    mzd_and_s256_256(mask_b, mask_192_192_64_b, tmp);
    // c
    mzd_and_s256_256(mask_c, mask_192_192_64_c, tmp);
    mzd_shift_left_uint64_192(mask_a, mask_a, 2);
    mzd_shift_left_uint64_192(mask_b, mask_b, 1);

    // make a mzd_local from tape[i] for and_helper
    mzd_local_t and_helper_ab[1], and_helper_bc[1], and_helper_ca[1];
    mzd_from_char_array(tmp, &tapes->tape[i][tapes->pos / 8 + 24], 24);
    // a
    mzd_and_s256_256(and_helper_ab, mask_192_192_64_c, tmp);
    // b
    mzd_and_s256_256(and_helper_bc, mask_192_192_64_b, tmp);
    // c
    mzd_and_s256_256(and_helper_ca, mask_192_192_64_a, tmp);
    mzd_shift_left_uint64_192(and_helper_ca, and_helper_ca, 2);
    mzd_shift_left_uint64_192(and_helper_bc, and_helper_bc, 1);

    // s_ab
    mzd_and_s256_256(t0, a, mask_b);
    mzd_and_s256_256(t1, b, mask_a);
    mzd_xor_s256_256(t0, t0, t1);
    mzd_xor_s256_256(tmp, t0, and_helper_ab);
    mzd_xor_s256_256(s_ab, tmp, s_ab);
    // s_bc
    mzd_and_s256_256(t0, b, mask_c);
    mzd_and_s256_256(t1, c, mask_b);
    mzd_xor_s256_256(t0, t0, t1);
    mzd_xor_s256_256(t0, t0, and_helper_bc);
    mzd_xor_s256_256(s_bc, t0, s_bc);

    mzd_shift_right_uint64_192(t0, t0, 1);
    mzd_xor_s256_256(tmp, tmp, t0);
    // s_ca
    mzd_and_s256_256(t0, c, mask_a);
    mzd_and_s256_256(t1, a, mask_c);
    mzd_xor_s256_256(t0, t0, t1);
    mzd_xor_s256_256(t0, t0, and_helper_ca);
    mzd_xor_s256_256(s_ca, t0, s_ca);

    mzd_shift_right_uint64_192(t0, t0, 2);
    mzd_xor_s256_256(tmp, tmp, t0);
    mzd_to_char_array(&msgs->msgs[i][msgs->pos / 8], tmp, 24);
  }
  tapes->pos += 192;
  tapes->pos += 192;
  msgs->pos += 192;
  // b & c
  mzd_and_s256_256(t0, b, c);
  // c & a
  mzd_and_s256_256(t1, c, a);
  // a & b
  mzd_and_s256_256(t2, a, b);
  // (b & c) ^ s_bc
  mzd_xor_s256_256(t0, t0, s_bc);
  // (c & a) ^ s_ca
  mzd_xor_s256_256(t1, t1, s_ca);
  // (a & b) ^ s_ab
  mzd_xor_s256_256(t2, t2, s_ab);

  // (b & c) ^ a
  mzd_xor_s256_256(t0, t0, a);

  // (c & a) ^ a ^ b
  mzd_xor_s256_256(t1, t1, a);
  mzd_xor_s256_256(t1, t1, b);

  // (a & b) ^ a ^ b ^c
  mzd_xor_s256_256(t2, t2, a);
  mzd_xor_s256_256(t2, t2, b);
  mzd_xor_s256_256(t2, t2, c);

  mzd_shift_right_uint64_192(t0, t0, 2);
  mzd_shift_right_uint64_192(t1, t1, 1);

  mzd_xor_s256_256(t2, t2, t1);
  mzd_xor_s256_256(statein, t2, t0);
}

#if defined(WITH_LOWMC_129_129_4)
#include "lowmc_129_129_4.h"
#endif
#if defined(WITH_LOWMC_192_192_4)
#include "lowmc_192_192_4.h"
#endif
#if defined(WITH_LOWMC_255_255_4)
#include "lowmc_255_255_4.h"
#endif

#if !defined(NO_UINT64_FALLBACK)
/* PICNIC3_L1_FS */
#include "lowmc_129_129_4_fns_uint64.h"
#undef SIM_ONLINE
#define SIM_ONLINE lowmc_simulate_online_uint64_129_43
#include "picnic3_simulate.c.i"

/* PICNIC3_L3_FS */
#include "lowmc_192_192_4_fns_uint64.h"
#undef SIM_ONLINE
#define SIM_ONLINE lowmc_simulate_online_uint64_192_64
#include "picnic3_simulate.c.i"

/* PICNIC3_L5_FS */
#include "lowmc_255_255_4_fns_uint64.h"
#undef SIM_ONLINE
#define SIM_ONLINE lowmc_simulate_online_uint64_255_85
#include "picnic3_simulate.c.i"
#endif

#if defined(WITH_OPT)
#if defined(WITH_SSE2) || defined(WITH_NEON)
#undef FN_ATTR
#if defined(WITH_SSE2)
#define FN_ATTR ATTR_TARGET_SSE2
#endif
/* PICNIC3_L1_FS */
#include "lowmc_129_129_4_fns_s128.h"
#undef SIM_ONLINE
#define SIM_ONLINE lowmc_simulate_online_s128_129_43
#include "picnic3_simulate.c.i"

/* PICNIC3_L3_FS */
#include "lowmc_192_192_4_fns_s128.h"
#undef SIM_ONLINE
#define SIM_ONLINE lowmc_simulate_online_s128_192_64
#include "picnic3_simulate.c.i"

/* PICNIC3_L5_FS */
#include "lowmc_255_255_4_fns_s128.h"
#undef SIM_ONLINE
#define SIM_ONLINE lowmc_simulate_online_s128_255_85
#include "picnic3_simulate.c.i"

#endif // SSE/NEON

#if defined(WITH_AVX2)
#undef FN_ATTR
#define FN_ATTR ATTR_TARGET_AVX2
/* PICNIC3_L1_FS */
#include "lowmc_129_129_4_fns_s256.h"
#undef SIM_ONLINE
#define SIM_ONLINE lowmc_simulate_online_s256_129_43
#include "picnic3_simulate.c.i"

/* PICNIC3_L3_FS */
#include "lowmc_192_192_4_fns_s256.h"
#undef SIM_ONLINE
#define SIM_ONLINE lowmc_simulate_online_s256_192_64
#include "picnic3_simulate.c.i"

/* PICNIC3_L5_FS */
#include "lowmc_255_255_4_fns_s256.h"
#undef SIM_ONLINE
#define SIM_ONLINE lowmc_simulate_online_s256_255_85
#include "picnic3_simulate.c.i"

#endif // AVX2
#endif // WITH_OPT

lowmc_simulate_online_f lowmc_simulate_online_get_implementation(const lowmc_parameters_t* lowmc) {
  assert((lowmc->m == 43 && lowmc->n == 129) || (lowmc->m == 64 && lowmc->n == 192) ||
         (lowmc->m == 85 && lowmc->n == 255));

#if defined(WITH_OPT)
#if defined(WITH_AVX2)
  if (CPU_SUPPORTS_AVX2) {
#if defined(WITH_LOWMC_129_129_4)
    if (lowmc->n == 129 && lowmc->m == 43)
      return lowmc_simulate_online_s256_129_43;
#endif
#if defined(WITH_LOWMC_192_192_4)
    if (lowmc->n == 192 && lowmc->m == 64)
      return lowmc_simulate_online_s256_192_64;
#endif
#if defined(WITH_LOWMC_255_255_4)
    if (lowmc->n == 255 && lowmc->m == 85)
      return lowmc_simulate_online_s256_255_85;
#endif
  }
#endif

#if defined(WITH_SSE2) || defined(WITH_NEON)
  if (CPU_SUPPORTS_SSE2 || CPU_SUPPORTS_NEON) {
#if defined(WITH_LOWMC_129_129_4)
    if (lowmc->n == 129 && lowmc->m == 43)
      return lowmc_simulate_online_s128_129_43;
#endif
#if defined(WITH_LOWMC_192_192_4)
    if (lowmc->n == 192 && lowmc->m == 64)
      return lowmc_simulate_online_s128_192_64;
#endif
#if defined(WITH_LOWMC_255_255_4)
    if (lowmc->n == 255 && lowmc->m == 85)
      return lowmc_simulate_online_s128_255_85;
#endif
  }
#endif
#endif

#if !defined(NO_UINT64_FALLBACK)
#if defined(WITH_LOWMC_129_129_4)
  if (lowmc->n == 129 && lowmc->m == 43)
    return lowmc_simulate_online_uint64_129_43;
#endif
#if defined(WITH_LOWMC_192_192_4)
  if (lowmc->n == 192 && lowmc->m == 64)
    return lowmc_simulate_online_uint64_192_64;
#endif
#if defined(WITH_LOWMC_255_255_4)
  if (lowmc->n == 255 && lowmc->m == 85)
    return lowmc_simulate_online_uint64_255_85;
#endif
#endif

  return NULL;
}
