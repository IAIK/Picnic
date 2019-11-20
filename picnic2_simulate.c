/*! @file picnic2_impl.c
 *  @brief This is the main file of the signature scheme for the Picnic2
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
#include "picnic2_simulate.h"
#include "picnic2_simulate_mul.h"

static void wordToMsgsNoTranspose(uint64_t w, msgs_t* msgs) {
  ((uint64_t*)msgs->msgs[msgs->pos % 64])[msgs->pos / 64] = w;
  msgs->pos++;
}

static void msgsTranspose(msgs_t* msgs) {
  uint64_t buffer[64] ATTR_ALIGNED(32);
  size_t pos;
  for (pos = 0; pos < msgs->pos / 64; pos++) {
    for (size_t i = 0; i < 64; i++) {
      buffer[i] = ((uint64_t*)msgs->msgs[i])[pos];
    }
    transpose_64_64(buffer, buffer);
    for (size_t i = 0; i < 64; i++) {
      ((uint64_t*)msgs->msgs[i])[pos] = buffer[i];
    }
  }
  if (msgs->pos % 64) {
    memset(buffer, 0, 64 * sizeof(uint64_t));
    for (size_t i = 0; i < msgs->pos % 64; i++) {
      buffer[i] = ((uint64_t*)msgs->msgs[i])[pos];
    }
    transpose_64_64(buffer, buffer);
    for (size_t i = 0; i < 64; i++) {
      ((uint64_t*)msgs->msgs[i])[pos] = buffer[i];
    }
  }
}

/* For each word in shares; write player i's share to their stream of msgs */
static void broadcast(shares_t* shares, msgs_t* msgs) {
  for (size_t w = 0; w < shares->numWords; w++) {
    wordToMsgsNoTranspose(shares->shares[w], msgs);
  }
}

/* For an input bit b = 0 or 1, return the word of all b bits, i.e.,
 * extend(1) = 0xFFFFFFFFFFFFFFFF
 * extend(0) = 0x0000000000000000
 * Assumes inputs are always 0 or 1.  If this doesn't hold, add "& 1" to the
 * input.
 */
static inline uint64_t extend(uint64_t bit) {
  return ~(bit - 1);
}

static uint8_t mpc_AND(uint8_t a, uint8_t b, uint64_t mask_a, uint64_t mask_b, randomTape_t* tapes,
                       msgs_t* msgs, uint8_t* unopened_msg) {
  uint64_t and_helper =
      tapesToWord(tapes); // The special mask value setup during preprocessing for each AND gate
  uint64_t s_shares = (extend(a) & mask_b) ^ (extend(b) & mask_a) ^ and_helper;

  if (msgs->unopened >= 0) {
    uint8_t unopenedPartyBit = getBit(unopened_msg, msgs->pos);
    setBit((uint8_t*)&s_shares, msgs->unopened, unopenedPartyBit);
  }

  // Broadcast each share of s
  wordToMsgsNoTranspose(s_shares, msgs);

  return (uint8_t)(parity64_uint64(s_shares) ^ (a & b));
}

static void mpc_sbox(mzd_local_t* statein, shares_t* state_masks, randomTape_t* tapes, msgs_t* msgs,
                     uint8_t* unopenened_msg, const picnic_instance_t* params) {
  uint8_t state[MAX_LOWMC_BLOCK_SIZE];
  mzd_to_char_array(state, statein, params->output_size);
  for (size_t i = 0; i < params->lowmc->m * 3; i += 3) {
    uint8_t a       = getBit((uint8_t*)state, i + 2);
    uint64_t mask_a = state_masks->shares[i + 2];

    uint8_t b       = getBit((uint8_t*)state, i + 1);
    uint64_t mask_b = state_masks->shares[i + 1];

    uint8_t c       = getBit((uint8_t*)state, i);
    uint64_t mask_c = state_masks->shares[i];

    uint8_t ab = mpc_AND(a, b, mask_a, mask_b, tapes, msgs, unopenened_msg);
    uint8_t bc = mpc_AND(b, c, mask_b, mask_c, tapes, msgs, unopenened_msg);
    uint8_t ca = mpc_AND(c, a, mask_c, mask_a, tapes, msgs, unopenened_msg);

    setBit((uint8_t*)state, i + 2, a ^ bc);
    setBit((uint8_t*)state, i + 1, a ^ b ^ ca);
    setBit((uint8_t*)state, i, a ^ b ^ c ^ ab);
  }
  size_t filler = params->input_size * 8 - params->lowmc->n;
  tapes->pos += filler;
  mzd_from_char_array(statein, state, params->output_size);
}

static void mpc_xor_masks(shares_t* out, const shares_t* a, const shares_t* b) {
  assert(out->numWords == a->numWords && a->numWords == b->numWords);

  for (size_t i = 0; i < out->numWords; i++) {
    out->shares[i] = a->shares[i] ^ b->shares[i];
  }
}

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
#include "lowmc_fns_uint64_L1.h"
/* PICNIC2_L1_FS */
#undef SIM_ONLINE
#define SIM_ONLINE lowmc_simulate_online_uint64_126_42
#include "picnic2_simulate.c.i"

/* PICNIC2_L3_FS */
#include "lowmc_fns_uint64_L3.h"
#undef SIM_ONLINE
#define SIM_ONLINE lowmc_simulate_online_uint64_192_64
#include "picnic2_simulate.c.i"

/* PICNIC2_L5_FS */
#include "lowmc_fns_uint64_L5.h"
#undef SIM_ONLINE
#define SIM_ONLINE lowmc_simulate_online_uint64_255_85
#include "picnic2_simulate.c.i"
#endif

#if defined(WITH_OPT)
#if defined(WITH_SSE2) || defined(WITH_NEON)
#undef FN_ATTR
#if defined(WITH_SSE2)
#define FN_ATTR ATTR_TARGET_SSE2
#endif
/* PICNIC2_L1_FS */
#include "lowmc_fns_s128_L1.h"
#undef SIM_ONLINE
#define SIM_ONLINE lowmc_simulate_online_s128_126_42
#include "picnic2_simulate.c.i"

/* PICNIC2_L3_FS */
#include "lowmc_fns_s128_L3.h"
#undef SIM_ONLINE
#define SIM_ONLINE lowmc_simulate_online_s128_192_64
#include "picnic2_simulate.c.i"

/* PICNIC2_L5_FS */
#include "lowmc_fns_s128_L3.h"
#undef SIM_ONLINE
#define SIM_ONLINE lowmc_simulate_online_s128_255_85
#include "picnic2_simulate.c.i"

#endif // SSE/NEON

#if defined(WITH_AVX2)
#undef FN_ATTR
#define FN_ATTR ATTR_TARGET_AVX2
/* PICNIC2_L1_FS */
#include "lowmc_fns_s256_L1.h"
#undef SIM_ONLINE
#define SIM_ONLINE lowmc_simulate_online_s256_126_42
#include "picnic2_simulate.c.i"

/* PICNIC2_L3_FS */
#include "lowmc_fns_s256_L3.h"
#undef SIM_ONLINE
#define SIM_ONLINE lowmc_simulate_online_s256_192_64
#include "picnic2_simulate.c.i"

/* PICNIC2_L5_FS */
#include "lowmc_fns_s256_L5.h"
#undef SIM_ONLINE
#define SIM_ONLINE lowmc_simulate_online_s256_255_85
#include "picnic2_simulate.c.i"

#endif // AVX2
#endif // WITH_OPT

lowmc_simulate_online_f lowmc_simulate_online_get_implementation(const lowmc_t* lowmc) {
  ASSUME(lowmc->m == 42 || lowmc->m == 64 || lowmc->m == 85);
  ASSUME(lowmc->n == 126 || lowmc->n == 192 || lowmc->n == 255);

#if defined(WITH_OPT)
#if defined(WITH_AVX2)
  if (CPU_SUPPORTS_AVX2) {
#if defined(WITH_LOWMC_126_126_4)
    if (lowmc->n == 126 && lowmc->m == 42)
      return lowmc_simulate_online_s256_126_42;
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
#if defined(WITH_LOWMC_126_126_4)
    if (lowmc->n == 126 && lowmc->m == 42)
      return lowmc_simulate_online_s128_126_42;
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
#if defined(WITH_LOWMC_126_126_4)
    if (lowmc->n == 126 && lowmc->m == 42)
      return lowmc_simulate_online_uint64_126_42;
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
