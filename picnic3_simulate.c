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
#include "picnic3_simulate_mul.h"

#define PACKING_FACTOR 4
#define PARTIES_LOG 4

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
//static void broadcast(shares_t* shares, msgs_t* msgs) {
//  for (size_t w = 0; w < shares->numWords; w++) {
//    wordToMsgsNoTranspose(shares->shares[w], msgs);
//  }
//}

static uint64_t mpc_AND(uint64_t a, uint64_t b, uint64_t mask_a, uint64_t mask_b, randomTape_t* tapes,
                       msgs_t* msgs, uint8_t** unopened_msg) {
  uint64_t and_helper =
      tapesToWord(tapes); // The special mask value setup during preprocessing for each AND gate
  uint64_t s_shares = (a & mask_b) ^ (b & mask_a) ^ and_helper;

  if (msgs->unopened != NULL) {
    for (uint32_t k = 0; k < PACKING_FACTOR; k++) {
      uint8_t unopenedPartyBit = getBit(unopened_msg[k], msgs->pos);
      setBit((uint8_t*)&s_shares, (64 / PACKING_FACTOR) * k + msgs->unopened[k], unopenedPartyBit);
    }
  }

  // Broadcast each share of s
  wordToMsgsNoTranspose(s_shares, msgs);
  for (uint32_t k = 0; k < PARTIES_LOG; k++) {
    s_shares ^= s_shares >> (1 << (PARTIES_LOG - 1 - k));
  }

  return (uint64_t)(s_shares ^ (a & b));
}

static void mpc_sbox(mzd_local_t** statein, shares_t* state_masks, randomTape_t* tapes,
                     msgs_t* msgs, uint8_t** unopenened_msg, const picnic_instance_t* params) {
  uint8_t state[PACKING_FACTOR][MAX_LOWMC_BLOCK_SIZE];
  for (uint32_t k = 0; k < PACKING_FACTOR; k++) {
    mzd_to_char_array(state[k], statein[k], params->output_size);
  }
  for (size_t i = 0; i < params->lowmc.m * 3; i += 3) {
    uint64_t a = 0;
    uint64_t b = 0;
    uint64_t c = 0;

    uint64_t mask_a = state_masks->shares[i + 2];
    uint64_t mask_b = state_masks->shares[i + 1];
    uint64_t mask_c = state_masks->shares[i];
    for (uint32_t k = 0; k < PACKING_FACTOR; k++) {
      a <<= (64 / PACKING_FACTOR);
      a |= getBit(state[PACKING_FACTOR - 1 - k], i + 2);
      b <<= (64 / PACKING_FACTOR);
      b |= getBit(state[PACKING_FACTOR - 1 - k], i + 1);
      c <<= (64 / PACKING_FACTOR);
      c |= getBit(state[PACKING_FACTOR - 1 - k], i);
    }
    for (size_t k = 0; k < PARTIES_LOG; k++) {
      a ^= (a << (1 << k));
      b ^= (b << (1 << k));
      c ^= (c << (1 << k));
    }

    uint64_t ab = mpc_AND(a, b, mask_a, mask_b, tapes, msgs, unopenened_msg);
    uint64_t bc = mpc_AND(b, c, mask_b, mask_c, tapes, msgs, unopenened_msg);
    uint64_t ca = mpc_AND(c, a, mask_c, mask_a, tapes, msgs, unopenened_msg);

    uint64_t d = a ^ bc;
    uint64_t e = a ^ b ^ ca;
    uint64_t f = a ^ b ^ c ^ ab;

    for (uint32_t k = 0; k < PACKING_FACTOR; k++) {
      setBit(state[k], i + 2, d & 1);
      d >>= (64 / PACKING_FACTOR);
      setBit(state[k], i + 1, e & 1);
      e >>= (64 / PACKING_FACTOR);
      setBit(state[k], i, f & 1);
      f >>= (64 / PACKING_FACTOR);
    }
  }
  for (uint32_t k = 0; k < PACKING_FACTOR; k++) {
    mzd_from_char_array(statein[k], state[k], params->output_size);
  }
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
