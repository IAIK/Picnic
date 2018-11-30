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

#include "bitstream.h"
#include "compat.h"
#include "io.h"
#include "kdf_shake.h"
#include "lowmc.h"
#include "lowmc_pars.h"
#include "mpc.h"
#include "mpc_lowmc.h"
#include "picnic_impl.h"
#include "randomness.h"

#include <limits.h>
#include <math.h>
#include <stdlib.h>

typedef struct {
  uint8_t* seeds[SC_PROOF];
  uint8_t* commitments[SC_PROOF];
  uint8_t* gs[SC_PROOF];
  uint8_t* input_shares[SC_PROOF];
  uint8_t* communicated_bits[SC_PROOF];
  uint8_t* output_shares[SC_PROOF];
} proof_round_t;

typedef struct {
  uint8_t* challenge;
  uint8_t* salt;
  proof_round_t round[];
} sig_proof_t;

// Prefix values for domain separation
static const uint8_t HASH_PREFIX_0 = 0;
static const uint8_t HASH_PREFIX_1 = 1;
static const uint8_t HASH_PREFIX_2 = 2;
static const uint8_t HASH_PREFIX_4 = 4;
static const uint8_t HASH_PREFIX_5 = 5;

#define LOWMC_UNSPECFIED_ARG UINT32_MAX

/**
 * Collapse challenge from one char per challenge to bit array.
 */
static void collapse_challenge(uint8_t* collapsed, const picnic_instance_t* pp,
                               const uint8_t* challenge) {
  bitstream_t bs;
  bs.buffer   = collapsed;
  bs.position = 0;

  for (unsigned int i = 0; i < pp->num_rounds; ++i) {
    // flip challenge bits according to spec
    bitstream_put_bits(&bs, (challenge[i] >> 1) | ((challenge[i] & 1) << 1), 2);
  }
}

/**
 * Expand challenge from bit array to one char per challenge.
 */
static bool expand_challenge(uint8_t* challenge, const picnic_instance_t* pp,
                             const uint8_t* collapsed) {
  bitstream_t bs;
  bs.cbuffer  = collapsed;
  bs.position = 0;

  for (unsigned int i = 0; i < pp->num_rounds; ++i) {
    uint8_t ch = bitstream_get_bits(&bs, 2);
    if (ch == 3) {
      return false;
    }
    // flip challenge bits according to spec
    challenge[i] = (ch & 1) << 1 | (ch >> 1);
  }

  const size_t remaining_bits = (pp->collapsed_challenge_size << 3) - bs.position;
  if (remaining_bits && bitstream_get_bits(&bs, remaining_bits)) {
    return false;
  }

  return true;
}

#define ALIGNT(s, t) (((s) + sizeof(t) - 1) & ~(sizeof(t) - 1))
#define ALIGNU64T(s) ALIGNT(s, uint64_t)

static sig_proof_t* proof_new(const picnic_instance_t* pp) {
  const size_t digest_size                    = pp->digest_size;
  const size_t seed_size                      = pp->seed_size;
  const size_t num_rounds                     = pp->num_rounds;
  const size_t input_size                     = pp->input_size;
  const size_t output_size                    = pp->output_size;
  const size_t view_size                      = ALIGNU64T(pp->view_size);
  const size_t unruh_with_input_bytes_size    = pp->unruh_with_input_bytes_size;
  const size_t unruh_without_input_bytes_size = pp->unruh_without_input_bytes_size;

  sig_proof_t* prf = calloc(1, sizeof(sig_proof_t) + num_rounds * sizeof(proof_round_t));

  size_t per_round_mem =
      SC_PROOF * (seed_size + digest_size + input_size + output_size + view_size);
  if (pp->transform == TRANSFORM_UR) {
    per_round_mem += (SC_PROOF - 1) * unruh_without_input_bytes_size + unruh_with_input_bytes_size;
  }

  // in memory:
  // - challenge (aligned to uint64_t)
  // - seeds
  // - salt
  // - commitments
  // - input shares
  // - communicated bits (aligned to uint64_t)
  // - output shares
  // - Gs
  //
  // Since seeds size, commitment size, input share size and output share size are all divisible by
  // the alignment of uint64_t, this means, that up to the memory of the Gs, everything is
  // uint64_t-aligned.
  uint8_t* slab  = calloc(1, num_rounds * per_round_mem + ALIGNU64T(num_rounds) + seed_size);
  prf->challenge = slab;
  slab += ALIGNU64T(num_rounds);

  for (uint32_t r = 0; r < num_rounds; ++r) {
    for (uint32_t i = 0; i < SC_PROOF; ++i) {
      prf->round[r].seeds[i] = slab;
      slab += seed_size;
    }
  }
  prf->salt = slab;
  slab += seed_size;

  for (uint32_t r = 0; r < num_rounds; ++r) {
    for (uint32_t i = 0; i < SC_PROOF; ++i) {
      prf->round[r].commitments[i] = slab;
      slab += digest_size;
    }
  }

  for (uint32_t r = 0; r < num_rounds; ++r) {
    for (uint32_t i = 0; i < SC_PROOF; ++i) {
      prf->round[r].input_shares[i] = slab;
      slab += input_size;
    }
  }

  for (uint32_t r = 0; r < num_rounds; ++r) {
    for (uint32_t i = 0; i < SC_PROOF; ++i) {
      prf->round[r].communicated_bits[i] = slab;
      slab += view_size;
    }
  }

  for (uint32_t r = 0; r < num_rounds; ++r) {
    for (uint32_t i = 0; i < SC_PROOF; ++i) {
      prf->round[r].output_shares[i] = slab;
      slab += output_size;
    }
  }

  if (pp->transform == TRANSFORM_UR) {
    for (uint32_t r = 0; r < num_rounds; ++r) {
      for (uint32_t i = 0; i < SC_PROOF - 1; ++i) {
        prf->round[r].gs[i] = slab;
        slab += unruh_without_input_bytes_size;
      }
      prf->round[r].gs[SC_PROOF - 1] = slab;
      slab += unruh_with_input_bytes_size;
    }
  }

  return prf;
}

static sig_proof_t* proof_new_verify(const picnic_instance_t* pp, uint8_t** rslab) {
  const size_t digest_size                 = pp->digest_size;
  const size_t num_rounds                  = pp->num_rounds;
  const size_t input_size                  = pp->input_size;
  const size_t output_size                 = pp->output_size;
  const size_t view_size                   = ALIGNU64T(pp->view_size);
  const size_t seed_size                   = pp->seed_size;
  const size_t unruh_with_input_bytes_size = pp->unruh_with_input_bytes_size;

  sig_proof_t* proof = calloc(1, sizeof(sig_proof_t) + num_rounds * sizeof(proof_round_t));

  size_t per_round_mem = SC_VERIFY * digest_size;
  if (pp->transform == TRANSFORM_UR) {
    // we don't know what we actually need, so allocate more than needed
    per_round_mem += SC_VERIFY * pp->unruh_with_input_bytes_size;
  }
  per_round_mem += SC_VERIFY * input_size + SC_PROOF * output_size + view_size;

  uint8_t* slab    = calloc(1, num_rounds * per_round_mem + ALIGNU64T(num_rounds) + seed_size);
  proof->challenge = slab;
  slab += ALIGNU64T(num_rounds);

  proof->salt = slab;
  slab += seed_size;

  for (uint32_t r = 0; r < num_rounds; ++r) {
    for (uint32_t i = 0; i < SC_VERIFY; ++i) {
      proof->round[r].commitments[i] = slab;
      slab += digest_size;
    }
  }

  for (uint32_t r = 0; r < num_rounds; ++r) {
    proof->round[r].communicated_bits[0] = slab;
    slab += view_size;
  }

  for (uint32_t r = 0; r < num_rounds; ++r) {
    proof->round[r].output_shares[0] = slab;
    slab += output_size;
    proof->round[r].output_shares[1] = slab;
    slab += output_size;
    proof->round[r].output_shares[2] = slab;
    slab += output_size;
  }

  if (pp->transform == TRANSFORM_UR) {
    for (uint32_t r = 0; r < num_rounds; ++r) {
      for (uint32_t i = 0; i < SC_VERIFY; ++i) {
        proof->round[r].gs[i] = slab;
        slab += unruh_with_input_bytes_size;
      }
    }
  }

  *rslab = slab;
  return proof;
}

static void proof_free(sig_proof_t* prf) {
  free(prf->challenge);
  free(prf);
}

static void kdf_shake_update_key_intLE(kdf_shake_t* kdf, uint16_t x) {
  const uint16_t x_le = htole16(x);
  kdf_shake_update_key(kdf, (const uint8_t*)&x_le, sizeof(x_le));
}

static void kdf_init_from_seed(kdf_shake_t* kdf, const uint8_t* seed, const uint8_t* salt, uint16_t round_number,
                               uint16_t player_number, bool include_input_size, const picnic_instance_t* pp) {
  // Hash the seed with H_2.
  kdf_shake_init_prefix(kdf, pp, HASH_PREFIX_2);
  kdf_shake_update_key(kdf, seed, pp->seed_size);
  kdf_shake_finalize_key(kdf);

  uint8_t tmp[MAX_DIGEST_SIZE];
  kdf_shake_get_randomness(kdf, tmp, pp->digest_size);
  kdf_shake_clear(kdf);

  // Initialize KDF with H_2(seed) || salt || round_number || player_number || output_size.
  kdf_shake_init(kdf, pp);
  kdf_shake_update_key(kdf, tmp, pp->digest_size);
  kdf_shake_update_key(kdf, salt, pp->seed_size);
  kdf_shake_update_key_intLE(kdf, round_number);
  kdf_shake_update_key_intLE(kdf, player_number);
  kdf_shake_update_key_intLE(kdf, pp->view_size + (include_input_size ? pp->input_size : 0));
  kdf_shake_finalize_key(kdf);
}

#if defined(WITH_CUSTOM_INSTANCES)
static void mzd_to_bitstream(bitstream_t* bs, const mzd_local_t* v, const size_t size) {
  const uint64_t* d = &CONST_FIRST_ROW(v)[v->width - 1];
  size_t bits       = size;
  for (; bits >= sizeof(uint64_t) * 8; bits -= sizeof(uint64_t) * 8, --d) {
    bitstream_put_bits(bs, *d, sizeof(uint64_t) * 8);
  }
  if (bits) {
    bitstream_put_bits(bs, *d >> (sizeof(uint64_t) * 8 - bits), bits);
  }
}

static void mzd_from_bitstream(bitstream_t* bs, mzd_local_t* v, const size_t size) {
  uint64_t* d = &FIRST_ROW(v)[v->width - 1];
  uint64_t* f = FIRST_ROW(v);

  size_t bits = size;
  for (; bits >= sizeof(uint64_t) * 8; bits -= sizeof(uint64_t) * 8, --d) {
    *d = bitstream_get_bits(bs, sizeof(uint64_t) * 8);
  }
  if (bits) {
    *d = bitstream_get_bits(bs, bits) << (sizeof(uint64_t) * 8 - bits);
    --d;
  }
  for (; d >= f; --d) {
    *d = 0;
  }
}
#endif

static void uint64_to_bitstream(bitstream_t* bs, const uint64_t v) {
  bitstream_put_bits(bs, v >> (64 - 30), 30);
}

static uint64_t uint64_from_bitstream(bitstream_t* bs) {
  return bitstream_get_bits(bs, 30) << (64 - 30);
}

static void compress_view(uint8_t* dst, const picnic_instance_t* pp, const view_t* views,
                          unsigned int idx) {
  const size_t num_views = pp->lowmc->r;

  bitstream_t bs;
  bs.buffer   = dst;
  bs.position = 0;

  const view_t* v = &views[0];
#if defined(WITH_CUSTOM_INSTANCES)
  if (pp->lowmc->m != 10) {
    const size_t view_round_size = pp->view_round_size;
    for (size_t i = 0; i < num_views; ++i, ++v) {
      mzd_to_bitstream(&bs, v->s[idx], view_round_size);
    }
    return;
  }
#endif

  for (size_t i = 0; i < num_views; ++i, ++v) {
    uint64_to_bitstream(&bs, v->t[idx]);
  }
}

static void decompress_view(view_t* views, const picnic_instance_t* pp, const uint8_t* src,
                            unsigned int idx) {
  const size_t num_views = pp->lowmc->r;

  bitstream_t bs;
  bs.buffer   = (uint8_t*)src;
  bs.position = 0;

  view_t* v = &views[0];
#if defined(WITH_CUSTOM_INSTANCES)
  if (pp->lowmc->m != 10) {
    const size_t view_round_size = pp->view_round_size;
    for (size_t i = 0; i < num_views; ++i, ++v) {
      mzd_from_bitstream(&bs, v->s[idx], view_round_size);
    }
    return;
  }
#endif

  for (size_t i = 0; i < num_views; ++i, ++v) {
    v->t[idx] = uint64_from_bitstream(&bs);
  }
}

static void decompress_random_tape(rvec_t* rvec, const picnic_instance_t* pp, const uint8_t* src,
                                   unsigned int idx) {
  const size_t num_views = pp->lowmc->r;

  bitstream_t bs;
  bs.buffer   = (uint8_t*)src;
  bs.position = 0;

  rvec_t* rv = &rvec[0];
#if defined(WITH_CUSTOM_INSTANCES)
  if (pp->lowmc->m != 10) {
    const size_t view_round_size = pp->view_round_size;
    for (size_t i = 0; i < num_views; ++i, ++rv) {
      mzd_from_bitstream(&bs, rv->s[idx], view_round_size);
    }
    return;
  }
#endif

  for (size_t i = 0; i < num_views; ++i, ++rv) {
    rv->t[idx] = uint64_from_bitstream(&bs);
  }
}

static void mzd_share(mzd_local_t* shared_value[SC_PROOF], const mzd_local_t* value) {
  mzd_xor(shared_value[2], shared_value[0], value);
  mzd_xor(shared_value[2], shared_value[1], shared_value[2]);
}

static void mzd_unshare(mzd_local_t* shared_value[SC_PROOF], const mzd_local_t* value) {
  mzd_xor(shared_value[2], shared_value[0], shared_value[1]);
  mzd_xor(shared_value[2], shared_value[2], value);
}

/**
 * Compute commitment to a view.
 */
static void hash_commitment(const picnic_instance_t* pp, proof_round_t* prf_round,
                            unsigned int vidx) {
  const size_t hashlen = pp->digest_size;

  hash_context ctx;
  // hash the seed
  hash_init_prefix(&ctx, pp, HASH_PREFIX_4);
  hash_update(&ctx, prf_round->seeds[vidx], pp->seed_size);
  hash_final(&ctx);
  uint8_t tmp[MAX_DIGEST_SIZE];
  hash_squeeze(&ctx, tmp, hashlen);

  // compute H_0(H_4(seed), view)
  hash_init_prefix(&ctx, pp, HASH_PREFIX_0);
  hash_update(&ctx, tmp, hashlen);
  // hash input share
  hash_update(&ctx, prf_round->input_shares[vidx], pp->input_size);
  // hash communicated bits
  hash_update(&ctx, prf_round->communicated_bits[vidx], pp->view_size);
  // hash output share
  hash_update(&ctx, prf_round->output_shares[vidx], pp->output_size);
  hash_final(&ctx);
  hash_squeeze(&ctx, prf_round->commitments[vidx], hashlen);
}

/**
 * Compute challenge from transform dependent hash - outputs {1,2 or 3}^t
 */
static void H3_compute(const picnic_instance_t* pp, uint8_t* hash, uint8_t* ch) {
  const size_t digest_size      = pp->digest_size;
  const size_t digest_size_bits = digest_size << 3;

  // Pick bits from hash
  uint8_t* eof   = ch + pp->num_rounds;
  size_t bit_idx = 0;
  while (ch < eof) {
    if (bit_idx >= digest_size_bits) {
      hash_context ctx;
      hash_init_prefix(&ctx, pp, HASH_PREFIX_1);
      hash_update(&ctx, hash, digest_size);
      hash_final(&ctx);
      hash_squeeze(&ctx, hash, digest_size);
      bit_idx = 0;
    }

    const uint8_t twobits = (hash[bit_idx >> 3] >> ((6 - (bit_idx & 0x7)))) & 0x3;
    if (twobits != 0x3) {
      *ch++ = twobits;
    }
    bit_idx += 2;
  }
}

/**
 * Hash public key, salt and message
 */
static void H3_public_key_message(hash_context* ctx, const picnic_instance_t* pp,
                                  const uint8_t* salt, const uint8_t* circuit_output,
                                  const uint8_t* circuit_input, const uint8_t* m, size_t m_len) {
  // hash circuit out and input (public key)
  hash_update(ctx, circuit_output, pp->output_size);
  hash_update(ctx, circuit_input, pp->input_size);
  // hash salt
  hash_update(ctx, salt, pp->seed_size);
  // hash message
  hash_update(ctx, m, m_len);
}

/**
 * Re-compute challenge for verification
 */
static void H3_verify(const picnic_instance_t* pp, sig_proof_t* prf, const uint8_t* circuit_output,
                      const uint8_t* circuit_input, const uint8_t* m, size_t m_len, uint8_t* ch) {
  const size_t digest_size = pp->digest_size;
  const size_t num_rounds  = pp->num_rounds;
  const size_t output_size = pp->output_size;

  hash_context ctx;
  hash_init_prefix(&ctx, pp, HASH_PREFIX_1);

  // hash output shares
  proof_round_t* round = prf->round;
  for (size_t i = 0; i < num_rounds; ++i, ++round) {
    switch (prf->challenge[i]) {
    case 0: {
      hash_update(&ctx, round->output_shares[0], output_size);
      hash_update(&ctx, round->output_shares[1], output_size);
      hash_update(&ctx, round->output_shares[2], output_size);
      break;
    }
    case 1: {
      hash_update(&ctx, round->output_shares[2], output_size);
      hash_update(&ctx, round->output_shares[0], output_size);
      hash_update(&ctx, round->output_shares[1], output_size);
      break;
    }
    default: {
      hash_update(&ctx, round->output_shares[1], output_size);
      hash_update(&ctx, round->output_shares[2], output_size);
      hash_update(&ctx, round->output_shares[0], output_size);
      break;
    }
    }
  }

  // hash commitments
  round = prf->round;
  for (size_t i = 0; i < num_rounds; ++i, ++round) {
    switch (prf->challenge[i]) {
    case 0: {
      hash_update(&ctx, round->commitments[0], digest_size);
      hash_update(&ctx, round->commitments[1], digest_size);
      hash_update(&ctx, round->commitments[2], digest_size);
      break;
    }
    case 1: {
      hash_update(&ctx, round->commitments[2], digest_size);
      hash_update(&ctx, round->commitments[0], digest_size);
      hash_update(&ctx, round->commitments[1], digest_size);
      break;
    }
    default: {
      hash_update(&ctx, round->commitments[1], digest_size);
      hash_update(&ctx, round->commitments[2], digest_size);
      hash_update(&ctx, round->commitments[0], digest_size);
      break;
    }
    }
  }

  if (pp->transform == TRANSFORM_UR) {
    const size_t without_input_bytes_size = pp->unruh_without_input_bytes_size;
    const size_t with_input_bytes_size    = pp->unruh_with_input_bytes_size;

    // hash commitments
    round = prf->round;
    for (size_t i = 0; i < num_rounds; ++i, ++round) {
      switch (prf->challenge[i]) {
      case 0: {
        hash_update(&ctx, round->gs[0], without_input_bytes_size);
        hash_update(&ctx, round->gs[1], without_input_bytes_size);
        hash_update(&ctx, round->gs[2], with_input_bytes_size);
        break;
      }
      case 1: {
        hash_update(&ctx, round->gs[2], without_input_bytes_size);
        hash_update(&ctx, round->gs[0], without_input_bytes_size);
        hash_update(&ctx, round->gs[1], with_input_bytes_size);
        break;
      }
      default: {
        hash_update(&ctx, round->gs[1], without_input_bytes_size);
        hash_update(&ctx, round->gs[2], without_input_bytes_size);
        hash_update(&ctx, round->gs[0], with_input_bytes_size);
        break;
      }
      }
    }
  }
  // hash public key, salt, and message
  H3_public_key_message(&ctx, pp, prf->salt, circuit_output, circuit_input, m, m_len);
  hash_final(&ctx);

  uint8_t hash[MAX_DIGEST_SIZE];
  hash_squeeze(&ctx, hash, digest_size);
  H3_compute(pp, hash, ch);
}

/**
 * Compute challenge
 */
static void H3(const picnic_instance_t* pp, sig_proof_t* prf, const uint8_t* circuit_output,
               const uint8_t* circuit_input, const uint8_t* m, size_t m_len) {
  const size_t num_rounds = pp->num_rounds;

  hash_context ctx;
  hash_init_prefix(&ctx, pp, HASH_PREFIX_1);

  // hash output shares
  hash_update(&ctx, prf->round[0].output_shares[0], pp->output_size * num_rounds * SC_PROOF);
  // hash all commitments C
  hash_update(&ctx, prf->round[0].commitments[0], pp->digest_size * num_rounds * SC_PROOF);
  if (pp->transform == TRANSFORM_UR) {
    // hash all commitments G
    hash_update(&ctx, prf->round[0].gs[0],
                num_rounds * ((SC_PROOF - 1) * pp->unruh_without_input_bytes_size +
                              pp->unruh_with_input_bytes_size));
  }
  // hash public key, salt, and message
  H3_public_key_message(&ctx, pp, prf->salt, circuit_output, circuit_input, m, m_len);
  hash_final(&ctx);

  uint8_t hash[MAX_DIGEST_SIZE];
  hash_squeeze(&ctx, hash, pp->digest_size);
  H3_compute(pp, hash, prf->challenge);
}

/*
 * G permutation for Unruh transform
 */
static void unruh_G(const picnic_instance_t* pp, proof_round_t* prf_round, unsigned int vidx,
                    bool include_is) {
  hash_context ctx;

  const size_t outputlen =
      include_is ? pp->unruh_with_input_bytes_size : pp->unruh_without_input_bytes_size;
  const uint16_t size_le   = htole16(outputlen);
  const size_t digest_size = pp->digest_size;
  const size_t seedlen     = pp->seed_size;

  // Hash the seed with H_5, store digest in output
  hash_init_prefix(&ctx, pp, HASH_PREFIX_5);
  hash_update(&ctx, prf_round->seeds[vidx], seedlen);
  hash_final(&ctx);

  uint8_t tmp[MAX_DIGEST_SIZE];
  hash_squeeze(&ctx, tmp, digest_size);

  // Hash H_5(seed), the view, and the length
  hash_init(&ctx, pp);
  hash_update(&ctx, tmp, digest_size);
  if (include_is) {
    hash_update(&ctx, prf_round->input_shares[vidx], pp->input_size);
  }
  hash_update(&ctx, prf_round->communicated_bits[vidx], pp->view_size);
  hash_update(&ctx, (const uint8_t*)&size_le, sizeof(uint16_t));
  hash_final(&ctx);
  hash_squeeze(&ctx, prf_round->gs[vidx], outputlen);
}

// serilization helper functions
static int sig_proof_to_char_array(const picnic_instance_t* pp, const sig_proof_t* prf,
                                   uint8_t* result, size_t* siglen) {
  const uint32_t num_rounds                   = pp->num_rounds;
  const uint32_t seed_size                    = pp->seed_size;
  const uint32_t challenge_size               = pp->collapsed_challenge_size;
  const uint32_t digest_size                  = pp->digest_size;
  const transform_t transform                 = pp->transform;
  const size_t view_size                      = pp->view_size;
  const size_t input_size                     = pp->input_size;
  const size_t unruh_with_input_bytes_size    = pp->unruh_with_input_bytes_size;
  const size_t unruh_without_input_bytes_size = pp->unruh_without_input_bytes_size;

  uint8_t* tmp = result;

  // write challenge
  collapse_challenge(tmp, pp, prf->challenge);
  tmp += challenge_size;

  // write salt
  memcpy(tmp, prf->salt, pp->seed_size);
  tmp += seed_size;

  const proof_round_t* round = prf->round;
  for (unsigned i = 0; i < num_rounds; ++i, ++round) {
    const unsigned int a = prf->challenge[i];
    const unsigned int b = (a + 1) % 3;
    const unsigned int c = (a + 2) % 3;

    // write commitment
    memcpy(tmp, round->commitments[c], digest_size);
    tmp += digest_size;

    // write unruh G
    if (transform == TRANSFORM_UR) {
      const uint32_t unruh_g_size =
          a ? unruh_without_input_bytes_size : unruh_with_input_bytes_size;
      memcpy(tmp, round->gs[c], unruh_g_size);
      tmp += unruh_g_size;
    }

    // write views
    memcpy(tmp, round->communicated_bits[b], view_size);
    tmp += view_size;

    // write seeds
    memcpy(tmp, round->seeds[a], seed_size);
    tmp += seed_size;
    memcpy(tmp, round->seeds[b], seed_size);
    tmp += seed_size;

    if (a) {
      // write input share
      memcpy(tmp, round->input_shares[SC_PROOF - 1], input_size);
      tmp += input_size;
    }
  }

  *siglen = tmp - result;
  return 0;
}

static sig_proof_t* sig_proof_from_char_array(const picnic_instance_t* pp, const uint8_t* data,
                                              size_t len) {
  const size_t digest_size              = pp->digest_size;
  const size_t seed_size                = pp->seed_size;
  const size_t num_rounds               = pp->num_rounds;
  const size_t challenge_size           = pp->collapsed_challenge_size;
  const transform_t transform           = pp->transform;
  const size_t input_size               = pp->input_size;
  const size_t view_size                = pp->view_size;
  const size_t without_input_bytes_size = pp->unruh_without_input_bytes_size;
  const size_t with_input_bytes_size    = pp->unruh_with_input_bytes_size;

  uint8_t* slab      = NULL;
  sig_proof_t* proof = proof_new_verify(pp, &slab);
  if (!proof) {
    return NULL;
  }

  size_t remaining_len = len;
  const uint8_t* tmp   = data;

  // read and process challenge
  if (sub_overflow_size_t(remaining_len, challenge_size, &remaining_len)) {
    goto err;
  }
  if (!expand_challenge(proof->challenge, pp, tmp)) {
    goto err;
  }
  tmp += challenge_size;

  // read salt
  if(sub_overflow_size_t(remaining_len, seed_size, &remaining_len)) {
      goto err;
  }
  memcpy(proof->salt, tmp, seed_size);
  tmp += seed_size;

  const size_t base_size = digest_size + view_size + 2 * seed_size;
  proof_round_t* round   = proof->round;
  for (unsigned int i = 0; i < num_rounds; ++i, ++round) {
    const unsigned char ch      = proof->challenge[i];
    const size_t unruh_g_len    = ch ? without_input_bytes_size : with_input_bytes_size;
    const size_t requested_size = base_size + unruh_g_len + (ch ? input_size : 0);
    if (sub_overflow_size_t(remaining_len, requested_size, &remaining_len)) {
      goto err;
    }

    // read commitments
    round->commitments[2] = (uint8_t*)tmp;
    tmp += digest_size;

    // read unruh G
    if (transform == TRANSFORM_UR) {
      round->gs[2] = (uint8_t*)tmp;
      tmp += unruh_g_len;
    }

    // read view
    round->communicated_bits[1] = (uint8_t*)tmp;
    tmp += view_size;

    // read seeds
    round->seeds[0] = (uint8_t*)tmp;
    tmp += seed_size;
    round->seeds[1] = (uint8_t*)tmp;
    tmp += seed_size;

    // read input shares
    switch (ch) {
    case 0:
      round->input_shares[0] = slab;
      slab += input_size;
      round->input_shares[1] = slab;
      slab += input_size;
      break;
    case 1:
      round->input_shares[0] = slab;
      slab += input_size;
      round->input_shares[1] = (uint8_t*)tmp;
      tmp += input_size;
      break;
    default:
      round->input_shares[0] = (uint8_t*)tmp;
      tmp += input_size;
      round->input_shares[1] = slab;
      slab += input_size;
    }
  }

  if (remaining_len) {
    goto err;
  }

  return proof;

err:
  proof_free(proof);
  return NULL;
}

static void generate_seeds(const picnic_instance_t* pp, const uint8_t* private_key,
                           const uint8_t* plaintext, const uint8_t* public_key, const uint8_t* m,
                           size_t m_len, uint8_t* seeds, uint8_t* salt) {
  const lowmc_t* lowmc     = pp->lowmc;
  const size_t seed_size   = pp->seed_size;
  const size_t num_rounds  = pp->num_rounds;
  const size_t input_size  = pp->input_size;
  const size_t output_size = pp->output_size;
  const size_t lowmc_n     = lowmc->n;

  kdf_shake_t ctx;
  kdf_shake_init(&ctx, pp);
  // sk || m || C || p
  kdf_shake_update_key(&ctx, private_key, input_size);
  kdf_shake_update_key(&ctx, m, m_len);
  kdf_shake_update_key(&ctx, public_key, output_size);
  kdf_shake_update_key(&ctx, plaintext, output_size);
  // N as 16 bit LE integer
  const uint16_t size_le = htole16(lowmc_n);
  kdf_shake_update_key(&ctx, (const uint8_t*)&size_le, sizeof(size_le));
#if defined(WITH_EXTRA_RANDOMNESS)
  // Add extra randomn bytes for fault attack mitigation
  unsigned char buffer[2 * MAX_DIGEST_SIZE];
  rand_bytes(buffer, 2 * seed_size);
  kdf_shake_update_key(&ctx, buffer, 2 * seed_size);
#endif
  kdf_shake_finalize_key(&ctx);

  // Generate seeds and salt
  kdf_shake_get_randomness(&ctx, seeds, seed_size * num_rounds * SC_PROOF);
  kdf_shake_get_randomness(&ctx, salt, seed_size);
  kdf_shake_clear(&ctx);
}

static int sign_impl(const picnic_instance_t* pp, const uint8_t* private_key,
                     const lowmc_key_t* lowmc_key, const uint8_t* plaintext, const mzd_local_t* p,
                     const uint8_t* public_key, const uint8_t* m, size_t m_len, uint8_t* sig,
                     size_t* siglen) {
  const lowmc_t* lowmc                                = pp->lowmc;
  const zkbpp_lowmc_implementation_f lowmc_impl       = pp->zkbpp_lowmc_impl;
  const lowmc_store_implementation_f lowmc_store_impl = pp->lowmc_store_impl;
  const size_t num_rounds                             = pp->num_rounds;
  const transform_t transform                         = pp->transform;
  const size_t input_size                             = pp->input_size;
  const size_t output_size                            = pp->output_size;
  const size_t view_count                             = lowmc->r;
  const size_t lowmc_k                                = lowmc->k;
  const size_t lowmc_n                                = lowmc->n;
  const size_t lowmc_r                                = lowmc->r;
  const size_t view_size                              = pp->view_size;

  // Perform LowMC evaluation and record state before AND gates
  recorded_state_t recorded_state;
  recorded_state.state = calloc(lowmc_r + 1, sizeof(mzd_local_t*));
  mzd_local_init_multiple_ex(recorded_state.state, lowmc_r + 1, 1, lowmc_n, false);
  lowmc_store_impl(lowmc, lowmc_key, p, &recorded_state);

  sig_proof_t* prf = proof_new(pp);
  view_t* views    = calloc(sizeof(view_t), view_count);
#if defined(WITH_CUSTOM_INSTANCES)
  if (lowmc->m != 10) {
    for (size_t i = 0; i < view_count; ++i) {
      mzd_local_init_multiple_ex(views[i].s, SC_PROOF, 1, lowmc_n, false);
    }
  }
#endif

  in_out_shares_t in_out_shares[2];
  mzd_local_init_multiple_ex(in_out_shares[0].s, SC_PROOF, 1, lowmc_k, false);
  mzd_local_init_multiple_ex(in_out_shares[1].s, SC_PROOF, 1, lowmc_n, false);

  // Generate seeds
  generate_seeds(pp, private_key, plaintext, public_key, m, m_len, prf->round[0].seeds[0], prf->salt);

  mzd_local_t* shared_key[SC_PROOF];
  mzd_local_init_multiple(shared_key, SC_PROOF, 1, lowmc_k);

  rvec_t* rvec = calloc(sizeof(rvec_t), lowmc_r); // random tapes for AND-gates
#if defined(WITH_CUSTOM_INSTANCES)
  if (lowmc->m != 10) {
    for (unsigned int i = 0; i < lowmc_r; ++i) {
      mzd_local_init_multiple_ex(rvec[i].s, SC_PROOF, 1, lowmc_n, false);
    }
  }
#endif

  uint8_t* tape_bytes  = malloc(view_size);
  proof_round_t* round = prf->round;
  for (unsigned int i = 0; i < num_rounds; ++i, ++round) {
    kdf_shake_t kdfs[SC_PROOF];
    for (unsigned int j = 0; j < SC_PROOF; ++j) {
      bool include_input_size = (j != SC_PROOF - 1);
      kdf_init_from_seed(&kdfs[j], round->seeds[j], prf->salt, i, j, include_input_size, pp);
    }

    // compute sharing
    for (unsigned int j = 0; j < SC_PROOF - 1; ++j) {
      kdf_shake_get_randomness(&kdfs[j], round->input_shares[j], input_size);
      mzd_from_char_array(shared_key[j], round->input_shares[j], input_size);
    }
    mzd_share(shared_key, lowmc_key);
    mzd_to_char_array(round->input_shares[SC_PROOF - 1], shared_key[SC_PROOF - 1], input_size);

    // compute random tapes
    for (unsigned int j = 0; j < SC_PROOF; ++j) {
      kdf_shake_get_randomness(&kdfs[j], tape_bytes, view_size);
      decompress_random_tape(rvec, pp, tape_bytes, j);
    }

    for (unsigned int j = 0; j < SC_PROOF; ++j) {
      kdf_shake_clear(&kdfs[j]);
    }

    // perform ZKB++ LowMC evaluation
    lowmc_impl(lowmc, shared_key, p, views, in_out_shares, rvec, &recorded_state);

    // commitments
    for (unsigned int j = 0; j < SC_PROOF; ++j) {
      mzd_to_char_array(round->output_shares[j], in_out_shares[1].s[j], output_size);
      compress_view(round->communicated_bits[j], pp, views, j);
      hash_commitment(pp, round, j);
    }

    // unruh G
    if (transform == TRANSFORM_UR) {
      for (unsigned int j = 0; j < SC_PROOF; ++j) {
        unruh_G(pp, round, j, j == SC_PROOF - 1);
      }
    }
  }
  H3(pp, prf, public_key, plaintext, m, m_len);

  const int ret = sig_proof_to_char_array(pp, prf, sig, siglen);

  // clean up
  free(tape_bytes);
#if defined(WITH_CUSTOM_INSTANCES)
  if (lowmc->m != 10) {
    for (unsigned int n = 0; n < view_count; ++n) {
      mzd_local_free_multiple(rvec[n].s);
    }
    for (unsigned int n = 0; n < view_count; ++n) {
      mzd_local_free_multiple(views[n].s);
    }
  }
#endif
  free(rvec);
  free(views);
  mzd_local_free_multiple(shared_key);
  mzd_local_free_multiple(in_out_shares[1].s);
  mzd_local_free_multiple(in_out_shares[0].s);
  proof_free(prf);

  mzd_local_free_multiple(recorded_state.state);
  free(recorded_state.state);

  return ret;
}

static int verify_impl(const picnic_instance_t* pp, const uint8_t* plaintext, mzd_local_t const* p,
                       const uint8_t* ciphertext, mzd_local_t const* c, const uint8_t* m,
                       size_t m_len, const uint8_t* sig, size_t siglen) {
  const size_t num_rounds                               = pp->num_rounds;
  const lowmc_t* lowmc                                  = pp->lowmc;
  const transform_t transform                           = pp->transform;
  zkbpp_lowmc_verify_implementation_f lowmc_verify_impl = pp->zkbpp_lowmc_verify_impl;
  const size_t input_size                               = pp->input_size;
  const size_t output_size                              = pp->output_size;
  const size_t view_count                               = lowmc->r;
  const size_t lowmc_k                                  = lowmc->k;
  const size_t lowmc_n                                  = lowmc->n;
  const size_t lowmc_r                                  = lowmc->r;
  const size_t view_size                                = pp->view_size;

  sig_proof_t* prf = sig_proof_from_char_array(pp, sig, siglen);
  if (!prf) {
    return -1;
  }

  in_out_shares_t in_out_shares[2];
  mzd_local_init_multiple_ex(in_out_shares[0].s, SC_VERIFY, 1, lowmc_k, false);
  mzd_local_init_multiple_ex(in_out_shares[1].s, SC_PROOF, 1, lowmc_n, false);
  view_t* views = calloc(sizeof(view_t), view_count);
#if defined(WITH_CUSTOM_INSTANCES)
  if (lowmc->m != 10) {
    for (size_t i = 0; i < view_count; ++i) {
      mzd_local_init_multiple_ex(views[i].s, SC_VERIFY, 1, lowmc_n, false);
    }
  }
#endif

  rvec_t* rvec = calloc(sizeof(rvec_t), lowmc_r); // random tapes for and-gates
#if defined(WITH_CUSTOM_INSTANCES)
  if (lowmc->m != 10) {
    for (unsigned int i = 0; i < lowmc_r; ++i) {
      mzd_local_init_multiple_ex(rvec[i].s, SC_VERIFY, 1, lowmc_n, false);
    }
  }
#endif
  uint8_t* tape_bytes = malloc(view_size);

  proof_round_t* round = prf->round;
  for (unsigned int i = 0; i < num_rounds; ++i, ++round) {
    const unsigned int a_i = prf->challenge[i];
    const unsigned int b_i = (a_i + 1) % 3;
    const unsigned int c_i = (a_i + 2) % 3;

    kdf_shake_t kdfs[SC_VERIFY];
    for (unsigned int j = 0; j < SC_VERIFY; ++j) {
      bool include_input_size = (j == 0 && b_i) || (j == 1 && c_i);
      unsigned int player_number = (j == 0)? a_i : b_i;
      kdf_init_from_seed(&kdfs[j], round->seeds[j], prf->salt, i, player_number, include_input_size, pp);
    }

    // compute input shares if necessary
    if (b_i) {
      kdf_shake_get_randomness(&kdfs[0], round->input_shares[0], input_size);
    }
    if (c_i) {
      kdf_shake_get_randomness(&kdfs[1], round->input_shares[1], input_size);
    }

    mzd_from_char_array(in_out_shares[0].s[0], round->input_shares[0], input_size);
    mzd_from_char_array(in_out_shares[0].s[1], round->input_shares[1], input_size);

    // compute random tapes
    for (unsigned int j = 0; j < SC_VERIFY; ++j) {
      kdf_shake_get_randomness(&kdfs[j], tape_bytes, view_size);
      decompress_random_tape(rvec, pp, tape_bytes, j);
    }

    for (unsigned int j = 0; j < SC_VERIFY; ++j) {
      kdf_shake_clear(&kdfs[j]);
    }

    decompress_view(views, pp, round->communicated_bits[1], 1);
    // perform ZKB++ LowMC evaluation
    lowmc_verify_impl(lowmc, p, views, in_out_shares, rvec, a_i);
    compress_view(round->communicated_bits[0], pp, views, 0);

    mzd_unshare(in_out_shares[1].s, c);
    // recompute commitments
    for (unsigned int j = 0; j < SC_VERIFY; ++j) {
      mzd_to_char_array(round->output_shares[j], in_out_shares[1].s[j], output_size);
      hash_commitment(pp, round, j);
    }
    mzd_to_char_array(round->output_shares[SC_VERIFY], in_out_shares[1].s[SC_VERIFY], output_size);

    if (transform == TRANSFORM_UR) {
      // apply Unruh G permutation
      for (unsigned int j = 0; j < SC_VERIFY; ++j) {
        unruh_G(pp, round, j, (a_i == 1 && j == 1) || (a_i == 2 && j == 0));
      }
    }
  }

  unsigned char challenge[MAX_NUM_ROUNDS] = {0};
  H3_verify(pp, prf, ciphertext, plaintext, m, m_len, challenge);
  const int success_status = memcmp(challenge, prf->challenge, pp->num_rounds);

  // clean up
  free(tape_bytes);
#if defined(WITH_CUSTOM_INSTANCES)
  if (lowmc->m != 10) {
    for (unsigned int n = 0; n < view_count; ++n) {
      mzd_local_free_multiple(rvec[n].s);
    }
    for (unsigned int n = 0; n < view_count; ++n) {
      mzd_local_free_multiple(views[n].s);
    }
  }
#endif
  free(rvec);
  free(views);
  mzd_local_free_multiple(in_out_shares[1].s);
  mzd_local_free_multiple(in_out_shares[0].s);

  proof_free(prf);

  return success_status;
}

int impl_sign(const picnic_instance_t* pp, const uint8_t* plaintext, const uint8_t* private_key,
              const uint8_t* public_key, const uint8_t* msg, size_t msglen, uint8_t* sig,
              size_t* siglen) {
  mzd_local_t* m_plaintext  = mzd_local_init_ex(1, pp->lowmc->n, false);
  mzd_local_t* m_privatekey = mzd_local_init_ex(1, pp->lowmc->k, false);

  mzd_from_char_array(m_plaintext, plaintext, pp->output_size);
  mzd_from_char_array(m_privatekey, private_key, pp->input_size);

  const int result = sign_impl(pp, private_key, m_privatekey, plaintext, m_plaintext, public_key,
                               msg, msglen, sig, siglen);

  mzd_local_free(m_privatekey);
  mzd_local_free(m_plaintext);

  return result;
}

int impl_verify(const picnic_instance_t* pp, const uint8_t* plaintext, const uint8_t* public_key,
                const uint8_t* msg, size_t msglen, const uint8_t* sig, size_t siglen) {
  mzd_local_t* m_plaintext = mzd_local_init_ex(1, pp->lowmc->n, false);
  mzd_local_t* m_publickey = mzd_local_init_ex(1, pp->lowmc->n, false);

  mzd_from_char_array(m_plaintext, plaintext, pp->output_size);
  mzd_from_char_array(m_publickey, public_key, pp->output_size);

  const int result =
      verify_impl(pp, plaintext, m_plaintext, public_key, m_publickey, msg, msglen, sig, siglen);

  mzd_local_free(m_publickey);
  mzd_local_free(m_plaintext);

  return result;
}

#if defined(PICNIC_STATIC)
void visualize_signature(FILE* out, const picnic_instance_t* pp, const uint8_t* msg, size_t msglen,
                         const uint8_t* sig, size_t siglen) {
  const size_t digest_size    = pp->digest_size;
  const size_t seed_size      = pp->seed_size;
  const size_t num_rounds     = pp->num_rounds;
  const size_t challenge_size = pp->collapsed_challenge_size;
  const transform_t transform = pp->transform;
  const size_t input_size     = pp->input_size;
  const size_t view_size      = pp->view_size;

  sig_proof_t* proof = sig_proof_from_char_array(pp, sig, siglen);

  fprintf(out, "message: ");
  print_hex(out, msg, msglen);
  fprintf(out, "\nsignature: ");
  print_hex(out, sig, siglen);
  fprintf(out, "\n\n");

  fprintf(out, "challenge: ");
  print_hex(out, sig, challenge_size);
  fprintf(out, "\n\n");

  fprintf(out, "salt: ");
  print_hex(out, proof->salt, seed_size);
  fprintf(out, "\n\n");


  proof_round_t* round = proof->round;
  for (unsigned int i = 0; i < num_rounds; ++i, ++round) {
    fprintf(out, "Iteration t: %d\n", i);

    // print challenge
    const unsigned char ch = proof->challenge[i];
    fprintf(out, "e_%d: %u\n", i, (unsigned int)ch);

    // print commitment
    fprintf(out, "b_%d: ", i);
    print_hex(out, round->commitments[2], digest_size);
    fprintf(out, "\n");

    // print unruh G
    if (transform == TRANSFORM_UR) {
      const size_t unruh_g_len =
          ch ? pp->unruh_without_input_bytes_size : pp->unruh_with_input_bytes_size;

      fprintf(out, "G_%d: ", i);
      print_hex(out, round->gs[2], unruh_g_len);
      fprintf(out, "\n");
    }

    // print view
    fprintf(out, "transcript: ");
    print_hex(out, round->communicated_bits[1], view_size);
    fprintf(out, "\n");

    // print seeds
    fprintf(out, "seed1: ");
    print_hex(out, round->seeds[0], seed_size);
    fprintf(out, "\nseed2: ");
    print_hex(out, round->seeds[1], seed_size);
    fprintf(out, "\n");

    // print input shares
    if (ch == 1) {
      fprintf(out, "inputShare: ");
      print_hex(out, round->input_shares[1], input_size);
      fprintf(out, "\n");
    } else if (ch == 2) {
      fprintf(out, "inputShare: ");
      print_hex(out, round->input_shares[0], input_size);
      fprintf(out, "\n");
    }
    fprintf(out, "\n");
  }

  proof_free(proof);
}
#endif

// instance handling

// L1, L3, and L5 LowMC instances
#if defined(WITH_LOWMC_128_128_20)
#include "lowmc_128_128_20.h"
#define LOWMC_L1_OR_NULL &lowmc_128_128_20
#else
#define LOWMC_L1_OR_NULL NULL
#endif
#if defined(WITH_LOWMC_192_192_30)
#include "lowmc_192_192_30.h"
#define LOWMC_L3_OR_NULL &lowmc_192_192_30
#else
#define LOWMC_L3_OR_NULL NULL
#endif
#if defined(WITH_LOWMC_256_256_38)
#include "lowmc_256_256_38.h"
#define LOWMC_L5_OR_NULL &lowmc_256_256_38
#else
#define LOWMC_L5_OR_NULL NULL
#endif

#if defined(MUL_M4RI)
static lowmc_t* const lowmc_instances[3] = {
#else
static const lowmc_t* const lowmc_instances[3] = {
#endif
    LOWMC_L1_OR_NULL, LOWMC_L3_OR_NULL, LOWMC_L5_OR_NULL};
static bool lowmc_instances_initialized[3];

static picnic_instance_t instances[PARAMETER_SET_MAX_INDEX] = {
    {0},
    {LOWMC_L1_OR_NULL, NULL, NULL, NULL, NULL, 32, 16, 219, 16, 16, 75, 30, 55, 0, 0,
     PICNIC_SIGNATURE_SIZE_Picnic_L1_FS, Picnic_L1_FS, TRANSFORM_FS},
    {LOWMC_L1_OR_NULL, NULL, NULL, NULL, NULL, 32, 16, 219, 16, 16, 75, 30, 55, 91, 107,
     PICNIC_SIGNATURE_SIZE_Picnic_L1_UR, Picnic_L1_UR, TRANSFORM_UR},
    {LOWMC_L3_OR_NULL, NULL, NULL, NULL, NULL, 48, 24, 329, 24, 24, 113, 30, 83, 0, 0,
     PICNIC_SIGNATURE_SIZE_Picnic_L3_FS, Picnic_L3_FS, TRANSFORM_FS},
    {LOWMC_L3_OR_NULL, NULL, NULL, NULL, NULL, 48, 24, 329, 24, 24, 113, 30, 83, 137, 161,
     PICNIC_SIGNATURE_SIZE_Picnic_L3_UR, Picnic_L3_UR, TRANSFORM_UR},
    {LOWMC_L5_OR_NULL, NULL, NULL, NULL, NULL, 64, 32, 438, 32, 32, 143, 30, 110, 0, 0,
     PICNIC_SIGNATURE_SIZE_Picnic_L5_FS, Picnic_L5_FS, TRANSFORM_FS},
    {LOWMC_L5_OR_NULL, NULL, NULL, NULL, NULL, 64, 32, 438, 32, 32, 143, 30, 110, 175, 207,
     PICNIC_SIGNATURE_SIZE_Picnic_L5_UR, Picnic_L5_UR, TRANSFORM_UR}};
static bool instance_initialized[PARAMETER_SET_MAX_INDEX];

static const lowmc_t* lowmc_get_instance(unsigned int idx) {
  if (!lowmc_instances_initialized[idx]) {
#if defined(MUL_M4RI)
    if (lowmc_init(lowmc_instances[idx])) {
      lowmc_instances_initialized[idx] = true;
      return lowmc_instances[idx];
    }
#else
    lowmc_instances_initialized[idx] = true;
    return lowmc_instances[idx];
#endif
  } else {
    return lowmc_instances[idx];
  }

  return NULL;
}

static void clear_lowmc_instance(unsigned int idx) {
  if (lowmc_instances_initialized[idx]) {
#if defined(MUL_M4RI)
    lowmc_clear(lowmc_instances[idx]);
#endif
    lowmc_instances_initialized[idx] = false;
  }
}

/**
 * For trying other LowMC instances, call this function with
 * PARAMETER_SET_INVALID and a custom LowMC instance and the transform as last arguments.
 */
#if defined(WITH_CUSTOM_INSTANCES)
static bool create_instance(picnic_instance_t* pp, picnic_params_t param, const lowmc_t* lowmc,
                            transform_t transform) {
#else
static bool create_instance(picnic_instance_t* pp, picnic_params_t param) {
#endif
  const lowmc_t* lowmc_instance = NULL;

#if defined(WITH_CUSTOM_INSTANCES)
  uint32_t pq_security_level, num_rounds;
#endif
  switch (param) {
  case Picnic_L1_FS:
  case Picnic_L1_UR:
    lowmc_instance = lowmc_get_instance(0);
    break;

  case Picnic_L3_FS:
  case Picnic_L3_UR:
    lowmc_instance = lowmc_get_instance(1);
    break;

  case Picnic_L5_FS:
  case Picnic_L5_UR:
    lowmc_instance = lowmc_get_instance(2);
    break;

#if defined(WITH_CUSTOM_INSTANCES)
  case PARAMETER_SET_INVALID:
    if (!lowmc) {
      return false;
    }

    pq_security_level = lowmc->n / 2;
    num_rounds        = ceil(lowmc->n / 0.5849625007211562);
    break;
#endif

  default:
    return false;
  }

  if (!lowmc_instance) {
#if defined(WITH_CUSTOM_INSTANCES)
    pp->lowmc = lowmc;
#else
    return false;
#endif
  }

  pp->lowmc_impl              = lowmc_get_implementation(pp->lowmc);
  pp->lowmc_store_impl        = lowmc_store_get_implementation(pp->lowmc);
  pp->zkbpp_lowmc_impl        = get_zkbpp_lowmc_implementation(pp->lowmc);
  pp->zkbpp_lowmc_verify_impl = get_zkbpp_lowmc_verify_implementation(pp->lowmc);
#if defined(WITH_CUSTOM_INSTANCES)
  if (!lowmc_instance) {
    pp->params    = param;
    pp->transform = transform;

    const uint32_t digest_size = MAX(32, (4 * pq_security_level + 7) / 8);
    const uint32_t seed_size   = (2 * pq_security_level + 7) / 8;
    pp->digest_size            = digest_size;
    pp->seed_size              = seed_size;
    pp->num_rounds             = num_rounds;

    // bytes required to store one input share
    pp->input_size = (pp->lowmc->k + 7) >> 3;
    // bytes required to store one output share
    pp->output_size = (pp->lowmc->n + 7) >> 3;
    // number of bits per view per LowMC round
    pp->view_round_size = pp->lowmc->m * 3;
    // bytes required to store communicated bits (i.e. views) of one round
    pp->view_size = (pp->view_round_size * pp->lowmc->r + 7) >> 3;

    // bytes required to store collapsed challenge
    pp->collapsed_challenge_size = (pp->num_rounds + 3) >> 2;

    if (pp->transform == TRANSFORM_UR) {
      pp->unruh_without_input_bytes_size = pp->seed_size + pp->view_size;
      pp->unruh_with_input_bytes_size    = pp->unruh_without_input_bytes_size + pp->input_size;
    } else {
      pp->unruh_without_input_bytes_size = pp->unruh_with_input_bytes_size = 0;
    }

    // we can use unruh_without_input_bytes_size here. In call cases where we need
    // to write more, we do not need to write the input share
    const uint32_t per_round_size = pp->input_size + pp->view_size + pp->digest_size +
                                    2 * pp->seed_size + pp->unruh_without_input_bytes_size;
    pp->max_signature_size = pp->collapsed_challenge_size + pp->num_rounds * per_round_size;
  }
#endif

  return true;
}

const picnic_instance_t* picnic_instance_get(picnic_params_t param) {
  if (param <= PARAMETER_SET_INVALID || param >= PARAMETER_SET_MAX_INDEX) {
    return NULL;
  }

  if (!instance_initialized[param]) {
#if defined(WITH_CUSTOM_INSTANCES)
    if (!create_instance(&instances[param], param, NULL, TRANSFORM_FS)) {
#else
    if (!create_instance(&instances[param], param)) {
#endif
      return NULL;
    }
    instance_initialized[param] = true;
  }

  return &instances[param];
}

ATTR_DTOR static void clear_instances(void) {
  for (unsigned int p = PARAMETER_SET_INVALID + 1; p < PARAMETER_SET_MAX_INDEX; ++p) {
    if (instance_initialized[p]) {
      instance_initialized[p] = false;
    }
  }

  for (unsigned int i = 0;
       i < sizeof(lowmc_instances_initialized) / sizeof(lowmc_instances_initialized[0]); ++i) {
    clear_lowmc_instance(i);
  }
}
