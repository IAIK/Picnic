/*
 *  This file is part of the optimized implementation of the Picnic signature scheme.
 *  See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */


#if defined(FN_ATTR)
FN_ATTR
#endif
#if defined(PICNIC2_AUX_COMPUTATION)
static void N_LOWMC(lowmc_key_t const* lowmc_key, randomTape_t* tapes) {
#else
#if defined(RECORD_STATE)
static void N_LOWMC(lowmc_key_t const* lowmc_key, mzd_local_t const* p, recorded_state_t* state) {
#else
static void N_LOWMC(lowmc_key_t const* lowmc_key, mzd_local_t const* p, mzd_local_t* c) {
#endif
#endif
  mzd_local_t x[((LOWMC_N) + 255) / 256];
  mzd_local_t y[((LOWMC_N) + 255) / 256];
mzd_local_t nl_part[(LOWMC_R * 32 + 255) / 256];

#if defined(PICNIC2_AUX_COMPUTATION)
  MUL(x, lowmc_key, LOWMC_INSTANCE.k0_matrix);
  MUL_MC(nl_part, lowmc_key, LOWMC_INSTANCE.precomputed_non_linear_part_matrix);
#else
  XOR(x, p, LOWMC_INSTANCE.precomputed_constant_linear);
  ADDMUL(x, lowmc_key, LOWMC_INSTANCE.k0_matrix);
  MUL_MC(nl_part, lowmc_key, LOWMC_INSTANCE.precomputed_non_linear_part_matrix);
  XOR_MC(nl_part, nl_part, LOWMC_INSTANCE.precomputed_constant_non_linear);
#endif

  // multiply non-linear part of state with Z0 matrix

  lowmc_round_t const* round = LOWMC_INSTANCE.rounds;
  for (unsigned i = 0; i < LOWMC_R - 1; ++i, ++round) {
#if defined(RECORD_STATE)
    COPY(state->state[i], x);
#endif
#if defined(PICNIC2_AUX_COMPUTATION)
    SBOX(x, tapes);
#else
    SBOX(x);
#endif

    const word nl = CONST_BLOCK(nl_part, i >> 3)->w64[(i & 0x7) >> 1];
    BLOCK(x, 0)->w64[(LOWMC_N) / (sizeof(word) * 8) - 1] ^=
        (nl << (1 - (i & 1)) * 32) & WORD_C(0xFFFFFFFF00000000);

    MUL_Z(y, x, round->z_matrix);
    SHUFFLE(x, round->r_mask);
    ADDMUL_R(y, x, round->r_matrix);

    BLOCK(x, 0)->w64[(LOWMC_N) / (sizeof(word) * 8) - 1] &=
        WORD_C(0x00000003FFFFFFFF); // clear nl part
    XOR(x, y, x);
  }
#if defined(RECORD_STATE)
  COPY(state->state[LOWMC_R - 1], x);
#endif
#if defined(PICNIC2_AUX_COMPUTATION)
  SBOX(x, tapes);
#else
  SBOX(x);

  unsigned int i = (LOWMC_R - 1);
  const word nl  = CONST_BLOCK(nl_part, i >> 3)->w64[(i & 0x7) >> 1];
  BLOCK(x, 0)->w64[(LOWMC_N) / (sizeof(word) * 8) - 1] ^=
      (nl << (1 - (i & 1)) * 32) & WORD_C(0xFFFFFFFF00000000);
  MUL(y, x, LOWMC_INSTANCE.zr_matrix);
  COPY(x, y);
#endif

#if !defined(PICNIC2_AUX_COMPUTATION)
#if defined(RECORD_STATE)
  COPY(state->state[LOWMC_R], x);
#else
  COPY(c, x);
#endif
#endif
}

// vim: ft=c
