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

#if defined(PICNIC2_AUX_COMPUTATION)
  MUL(x, lowmc_key, LOWMC_INSTANCE.k0_matrix);
#else
  COPY(x, p);
  ADDMUL(x, lowmc_key, LOWMC_INSTANCE.k0_matrix);
#endif

  lowmc_round_t const* round = LOWMC_INSTANCE.rounds;
  for (unsigned i = 0; i < LOWMC_R; ++i, ++round) {
#if defined(RECORD_STATE)
    COPY(state->state[i], x);
#endif
#if defined(PICNIC2_AUX_COMPUTATION)
    SBOX(x, tapes);
#else
    SBOX(x);
#endif

    MUL(y, x, round->l_matrix);
#if !defined(PICNIC2_AUX_COMPUTATION)
    XOR(x, y, round->constant);
#else
    COPY(x, y);
#endif
    ADDMUL(x, lowmc_key, round->k_matrix);
  }

#if !defined(PICNIC2_AUX_COMPUTATION)
#if defined(RECORD_STATE)
  COPY(state->state[LOWMC_R], x);
#else
  COPY(c, x);
#endif
#endif
}

// vim: ft=c
