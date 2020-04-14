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
static void N_LOWMC(lowmc_key_t* lowmc_key, randomTape_t* tapes) {
  mzd_local_t x[((LOWMC_N) + 255) / 256];
  mzd_local_t y[((LOWMC_N) + 255) / 256];
  mzd_local_t key0[((LOWMC_N) + 255) / 256];
  uint8_t temp[32] = { 0 };

  const size_t state_size_bytes = (LOWMC_N + 7) / 8;

  mzd_from_char_array(x, temp, state_size_bytes);
  COPY(key0, lowmc_key);
  MUL(lowmc_key, key0, LOWMC_INSTANCE.ki0_matrix);

  lowmc_round_t const* round = &LOWMC_INSTANCE.rounds[LOWMC_R - 1];
  for (unsigned r = 0; r < LOWMC_R; ++r, round--) {
    ADDMUL(x, lowmc_key, round->k_matrix);
    MUL(y, x, round->li_matrix);

    // recover input masks from tapes, only in first round we use the key as input
    if (r == LOWMC_R - 1) {
      COPY(x, key0);
    } else {
      tapes->pos = LOWMC_N * 2 * (LOWMC_R - 1 - r);
      memset(temp, 0, 32);
      for (size_t j = 0; j < LOWMC_N; j++) {
        setBit(temp, j, getBit(tapes->parity_tapes, tapes->pos + j));
      }
      mzd_from_char_array(x, temp, state_size_bytes);
    }
    tapes->pos = LOWMC_N * 2 * (LOWMC_R - 1 - r) + LOWMC_N;
    SBOX(x, y, tapes);
  }
}

// vim: ft=c
