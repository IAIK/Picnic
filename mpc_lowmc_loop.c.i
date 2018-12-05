/*
 *  This file is part of the optimized implementation of the Picnic signature scheme.
 *  See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */

lowmc_round_t const* round = lowmc->rounds;
#if defined(REDUCED_ROUND_KEY_COMPUTATION)
  mzd_local_t* nl_part[reduced_shares];
  mzd_local_init_multiple_ex(nl_part, reduced_shares, 1, (LOWMC_R)*32, false);
#if defined(OPTIMIZED_LINEAR_LAYER_EVALUATION)
  MPC_LOOP_CONST_C(XOR, x, x, lowmc->precomputed_constant_linear, reduced_shares, ch);
  MPC_LOOP_CONST(MUL_MC, nl_part, lowmc_key,
                 CONCAT(lowmc->precomputed_non_linear_part, matrix_postfix), reduced_shares);
  MPC_LOOP_CONST_C(XOR_MC, nl_part, nl_part, lowmc->precomputed_constant_non_linear, reduced_shares, ch);
  for (unsigned i = 0; i < (LOWMC_R-1); ++i, ++views, ++round) {
    RANDTAPE;
#if defined(RECOVER_FROM_STATE)
    RECOVER_FROM_STATE(x, i);
#endif
    SBOX(SBOX_ARGS, sbox, y, x, views, r, &lowmc->mask, LOWMC_N, shares);
    for (unsigned int k = 0; k < reduced_shares; ++k) {
#if defined(M_FIXED_10)
      const word nl = CONST_FIRST_ROW(nl_part[k])[i >> 1];
      FIRST_ROW(y[k])[(LOWMC_N) / (sizeof(word) * 8) - 1] ^=
        (i & 1) ? (nl & WORD_C(0xFFFFFFFF00000000)) : (nl << 32);
#elif defined(M_FIXED_1)
      const word nl = CONST_FIRST_ROW(nl_part[k])[i / 21];
      FIRST_ROW(y[k])[(LOWMC_N) / (sizeof(word) * 8) - 1] ^= (nl << ((20-(i%21))*3)) & WORD_C(0xE000000000000000);
#else
#error "RLL only works with 1 or 10 Sboxes atm"
#endif
    }
    MPC_LOOP_CONST(MUL_Z, x, y, CONCAT(round->z, matrix_postfix), reduced_shares);

#if defined(WITH_AVX2)
    for(unsigned int k = 0; k < reduced_shares; ++k) {
#if defined(M_FIXED_10)
      mzd_shuffle_pext_30(y[k], round->r_mask);
#elif defined(M_FIXED_1)
      mzd_shuffle_pext_3(y[k], round->r_mask);
#else
#error "RLL only works with 1 or 10 Sboxes atm"
#endif
    }
#else
    //shuffle x correctly (in-place), slow and probably stupid version
    for(unsigned j = round->num_fixes; j; j--) {
      for(unsigned l = round->r_cols[j-1]; l < LOWMC_N - 1 - (3*LOWMC_M-j); l++) {
        for(unsigned int k = 0; k < reduced_shares; ++k) {
          //swap bits
          word a = (FIRST_ROW(y[k])[l / (sizeof(word) * 8)] >> (l % (sizeof(word) * 8))) & WORD_C(0x1);
          word b = (FIRST_ROW(y[k])[(l+1) / (sizeof(word) * 8)] >> ((l+1) % (sizeof(word) * 8))) & WORD_C(0x1);
          word xx = a ^ b;
          FIRST_ROW(y[k])[l / (sizeof(word) * 8)] ^=  xx << (l % (sizeof(word) * 8));
          FIRST_ROW(y[k])[(l+1) / (sizeof(word) * 8)] ^=  xx << ((l+1) % (sizeof(word) * 8));
        }
      }
    }
#endif

    MPC_LOOP_CONST(MUL_R, x, y, CONCAT(round->r, matrix_postfix), reduced_shares);
    for(unsigned int k = 0; k < reduced_shares; ++k) {
#if defined(M_FIXED_10)
      FIRST_ROW(y[k])[(LOWMC_N) / (sizeof(word) * 8) - 1] &= WORD_C(0x00000003FFFFFFFF); //clear nl part
#elif defined(M_FIXED_1)
      FIRST_ROW(y[k])[(LOWMC_N) / (sizeof(word) * 8) - 1] &= WORD_C(0x1FFFFFFFFFFFFFFF); //clear nl part
#else
#error "RLL only works with 1 or 10 Sboxes atm"
#endif
    }
    MPC_LOOP_SHARED(XOR, x, x, y, reduced_shares);
  }
  unsigned i = (LOWMC_R-1);
  RANDTAPE;
#if defined(RECOVER_FROM_STATE)
  RECOVER_FROM_STATE(x, i);
#endif
  SBOX(SBOX_ARGS, sbox, y, x, views, r, &lowmc->mask, LOWMC_N, shares);

  for (unsigned int k = 0; k < reduced_shares; ++k) {
#if defined(M_FIXED_10)
    const word nl = CONST_FIRST_ROW(nl_part[k])[i >> 1];
    FIRST_ROW(y[k])[(LOWMC_N) / (sizeof(word) * 8) - 1] ^=
        (i & 1) ? (nl & WORD_C(0xFFFFFFFF00000000)) : (nl << 32);
#elif defined(M_FIXED_1)
    const word nl = CONST_FIRST_ROW(nl_part[k])[i / 21];
    FIRST_ROW(y[k])[(LOWMC_N) / (sizeof(word) * 8) - 1] ^= (nl << ((20-(i%21))*3)) & WORD_C(0xE000000000000000);
#else
#error "RLL only works with 1 or 10 Sboxes atm"
#endif
  }
  MPC_LOOP_CONST(MUL, x, y, CONCAT(lowmc->zr, matrix_postfix), reduced_shares);
#else
  MPC_LOOP_CONST_C(XOR, x, x, lowmc->precomputed_constant_linear, reduced_shares, ch);
  MPC_LOOP_CONST(MUL_MC, nl_part, lowmc_key,
                 CONCAT(lowmc->precomputed_non_linear_part, matrix_postfix), reduced_shares);
  MPC_LOOP_CONST_C(XOR_MC, nl_part, nl_part, lowmc->precomputed_constant_non_linear, reduced_shares, ch);
  for (unsigned i = 0; i < (LOWMC_R); ++i, ++views, ++round) {
    RANDTAPE;
#if defined(RECOVER_FROM_STATE)
    RECOVER_FROM_STATE(x, i);
#endif
    SBOX(SBOX_ARGS, sbox, y, x, views, r, &lowmc->mask, LOWMC_N, shares);
    for (unsigned int k = 0; k < reduced_shares; ++k) {
#if defined(M_FIXED_10)
      const word nl = CONST_FIRST_ROW(nl_part[k])[i >> 1];
      FIRST_ROW(y[k])[(LOWMC_N) / (sizeof(word) * 8) - 1] ^=
        (i & 1) ? (nl & WORD_C(0xFFFFFFFF00000000)) : (nl << 32);
#elif defined(M_FIXED_1)
      const word nl = CONST_FIRST_ROW(nl_part[k])[i / 21];
      FIRST_ROW(y[k])[(LOWMC_N) / (sizeof(word) * 8) - 1] ^= (nl << ((20-(i%21))*3)) & WORD_C(0xE000000000000000);
#else
#error "RLL only works with 1 or 10 Sboxes atm"
#endif
    }
    MPC_LOOP_CONST(MUL, x, y, CONCAT(round->l, matrix_postfix), reduced_shares);
  }
#endif
  mzd_local_free_multiple(nl_part);
#else
for (unsigned i = 0; i < (LOWMC_R); ++i, ++views, ++round) {
  RANDTAPE;
#if defined(RECOVER_FROM_STATE)
  RECOVER_FROM_STATE(x, i);
#endif
  SBOX(SBOX_ARGS, sbox, y, x, views, r, &lowmc->mask, LOWMC_N, shares);
  MPC_LOOP_CONST(MUL, x, y, CONCAT(round->l, matrix_postfix), reduced_shares);
  MPC_LOOP_CONST_C(XOR, x, x, round->constant, reduced_shares, ch);
  MPC_LOOP_CONST(ADDMUL, x, lowmc_key, CONCAT(round->k, matrix_postfix), reduced_shares);
}
#endif
#if defined(RECOVER_FROM_STATE)
RECOVER_FROM_STATE(x, LOWMC_R);
#endif

// vim: ft=c
