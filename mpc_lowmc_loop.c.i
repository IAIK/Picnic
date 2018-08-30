/*
 *  This file is part of the optimized implementation of the Picnic signature scheme.
 *  See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */

lowmc_round_t const* round = lowmc->rounds;
#if defined(REDUCED_LINEAR_LAYER)
  mzd_local_t* nl_part[reduced_shares];
  mzd_local_init_multiple_ex(nl_part, reduced_shares, 1, (LOWMC_R)*32, false);
#if defined(REDUCED_LINEAR_LAYER_NEXT)
  MPC_LOOP_CONST_C(XOR, x, x, lowmc->precomputed_constant_linear, reduced_shares, ch);
  MPC_LOOP_CONST(MUL_MC, nl_part, lowmc_key,
                 CONCAT(lowmc->precomputed_non_linear_part, matrix_postfix), reduced_shares);
  MPC_LOOP_CONST_C(XOR_MC, nl_part, nl_part, lowmc->precomputed_constant_non_linear, reduced_shares, ch);
  MPC_LOOP_CONST(MUL, y, x, CONCAT(lowmc->z0, matrix_postfix), reduced_shares);
  mpc_copy(x, y, reduced_shares); //TODO: avoid copy?
  for (unsigned i = 0; i < (LOWMC_R); ++i, ++views, ++round) {
    RANDTAPE;
#if defined(RECOVER_FROM_STATE)
    RECOVER_FROM_STATE(x, i);
#endif
    SBOX(SBOX_ARGS, sbox, y, x, views, r, &lowmc->mask, &vars, LOWMC_N, shares);
    for (unsigned int k = 0; k < reduced_shares; ++k) {
      const word nl = CONST_FIRST_ROW(nl_part[k])[i >> 1];
      FIRST_ROW(y[k])
      [(LOWMC_N) / (sizeof(word) * 8) - 1] ^=
          (i & 1) ? (nl & WORD_C(0xFFFFFFFF00000000)) : (nl << 32);
    }
    MPC_LOOP_CONST(MUL_Z, x, y, CONCAT(round->z, matrix_postfix), reduced_shares);
    MPC_LOOP_CONST(MUL_A, y, y, CONCAT(round->aT, matrix_postfix), reduced_shares);
    MPC_LOOP_SHARED(XOR, x, x, y, shares);
  }
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
    SBOX(SBOX_ARGS, sbox, y, x, views, r, &lowmc->mask, &vars, LOWMC_N, shares);
    for (unsigned int k = 0; k < reduced_shares; ++k) {
      const word nl = CONST_FIRST_ROW(nl_part[k])[i >> 1];
      FIRST_ROW(y[k])
      [(LOWMC_N) / (sizeof(word) * 8) - 1] ^=
          (i & 1) ? (nl & WORD_C(0xFFFFFFFF00000000)) : (nl << 32);
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
  SBOX(SBOX_ARGS, sbox, y, x, views, r, &lowmc->mask, &vars, LOWMC_N, shares);
  MPC_LOOP_CONST(MUL, x, y, CONCAT(round->l, matrix_postfix), reduced_shares);
  MPC_LOOP_CONST_C(XOR, x, x, round->constant, reduced_shares, ch);
  MPC_LOOP_CONST(ADDMUL, x, lowmc_key, CONCAT(round->k, matrix_postfix), reduced_shares);
}
#endif
#if defined(RECOVER_FROM_STATE)
RECOVER_FROM_STATE(x, LOWMC_R);
#endif

// vim: ft=c
