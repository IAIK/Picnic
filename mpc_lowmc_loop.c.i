/*
 *  This file is part of the optimized implementation of the Picnic signature scheme.
 *  See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */

lowmc_round_t const* round = LOWMC_INSTANCE.rounds;
for (unsigned i = 0; i < (LOWMC_R); ++i, ++views, ++round) {
  RANDTAPE;
#if defined(RECOVER_FROM_STATE)
  RECOVER_FROM_STATE(x, i);
#endif
  SBOX(sbox, y, x, views, r, LOWMC_N, shares, reduced_shares);
  MPC_LOOP_CONST(MUL, x, y, round->l_matrix, reduced_shares);
  MPC_LOOP_CONST_C(XOR, x, x, round->constant, reduced_shares, ch);
  MPC_LOOP_CONST(ADDMUL, x, lowmc_key, round->k_matrix, reduced_shares);
}
#if defined(RECOVER_FROM_STATE)
RECOVER_FROM_STATE(x, LOWMC_R);
#endif

// vim: ft=c
