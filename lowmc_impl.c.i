/*
 *  This file is part of the optimized implementation of the Picnic signature scheme.
 *  See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */

#if defined(LOWMC_INSTANCE)
#define lowmc LOWMC_INSTANCE
#else
#define lowmc lowmc_instance
#endif

#if defined(RECORD_STATE)
static void N_LOWMC(lowmc_t const* lowmc_instance, lowmc_key_t const* lowmc_key,
                            mzd_local_t const* p, recorded_state_t* state) {
#else
static mzd_local_t* N_LOWMC(lowmc_t const* lowmc_instance, lowmc_key_t const* lowmc_key,
                            mzd_local_t const* p) {
#endif
#if defined(LOWMC_INSTANCE)
  (void)lowmc_instance;
#endif
#if defined(REDUCED_LINEAR_LAYER)
#if defined(REDUCED_LINEAR_LAYER_NEXT)
  mzd_local_t* x       = mzd_local_init_ex(1, LOWMC_N, false);
  mzd_local_t* y       = mzd_local_init_ex(1, LOWMC_N, false);
  mzd_local_t* x_nl    = mzd_local_init_ex(1, 3*LOWMC_M, false);
//  mzd_local_t* x_l     = mzd_local_init_ex(1, LOWMC_N-3*LOWMC_M, false);
  mzd_local_t* nl_part = mzd_local_init_ex(1, LOWMC_R * 32, false);

  XOR(y, p, lowmc->precomputed_constant_linear);
  ADDMUL(y, lowmc_key, CONCAT(lowmc->k0, matrix_postfix));
  MUL_MC(nl_part, lowmc_key, CONCAT(lowmc->precomputed_non_linear_part, matrix_postfix));
  XOR_MC(nl_part, nl_part, lowmc->precomputed_constant_non_linear);

  //multiply non-linear part of state with Z0 matrix
  MUL(x, y, CONCAT(lowmc->z0, matrix_postfix));

  lowmc_round_t const* round = lowmc->rounds;
  for (unsigned i = 0; i < LOWMC_R; ++i, ++round) {
    SBOX(x, &lowmc->mask);

    const word nl = CONST_FIRST_ROW(nl_part)[i >> 1];
    FIRST_ROW(x)
    [(LOWMC_N) / (sizeof(word) * 8) - 1] ^=
        (i & 1) ? (nl & WORD_C(0xFFFFFFFF00000000)) : (nl << 32);

    //hardcoded masks for 10-sbox case
    const word nl_mask = WORD_C(0xFFFFFFFC00000000);
    const word inv_nl_mask = WORD_C(0x00000003FFFFFFFF);
    FIRST_ROW(x_nl)[x_nl->width - 1] = (FIRST_ROW(x)[x->width - 1] & nl_mask) >> (sizeof(word)*8-3*LOWMC_M);
    MUL(y, x_nl, CONCAT(round->z, matrix_postfix));

    //copy x to x_l
    FIRST_ROW(x)[x->width-1] &= inv_nl_mask;
//    for(unsigned j = 0; j < x_l->width; j++)
//        FIRST_ROW(x_l)[j] = FIRST_ROW(x)[j];

//    MUL(x_nl, x_l, CONCAT(round->a, matrix_postfix));
    mzd_mul_v_popcnt(x_nl, x, CONCAT(round->aT, matrix_postfix));
    FIRST_ROW(y)[y->width - 1] ^= (FIRST_ROW(x_nl)[x_nl->width - 1] << (sizeof(word)*8-3*LOWMC_M)) & nl_mask;
    XOR(x, x, y);
  }

  mzd_local_free(y);
//  mzd_local_free(x_l);
  mzd_local_free(x_nl);
  mzd_local_free(nl_part);
  return x;
#else
  mzd_local_t* x       = mzd_local_init_ex(1, LOWMC_N, false);
  mzd_local_t* y       = mzd_local_init_ex(1, LOWMC_N, false);
  mzd_local_t* nl_part = mzd_local_init_ex(1, LOWMC_R * 32, false);

  XOR(x, p, lowmc->precomputed_constant_linear);
  ADDMUL(x, lowmc_key, CONCAT(lowmc->k0, matrix_postfix));
  MUL_MC(nl_part, lowmc_key, CONCAT(lowmc->precomputed_non_linear_part, matrix_postfix));
  XOR_MC(nl_part, nl_part, lowmc->precomputed_constant_non_linear);

  lowmc_round_t const* round = lowmc->rounds;
  for (unsigned i = 0; i < LOWMC_R; ++i, ++round) {
#if defined(RECORD_STATE)
    mzd_local_copy(state->state[i], x);
#endif
    SBOX(x, &lowmc->mask);

    const word nl = CONST_FIRST_ROW(nl_part)[i >> 1];
    FIRST_ROW(x)
    [(LOWMC_N) / (sizeof(word) * 8) - 1] ^=
        (i & 1) ? (nl & WORD_C(0xFFFFFFFF00000000)) : (nl << 32);
    MUL(y, x, CONCAT(round->l, matrix_postfix));

    // swap x and y
    mzd_local_t* t = x;
    x              = y;
    y              = t;
  }

  mzd_local_free(nl_part);
  mzd_local_free(y);
#if defined(RECORD_STATE)
  mzd_local_copy(state->state[LOWMC_R], x);
  mzd_local_free(x);
#else
  return x;
#endif
#else
  mzd_local_t* x = mzd_local_init_ex(1, LOWMC_N, false);
  mzd_local_t* y = mzd_local_init_ex(1, LOWMC_N, false);

  mzd_local_copy(x, p);
  ADDMUL(x, lowmc_key, CONCAT(lowmc->k0, matrix_postfix));

  lowmc_round_t const* round = lowmc->rounds;
  for (unsigned i = 0; i < LOWMC_R; ++i, ++round) {
#if defined(RECORD_STATE)
    mzd_local_copy(state->state[i], x);
#endif
    SBOX(x, &lowmc->mask);

    MUL(y, x, CONCAT(round->l, matrix_postfix));
    XOR(x, y, round->constant);
    ADDMUL(x, lowmc_key, CONCAT(round->k, matrix_postfix));
  }

  mzd_local_free(y);
#if defined(RECORD_STATE)
  mzd_local_copy(state->state[LOWMC_R], x);
  mzd_local_free(x);
#else
  return x;
#endif
#endif
}

#undef lowmc

// vim: ft=c
