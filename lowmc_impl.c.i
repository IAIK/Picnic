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
  mzd_local_t* x       = mzd_local_init_ex(1, LOWMC_N, false);
  mzd_local_t* y       = mzd_local_init_ex(1, LOWMC_N, false);
  mzd_local_t* tmp     = mzd_local_init_ex(1, LOWMC_N, false);
#if defined(M_FIXED_10)
  mzd_local_t* nl_part = mzd_local_init_ex(1, LOWMC_R * 32, false);
#elif defined(M_FIXED_1)
  mzd_local_t* nl_part = mzd_local_init_ex(1, ((LOWMC_R + 20) / 21) * 64, false);
#else
  #error "RLL only works with 1 or 10 Sboxes atm"
#endif

#if defined(REDUCED_LINEAR_LAYER_NEXT)
  XOR(x, p, lowmc->precomputed_constant_linear);
  ADDMUL(x, lowmc_key, CONCAT(lowmc->k0, matrix_postfix));
  MUL_MC(nl_part, lowmc_key, CONCAT(lowmc->precomputed_non_linear_part, matrix_postfix));
  XOR_MC(nl_part, nl_part, lowmc->precomputed_constant_non_linear);

  //multiply non-linear part of state with Z0 matrix

  lowmc_round_t const* round = lowmc->rounds;
  for (unsigned i = 0; i < LOWMC_R-1; ++i, ++round) {
#if defined(RECORD_STATE)
    mzd_local_copy(state->state[i], x);
#endif
    SBOX(x, &lowmc->mask);

#if defined(M_FIXED_10)
    const word nl = CONST_FIRST_ROW(nl_part)[i >> 1];
    FIRST_ROW(x)
    [(LOWMC_N) / (sizeof(word) * 8) - 1] ^= (nl << (1-(i&1))*32) & WORD_C(0xFFFFFFFF00000000);
#elif defined(M_FIXED_1)
    const word nl = CONST_FIRST_ROW(nl_part)[i / 21];
    FIRST_ROW(x)[(LOWMC_N) / (sizeof(word) * 8) - 1] ^= (nl << ((20-(i%21))*3)) & WORD_C(0xE000000000000000);
#else
#error "RLL only works with 1 or 10 Sboxes atm"
#endif

    MUL_Z(y, x, CONCAT(round->z, matrix_postfix));
    //shuffle x correctly (in-place), slow and probably stupid version
    for(unsigned j = round->num_fixes; j; j--) {
        for(unsigned k = round->r_cols[j-1]; k < LOWMC_N - 1 - (3*LOWMC_M-j); k++) {
            //swap bits
            word a = (FIRST_ROW(x)[k / (sizeof(word) * 8)] >> (k % (sizeof(word) * 8))) & WORD_C(0x1);
            word b = (FIRST_ROW(x)[(k+1) / (sizeof(word) * 8)] >> ((k+1) % (sizeof(word) * 8))) & WORD_C(0x1);
            word xx = a ^ b;
            FIRST_ROW(x)[k / (sizeof(word) * 8)] ^=  xx << (k % (sizeof(word) * 8));
            FIRST_ROW(x)[(k+1) / (sizeof(word) * 8)] ^=  xx << ((k+1) % (sizeof(word) * 8));
        }
    }
    MUL_R(y, x, CONCAT(round->r, matrix_postfix));

#if defined(M_FIXED_10)
    FIRST_ROW(x)[(LOWMC_N) / (sizeof(word) * 8) - 1] &= WORD_C(0x00000003FFFFFFFF); //clear nl part
#elif defined(M_FIXED_1)
    FIRST_ROW(x)[(LOWMC_N) / (sizeof(word) * 8) - 1] &= WORD_C(0x1FFFFFFFFFFFFFFF); //clear nl part
#else
#error "RLL only works with 1 or 10 Sboxes atm"
#endif
    XOR(x, y, x);
//    mzd_local_t* t = x;
//    x              = y;
//    y              = t;

  }
#if defined(RECORD_STATE)
  mzd_local_copy(state->state[LOWMC_R-1], x);
#endif
  SBOX(x, &lowmc->mask);

  unsigned i = (LOWMC_R-1);
#if defined(M_FIXED_10)
  const word nl = CONST_FIRST_ROW(nl_part)[i >> 1];
  FIRST_ROW(x)
  [(LOWMC_N) / (sizeof(word) * 8) - 1] ^= (nl << (1-(i&1))*32) & WORD_C(0xFFFFFFFF00000000);
#elif defined(M_FIXED_1)
  const word nl = CONST_FIRST_ROW(nl_part)[i / 21];
  FIRST_ROW(x)[(LOWMC_N) / (sizeof(word) * 8) - 1] ^= (nl << ((20-(i%21))*3)) & WORD_C(0xE000000000000000);
#else
#error "RLL only works with 1 or 10 Sboxes atm"
#endif
  MUL(y,x,CONCAT(lowmc->zr, matrix_postfix));
  mzd_local_t* t = x;
  x              = y;
  y              = t;
#else
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

#if defined(M_FIXED_10)
    const word nl = CONST_FIRST_ROW(nl_part)[i >> 1];
    FIRST_ROW(x)
    [(LOWMC_N) / (sizeof(word) * 8) - 1] ^=
        (i & 1) ? (nl & WORD_C(0xFFFFFFFF00000000)) : (nl << 32);
#elif defined(M_FIXED_1)
    const word nl = CONST_FIRST_ROW(nl_part)[i / 21];
    FIRST_ROW(x)[(LOWMC_N) / (sizeof(word) * 8) - 1] ^= (nl << ((20-(i%21))*3)) & WORD_C(0xE000000000000000);
#else
#error "RLL only works with 1 or 10 Sboxes atm"
#endif
    MUL(y, x, CONCAT(round->l, matrix_postfix));
    // swap x and y
    mzd_local_t* t = x;
    x              = y;
    y              = t;
  }
#endif
  mzd_local_free(nl_part);
  mzd_local_free(tmp);
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
