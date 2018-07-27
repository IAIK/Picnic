  lowmc_round_t const* round = lowmc->rounds;

#if defined(REDUCED_LINEAR_LAYER)
  MPC_LOOP_CONST_C(XOR, x, x, lowmc->precomputed_constant_linear, shares, ch);
  mzd_local_t* nl_part[shares];
  mzd_local_init_multiple_ex(nl_part, shares, 1, (LOWMC_R)*32, false);
  MPC_LOOP_CONST(MUL_MC, nl_part, lowmc_key,
                 CONCAT(lowmc->precomputed_non_linear_part, matrix_postfix), shares);
  MPC_LOOP_CONST_C(XOR_MC, nl_part, nl_part, lowmc->precomputed_constant_non_linear, shares, ch);
  for (unsigned i = 0; i < (LOWMC_R); ++i, ++views, ++round) {
    RANDTAPE;
    SBOX(SBOX_ARGS, sbox, y, x, views, r, &lowmc->mask, &vars, LOWMC_N, shares);
    for (unsigned int k = 0; k < shares; ++k) {
      const word nl = CONST_FIRST_ROW(nl_part[k])[i >> 1];
      FIRST_ROW(y[k])
      [(LOWMC_N) / (sizeof(word) * 8) - 1] ^=
          (i & 1) ? (nl & WORD_C(0xFFFFFFFF00000000)) : (nl << 32);
    }
    MPC_LOOP_CONST(MUL, x, y, CONCAT(round->l, matrix_postfix), shares);
  }
  mzd_local_free_multiple(nl_part);
#else
  for (unsigned i = 0; i < (LOWMC_R); ++i, ++views, ++round) {
    RANDTAPE;
    SBOX(SBOX_ARGS, sbox, y, x, views, r, &lowmc->mask, &vars, LOWMC_N, shares);
    MPC_LOOP_CONST(MUL, x, y, CONCAT(round->l, matrix_postfix), shares);
    MPC_LOOP_CONST_C(XOR, x, x, round->constant, shares, ch);
    MPC_LOOP_CONST(ADDMUL, x, lowmc_key, CONCAT(round->k, matrix_postfix), shares);
  }
#endif

// vim: ft=c
