#if defined(REDUCED_LINEAR_LAYER)
#define loop_impl(sbox_args, sbox, sbox_selector, ch, shares)                                      \
  MPC_LOOP_CONST_C(XOR, x, x, lowmc->precomputed_constant_linear, shares, ch);                     \
  mzd_local_t* nl_part[shares];                                                                    \
  mzd_local_init_multiple_ex(nl_part, shares, 1, (LOWMC_R)*32, false);                             \
  MPC_LOOP_CONST(MUL_MC, nl_part, lowmc_key,                                                       \
                 CONCAT(lowmc->precomputed_non_linear_part, matrix_postfix), shares);              \
  MPC_LOOP_CONST_C(XOR_MC, nl_part, nl_part, lowmc->precomputed_constant_non_linear, shares, ch);  \
  for (unsigned i = 0; i < (LOWMC_R); ++i, ++views, ++round) {                                     \
    R(sbox_selector);                                                                              \
    SBOX(sbox_args, sbox, sbox_selector, y, x, views, r, &lowmc->mask, &vars, LOWMC_N, shares);    \
    for (unsigned int k = 0; k < shares; ++k) {                                                    \
      const word nl = CONST_FIRST_ROW(nl_part[k])[i >> 1];                                         \
      FIRST_ROW(y[k])                                                                              \
      [(LOWMC_N) / (sizeof(word) * 8) - 1] ^=                                                      \
          (i & 1) ? (nl & WORD_C(0xFFFFFFFF00000000)) : (nl << 32);                                \
    }                                                                                              \
    MPC_LOOP_CONST(MUL, x, y, CONCAT(round->l, matrix_postfix), shares);                           \
  }                                                                                                \
  mzd_local_free_multiple(nl_part);
#else
#define loop_impl(sbox_args, sbox, sbox_selector, ch, shares)                                      \
  for (unsigned i = 0; i < (LOWMC_R); ++i, ++views, ++round) {                                     \
    R(sbox_selector);                                                                              \
    SBOX(sbox_args, sbox, sbox_selector, y, x, views, r, &lowmc->mask, &vars, LOWMC_N, shares);    \
    MPC_LOOP_CONST(MUL, x, y, CONCAT(round->l, matrix_postfix), shares);                           \
    MPC_LOOP_CONST_C(XOR, x, x, round->constant, shares, ch);                                      \
    MPC_LOOP_CONST(ADDMUL, x, lowmc_key, CONCAT(round->k, matrix_postfix), shares);                \
  }
#endif

#if defined(M_FIXED_10)
#undef SBOX_ARGS
#undef SBOX_SIGN
#undef SBOX_VERIFY

#define SBOX_ARGS 5
#define SBOX_SIGN mpc_sbox_layer_bitsliced_uint64
#define SBOX_VERIFY mpc_sbox_layer_bitsliced_verify_uint64

#define sbox_selector uint64
#else
#define sbox_selector mzd
#endif

static inline void N_SIGN(lowmc_t const* lowmc, mpc_lowmc_key_t* lowmc_key, mzd_local_t const* p,
                          view_t* views, in_out_shares_t* in_out_shares, rvec_t* rvec) {
  mpc_copy(in_out_shares->s, lowmc_key, SC_PROOF);
  ++in_out_shares;
  CONCAT(VARS, SBOX_ARGS)(SC_PROOF, LOWMC_N);
  mzd_local_t** x = in_out_shares->s;
  mzd_local_t* y[SC_PROOF];
  mzd_local_init_multiple_ex(y, SC_PROOF, 1, (LOWMC_N), false);

  MPC_LOOP_CONST(MUL, x, lowmc_key, CONCAT(lowmc->k0, matrix_postfix), SC_PROOF);
  MPC_LOOP_CONST_C(XOR, x, x, p, SC_PROOF, 0);

  lowmc_round_t const* round = lowmc->rounds;
  loop_impl(SBOX_ARGS, SBOX_SIGN, sbox_selector, 0, SC_PROOF);

  mzd_local_free_multiple(y);
  CONCAT(VARS_FREE, SBOX_ARGS);
}

static inline void N_VERIFY(lowmc_t const* lowmc, mzd_local_t const* p, view_t* views,
                            in_out_shares_t* in_out_shares, rvec_t* rvec, unsigned int ch) {
  mzd_local_t* const* lowmc_key = &in_out_shares->s[0];

  ++in_out_shares;
  CONCAT(VARS, SBOX_ARGS)(SC_VERIFY, LOWMC_N);

  mzd_local_t* x[2 * SC_VERIFY];
  mzd_local_t** y = &x[SC_VERIFY];
  mzd_local_init_multiple_ex(x, 2 * SC_VERIFY, 1, LOWMC_N, false);

  MPC_LOOP_CONST(MUL, x, lowmc_key, CONCAT(lowmc->k0, matrix_postfix), SC_VERIFY);
  MPC_LOOP_CONST_C(XOR, x, x, p, SC_VERIFY, ch);

  lowmc_round_t const* round = lowmc->rounds;
  loop_impl(SBOX_ARGS, SBOX_VERIFY, sbox_selector, ch, SC_VERIFY);
  mpc_copy(in_out_shares->s, x, SC_VERIFY);

  mzd_local_free_multiple(x);
  CONCAT(VARS_FREE, SBOX_ARGS);
}

#if defined(M_FIXED_10)
#undef SBOX_SIGN
#undef SBOX_VERIFY
#undef SBOX_ARGS
#endif

#undef sbox_selector
#undef loop_impl
#undef N_SIGN
#undef N_VERIFY

// vim: ft=c
