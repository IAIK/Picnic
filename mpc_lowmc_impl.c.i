#if defined(M_FIXED_10)
#undef SBOX_ARGS
#undef SBOX_SIGN
#undef SBOX_VERIFY

#define SBOX_ARGS 5
#define SBOX_SIGN mpc_sbox_layer_bitsliced_uint64
#define SBOX_VERIFY mpc_sbox_layer_bitsliced_verify_uint64

#define RANDTAPE R_uint64
#define SBOX SBOX_uint64
#else
#define RANDTAPE R_mzd
#define SBOX SBOX_mzd
#endif

#if defined(LOWMC_INSTANCE)
#define lowmc LOWMC_INSTANCE
#else
#define lowmc lowmc_instance
#endif

static void N_SIGN(lowmc_t const* lowmc_instance, mpc_lowmc_key_t* lowmc_key, mzd_local_t const* p,
                          view_t* views, in_out_shares_t* in_out_shares, rvec_t* rvec) {
  (void) lowmc_instance;

  mpc_copy(in_out_shares->s, lowmc_key, SC_PROOF);
  ++in_out_shares;

  CONCAT(VARS, SBOX_ARGS)(SC_PROOF, LOWMC_N);
  mzd_local_t** x = in_out_shares->s;
  mzd_local_t* y[SC_PROOF];
  mzd_local_init_multiple_ex(y, SC_PROOF, 1, (LOWMC_N), false);

  MPC_LOOP_CONST(MUL, x, lowmc_key, CONCAT(lowmc->k0, matrix_postfix), SC_PROOF);
  MPC_LOOP_CONST_C(XOR, x, x, p, SC_PROOF, 0);

#define ch 0
#define shares SC_PROOF
#define sbox SBOX_SIGN
#include "mpc_lowmc_loop.c.i"
#undef ch
#undef shares
#undef sbox

  mzd_local_free_multiple(y);
  CONCAT(VARS_FREE, SBOX_ARGS);
}

static void N_VERIFY(lowmc_t const* lowmc_instance, mzd_local_t const* p, view_t* views,
                            in_out_shares_t* in_out_shares, rvec_t* rvec, unsigned int ch) {
  (void) lowmc_instance;

  mzd_local_t* const* lowmc_key = &in_out_shares->s[0];
  ++in_out_shares;

  CONCAT(VARS, SBOX_ARGS)(SC_VERIFY, LOWMC_N);
  mzd_local_t* x[2 * SC_VERIFY];
  mzd_local_t** y = &x[SC_VERIFY];
  mzd_local_init_multiple_ex(x, 2 * SC_VERIFY, 1, LOWMC_N, false);

  MPC_LOOP_CONST(MUL, x, lowmc_key, CONCAT(lowmc->k0, matrix_postfix), SC_VERIFY);
  MPC_LOOP_CONST_C(XOR, x, x, p, SC_VERIFY, ch);

#define shares SC_VERIFY
#define sbox SBOX_VERIFY
#include "mpc_lowmc_loop.c.i"
#undef shares
#undef sbox

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
#undef RANDTAPE
#undef SBOX
#undef lowmc

// vim: ft=c
