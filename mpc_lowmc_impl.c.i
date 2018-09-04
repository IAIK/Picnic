/*
 *  This file is part of the optimized implementation of the Picnic signature scheme.
 *  See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */

#if defined(M_FIXED_10)
#undef SBOX_ARGS
#undef SBOX_SIGN
#undef SBOX_VERIFY

#define LOWMC_M 10
#define SBOX_ARGS 5
#define SBOX_SIGN mpc_sbox_layer_bitsliced_uint64_10
#define SBOX_VERIFY mpc_sbox_layer_bitsliced_verify_uint64_10

#define RANDTAPE R_uint64
#define SBOX SBOX_uint64
#elif defined(M_FIXED_1)
#define LOWMC_M 1
#undef SBOX_ARGS
#undef SBOX_SIGN
#undef SBOX_VERIFY

#define SBOX_ARGS 5
#define SBOX_SIGN mpc_sbox_layer_bitsliced_uint64_1
#define SBOX_VERIFY mpc_sbox_layer_bitsliced_verify_uint64_1

#define RANDTAPE R_uint64
#define SBOX SBOX_uint64
#else
#define LOWMC_M (lowmc->m)
#define RANDTAPE R_mzd
#define SBOX SBOX_mzd
#endif

#if defined(LOWMC_INSTANCE)
#define lowmc LOWMC_INSTANCE
#else
#define lowmc lowmc_instance
#endif

static void N_SIGN(lowmc_t const* lowmc_instance, mpc_lowmc_key_t* lowmc_key, mzd_local_t const* p,
                   view_t* views, in_out_shares_t* in_out_shares, rvec_t* rvec,
                   recorded_state_t* recorded_state) {
#if defined(LOWMC_INSTANCE)
  (void)lowmc_instance;
#endif
  mpc_copy(in_out_shares->s, lowmc_key, SC_PROOF);
  ++in_out_shares;

  CONCAT(VARS, SBOX_ARGS)(SC_PROOF, LOWMC_N);
  mzd_local_t** x = in_out_shares->s;
  mzd_local_t* y[SC_PROOF];
  mzd_local_init_multiple_ex(y, SC_PROOF, 1, (LOWMC_N), false);

#define reduced_shares (SC_PROOF - 1)

  MPC_LOOP_CONST(MUL, x, lowmc_key, CONCAT(lowmc->k0, matrix_postfix), reduced_shares);
  MPC_LOOP_CONST_C(XOR, x, x, p, reduced_shares, 0);

#define RECOVER_FROM_STATE(x, i)                                                                   \
  XOR((x)[SC_PROOF - 1], (x)[0], (x)[1]);                                                          \
  XOR((x)[SC_PROOF - 1], (x)[SC_PROOF - 1], recorded_state->state[i])
#define ch 0
#define shares SC_PROOF
#define sbox SBOX_SIGN
#include "mpc_lowmc_loop.c.i"
#undef reduced_shares
#undef RECOVER_FROM_STATE
#undef ch
#undef shares
#undef sbox

  mzd_local_free_multiple(y);
  CONCAT(VARS_FREE, SBOX_ARGS);
}

static void N_VERIFY(lowmc_t const* lowmc_instance, mzd_local_t const* p, view_t* views,
                     in_out_shares_t* in_out_shares, rvec_t* rvec, unsigned int ch) {
#if defined(LOWMC_INSTANCE)
  (void)lowmc_instance;
#endif
  mzd_local_t* const* lowmc_key = &in_out_shares->s[0];
  ++in_out_shares;

  CONCAT(VARS, SBOX_ARGS)(SC_VERIFY, LOWMC_N);
  mzd_local_t* x[2 * SC_VERIFY];
  mzd_local_t** y = &x[SC_VERIFY];
  mzd_local_init_multiple_ex(x, 2 * SC_VERIFY, 1, LOWMC_N, false);

  MPC_LOOP_CONST(MUL, x, lowmc_key, CONCAT(lowmc->k0, matrix_postfix), SC_VERIFY);
  MPC_LOOP_CONST_C(XOR, x, x, p, SC_VERIFY, ch);

#define shares SC_VERIFY
#define reduced_shares shares
#define sbox SBOX_VERIFY
#include "mpc_lowmc_loop.c.i"
#undef sbox
#undef reduced_shares
#undef shares

  mpc_copy(in_out_shares->s, x, SC_VERIFY);

  mzd_local_free_multiple(x);
  CONCAT(VARS_FREE, SBOX_ARGS);
}

#if defined(M_FIXED_10) || defined(M_FIXED_1)
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
#undef LOWMC_M
#undef lowmc

// vim: ft=c
