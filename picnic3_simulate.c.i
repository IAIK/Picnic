/*! @file picnic3_impl.c
 *  @brief This is the main file of the signature scheme for the Picnic3
 *  parameter sets.
 *
 *  This file is part of the reference implementation of the Picnic signature scheme.
 *  See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */

#if defined(LOWMC_INSTANCE)
#if defined(FN_ATTR)
FN_ATTR
#endif
static int SIM_ONLINE(mzd_local_t* maskedKey, shares_t* mask_shares, randomTape_t* tapes,
                      msgs_t* msgs, const mzd_local_t* plaintext, const uint8_t* pubKey,
                      const picnic_instance_t* params) {
  int ret                = 0;
  mzd_local_t* state     = mzd_local_init_ex(1, LOWMC_N, false);
  mzd_local_t* temp      = mzd_local_init_ex(1, LOWMC_N, false);
  uint8_t* unopened_msgs = NULL;

  if (msgs->unopened != -1) { // We are in verify, save the unopenend parties msgs
    unopened_msgs = malloc(params->view_size);
    memcpy(unopened_msgs, msgs->msgs[msgs->unopened], params->view_size);
  }

  //  MPC_MUL(temp, maskedKey, LOWMC_INSTANCE.k0_matrix,
  //          mask_shares); // roundKey = maskedKey * KMatrix[0]
  MUL(temp, maskedKey, LOWMC_INSTANCE.k0_matrix);
  XOR(state, temp, plaintext);

  for (uint32_t r = 0; r < LOWMC_R; r++) {
    // MPC_MUL(roundKey, maskedKey, LOWMC_INSTANCE.rounds[r].k_matrix, round_key_masks);
    for (size_t i = 0; i < LOWMC_N; i++) {
      mask_shares->shares[i] = tapesToWord(tapes);
    }
    mpc_sbox(state, mask_shares, tapes, msgs, unopened_msgs, params);
    // MPC_MUL(state, state, LOWMC_INSTANCE.rounds[r].l_matrix,
    //        mask_shares); // state = state * LMatrix (r-1)
    MUL(temp, state, LOWMC_INSTANCE.rounds[r].l_matrix);
    XOR(state, temp, LOWMC_INSTANCE.rounds[r].constant);
    ADDMUL(state, maskedKey, LOWMC_INSTANCE.rounds[r].k_matrix);
  }

  /* check that the output is correct */
  uint8_t output[params->output_size];
  mzd_to_char_array(output, state, params->output_size);

  if (memcmp(output, pubKey, params->output_size) != 0) {
#if !defined(NDEBUG)
    printf("%s: output does not match pubKey\n", __func__);
    printf("pubKey: ");
    print_hex(stdout, pubKey, params->output_size);
    printf("\noutput: ");
    print_hex(stdout, output, params->output_size);
    printf("\n");
#endif
    ret = -1;
    goto Exit;
  }

  msgsTranspose(msgs);

Exit:
  if (msgs->unopened != -1) {
    free(unopened_msgs);
  }
  mzd_local_free(temp);
  mzd_local_free(state);
  return ret;
}
#endif
