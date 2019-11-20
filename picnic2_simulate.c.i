/*! @file picnic2_impl.c
 *  @brief This is the main file of the signature scheme for the Picnic2
 *  parameter sets.
 *
 *  This file is part of the reference implementation of the Picnic signature scheme.
 *  See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */

#if defined(FN_ATTR)
FN_ATTR
#endif
static int SIM_ONLINE(mzd_local_t* maskedKey, shares_t* mask_shares, randomTape_t* tapes, msgs_t* msgs,
                      const mzd_local_t* plaintext, const uint32_t* pubKey,
                      const picnic_instance_t* params) {
  int ret                 = 0;
  mzd_local_t state[((LOWMC_N) + 255) / 256];
  shares_t* key_masks = allocateShares(params->input_size*8); // Make a copy to use when computing each round key
  shares_t* mask2_shares = allocateShares(params->input_size*8);
  uint8_t* unopened_msgs = NULL;

  if (msgs->unopened >= 0) { // We are in verify, save the unopenend parties msgs
    unopened_msgs = malloc(params->view_size + params->input_size);
    memcpy(unopened_msgs, msgs->msgs[msgs->unopened], params->view_size + params->input_size);
  }

  copyShares(key_masks, mask_shares);

  mzd_local_t roundKey[((LOWMC_N) + 255) / 256];
  MPC_MUL(roundKey, maskedKey, LOWMC_INSTANCE.k0_matrix,
          mask_shares);                                       // roundKey = maskedKey * KMatrix[0]
  XOR(state, roundKey, plaintext);

  shares_t* round_key_masks = allocateShares(mask_shares->numWords);
  for (uint32_t r = 0; r < LOWMC_R; r++) {
    copyShares(round_key_masks, key_masks);
    MPC_MUL(roundKey, maskedKey, LOWMC_INSTANCE.rounds[r].k_matrix, round_key_masks);

    mpc_sbox(state, mask_shares, tapes, msgs, unopened_msgs, params);
    MPC_MUL(state, state, LOWMC_INSTANCE.rounds[r].l_matrix,
            mask_shares); // state = state * LMatrix (r-1)
    XOR(state, state, LOWMC_INSTANCE.rounds[r].constant);
    XOR(state, state, roundKey);
    mpc_xor_masks(mask_shares, mask_shares, round_key_masks);
  }
  freeShares(round_key_masks);

  /* Unmask the output, and check that it's correct */
  if (msgs->unopened >= 0) {
    /* During signature verification we have the shares of the output for
     * the unopened party already in msgs, but not in mask_shares. */
    for (size_t i = 0; i < LOWMC_N; i++) {
      uint8_t share = getBit(unopened_msgs, msgs->pos + i);
      setBit((uint8_t*)&mask_shares->shares[i], msgs->unopened, share);
    }
  }
  uint32_t output[params->output_size*8 / 32];
  uint32_t outstate[params->output_size*8 / 32];
  mzd_to_char_array((uint8_t*)outstate, state, params->output_size);
  reconstructShares(output, mask_shares);
  xor_word_array(output, output, outstate, (params->output_size*8 / 32));

  if (memcmp(output, pubKey, params->output_size) != 0) {
#if !defined(NDEBUG)
    printf("%s: output does not match pubKey\n", __func__);
    printf("pubKey: ");
    print_hex(stdout, (uint8_t*)pubKey, params->output_size);
    printf("\noutput: ");
    print_hex(stdout, (uint8_t*)output, params->output_size);
    printf("\n");
#endif
    ret = -1;
    goto Exit;
  }

  broadcast(mask_shares, msgs);
  msgsTranspose(msgs);

  free(unopened_msgs);
  freeShares(key_masks);
  freeShares(mask2_shares);

Exit:
  return ret;
}
