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
static int SIM_ONLINE(mzd_local_t** maskedKey, shares_t* mask_shares, randomTape_t* tapes,
                      msgs_t* msgs, const mzd_local_t* plaintext, const uint8_t* pubKey,
                      const picnic_instance_t* params) {
  int ret = 0;
  mzd_local_t* state[PACKING_FACTOR];
  mzd_local_init_multiple_ex(state, PACKING_FACTOR, 1, LOWMC_N, false);
  mzd_local_t* temp[PACKING_FACTOR];
  mzd_local_init_multiple_ex(temp, PACKING_FACTOR, 1, LOWMC_N, false);
//  mzd_local_t state_a[PACKING_FACTOR][((LOWMC_N) + 255) / 256];
//  mzd_local_t temp_a[PACKING_FACTOR][((LOWMC_N) + 255) / 256];
//  mzd_local_t* state[PACKING_FACTOR];
//  mzd_local_t* temp[PACKING_FACTOR];
//  for (uint32_t k = 0; k < PACKING_FACTOR; k++) {
//	state[k] = state_a[k];
//	temp[k] = temp_a[k];
//  }
  uint8_t* unopened_msgs[PACKING_FACTOR];

  if (msgs->unopened != NULL) { // We are in verify, save the unopenend parties msgs
    for (uint32_t k = 0; k < PACKING_FACTOR; k++) {
      unopened_msgs[k] = malloc(params->view_size);
      memcpy(unopened_msgs[k], msgs->msgs[(64 / PACKING_FACTOR) * k + msgs->unopened[k]],
             params->view_size);
    }
  }

//  MPC_MUL(temp, maskedKey, LOWMC_INSTANCE.k0_matrix,
//          mask_shares); // roundKey = maskedKey * KMatrix[0]
  for (uint32_t k = 0; k < PACKING_FACTOR; k++) {
    MUL(temp[k], maskedKey[k], LOWMC_INSTANCE.k0_matrix);
    XOR(state[k], temp[k], plaintext);
  }

  for (uint32_t r = 0; r < LOWMC_R; r++) {
    // MPC_MUL(roundKey, maskedKey, LOWMC_INSTANCE.rounds[r].k_matrix, round_key_masks);
    for (size_t i = 0; i < LOWMC_N; i++) {
      mask_shares->shares[i] = tapesToWord(tapes);
    }
    mpc_sbox(state, mask_shares, tapes, msgs, unopened_msgs, params);
    // MPC_MUL(state, state, LOWMC_INSTANCE.rounds[r].l_matrix,
    //        mask_shares); // state = state * LMatrix (r-1)
    for (uint32_t k = 0; k < PACKING_FACTOR; k++) {
      MUL(temp[k], state[k], LOWMC_INSTANCE.rounds[r].l_matrix);
      XOR(state[k], temp[k], LOWMC_INSTANCE.rounds[r].constant);
      ADDMUL(state[k], maskedKey[k], LOWMC_INSTANCE.rounds[r].k_matrix);
    }
  }

  for (size_t i = 0; i < LOWMC_N; i++) {
    mask_shares->shares[i] = tapesToWord(tapes);
  }

  /* check that the output is correct */
  for (uint32_t k = 0; k < PACKING_FACTOR; k++) {
    uint8_t output[params->output_size];
    mzd_to_char_array(output, state[k], params->output_size);

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
  }

  msgsTranspose(msgs);

Exit:
  if (msgs->unopened != NULL) {
    for (uint32_t k = 0; k < PACKING_FACTOR; k++) {
      free(unopened_msgs[k]);
    }
  }
  mzd_local_free_multiple(temp);
  mzd_local_free_multiple(state);
  return ret;
}
#endif
