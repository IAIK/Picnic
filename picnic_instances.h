/*
 *  This file is part of the optimized implementation of the Picnic signature scheme.
 *  See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */

#ifndef PICNIC_INSTANCES_H
#define PICNIC_INSTANCES_H

#include "lowmc.h"
#if defined(WITH_ZKBPP)
#include "mpc_lowmc.h"
#endif
#if defined(WITH_KKW)
#include "picnic3_simulate.h"
#endif
#include "picnic.h"

#define SALT_SIZE 32

/* max digest and seed size */
#if defined(WITH_LOWMC_255_255_4) || defined(WITH_LOWMC_256_256_38)
#define MAX_DIGEST_SIZE 64
#define MAX_SEED_SIZE 32
#elif defined(WITH_LOWMC_192_192_4) || defined(WITH_LOWMC_192_192_30)
#define MAX_DIGEST_SIZE 48
#define MAX_SEED_SIZE 24
#elif defined(WITH_LOWMC_129_129_4) || defined(WITH_LOWMC_128_128_20)
#define MAX_DIGEST_SIZE 32
#define MAX_SEED_SIZE 16
#endif

typedef struct picnic_instance_t {
  uint16_t num_rounds;                    // T
  uint8_t digest_size;                    // bytes
  uint8_t seed_size;                      // bytes
  uint8_t input_size;                     // bytes
  uint8_t output_size;                    // bytes
  uint8_t view_size;                      // bytes
  uint8_t num_opened_rounds;              // u (KKW only)
  uint8_t num_MPC_parties;                // N (KKW only)
  uint8_t unruh_without_input_bytes_size; // bytes (Unruh only)

  lowmc_parameters_t lowmc;
  lowmc_implementation_f impl_lowmc;
#if defined(WITH_ZKBPP)
  lowmc_store_implementation_f impl_lowmc_store;
  zkbpp_lowmc_implementation_f impl_zkbpp_lowmc;
  zkbpp_lowmc_verify_implementation_f impl_zkbpp_lowmc_verify;
  zkbpp_share_implementation_f impl_mzd_share;
#endif
#if defined(WITH_KKW)
  lowmc_compute_aux_implementation_f impl_lowmc_aux;
  lowmc_simulate_online_f impl_lowmc_simulate_online;
#endif
} picnic_instance_t;

const picnic_instance_t* picnic_instance_get(picnic_params_t param);

PICNIC_EXPORT size_t PICNIC_CALLING_CONVENTION picnic_get_lowmc_block_size(picnic_params_t param);

/* Prefix values for domain separation */
static const uint8_t HASH_PREFIX_0 = 0;
static const uint8_t HASH_PREFIX_1 = 1;
static const uint8_t HASH_PREFIX_2 = 2;
static const uint8_t HASH_PREFIX_3 = 3;
static const uint8_t HASH_PREFIX_4 = 4;
static const uint8_t HASH_PREFIX_5 = 5;

#endif
