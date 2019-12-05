/*
 *  This file is part of the optimized implementation of the Picnic signature scheme.
 *  See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */

#ifndef LOWMC_PARS_H
#define LOWMC_PARS_H

#include <stddef.h>

#include "mzd_additional.h"

typedef mzd_local_t lowmc_key_t;

#define MAX_LOWMC_BLOCK_SIZE 32
#define MAX_LOWMC_BLOCK_SIZE_BITS (MAX_LOWMC_BLOCK_SIZE * 8)
#define MAX_LOWMC_KEY_SIZE MAX_LOWMC_BLOCK_SIZE
#define MAX_LOWMC_KEY_SIZE_BITS (MAX_LOWMC_KEY_SIZE * 8)
#define MAX_LOWMC_ROUNDS 5

/**
 * LowMC instances
 */
#define LOWMC_L1_N 126
#define LOWMC_L1_M 42
#define LOWMC_L1_K LOWMC_L1_N
#define LOWMC_L1_R 4

#define LOWMC_L3_N 192
#define LOWMC_L3_M 64
#define LOWMC_L3_K LOWMC_L3_N
#define LOWMC_L3_R 4

#define LOWMC_L5_N 255
#define LOWMC_L5_M 85
#define LOWMC_L5_K LOWMC_L5_N
#define LOWMC_L5_R 4

typedef struct {
  const mzd_local_t* k_matrix;
  const mzd_local_t* l_matrix;
  const mzd_local_t* li_matrix;
  const mzd_local_t* constant;
} lowmc_round_t;

/**
 * LowMC definition
 */
typedef struct {
  uint32_t m;
  uint32_t n;
  uint32_t r;
  uint32_t k;

  const mzd_local_t* k0_matrix; // K_0 or K_0 + precomputed if reduced_linear_layer is set
  const mzd_local_t* ki0_matrix; // inverse of K_0 or K_0 + precomputed if reduced_linear_layer is set
  const lowmc_round_t* rounds;

} lowmc_t;

#endif
