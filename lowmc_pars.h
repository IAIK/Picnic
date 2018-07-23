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

typedef struct {
  mzd_local_t* x0;
  mzd_local_t* x1;
  mzd_local_t* x2;
  mzd_local_t* mask;
} mask_t;

/**
 * Masks for 10 S-boxes.
 */
#define MASK_X0I UINT64_C(0x2492492400000000)
#define MASK_X1I UINT64_C(0x4924924800000000)
#define MASK_X2I UINT64_C(0x9249249000000000)
#define MASK_MASK UINT64_C(0x00000003ffffffff)

/**
 * LowMC instances
 */
#define LOWMC_L1_N 128
#define LOWMC_L1_M 10
#define LOWMC_L1_K LOWMC_L1_N
#define LOWMC_L1_R 20

#define LOWMC_L3_N 192
#define LOWMC_L3_M 10
#define LOWMC_L3_K LOWMC_L3_N
#define LOWMC_L3_R 30

#define LOWMC_L5_N 256
#define LOWMC_L5_M 10
#define LOWMC_L5_K LOWMC_L5_N
#define LOWMC_L5_R 38

typedef struct {
#if !defined(REDUCED_LINEAR_LAYER)
  const mzd_local_t* k_matrix;
#endif
  const mzd_local_t* l_matrix;
#if !defined(REDUCED_LINEAR_LAYER)
  const mzd_local_t* constant;
#endif

#if defined(MUL_M4RI)
#if !defined(REDUCED_LINEAR_LAYER)
  mzd_local_t* k_lookup;
#endif
  mzd_local_t* l_lookup;
#endif
} lowmc_round_t;

/**
 * LowMC definition
 */
typedef struct {
  unsigned int m;
  unsigned int n;
  unsigned int r;
  unsigned int k;

  const mzd_local_t* k0_matrix; // K_0 or K_0 + precomputed if reduced_linear_layer is set
#if defined(MUL_M4RI)
  mzd_local_t* k0_lookup;
  lowmc_round_t* rounds;
#else
  const lowmc_round_t* rounds;
#endif

#if defined(REDUCED_LINEAR_LAYER)
  const mzd_local_t* precomputed_non_linear_part_matrix;
#if defined(MUL_M4RI)
  mzd_local_t* precomputed_non_linear_part_lookup;
#endif
  const mzd_local_t* precomputed_constant_linear;
  const mzd_local_t* precomputed_constant_non_linear;
#endif

#if defined(WITH_CUSTOM_INSTANCES)
  mask_t mask;
  bool needs_free;
#endif
} lowmc_t;

#if defined(MUL_M4RI)
/**
 * Initiaizes lookup tables of a LowMC instance
 *
 * \return parameters defining a LowMC instance
 */
bool lowmc_init(lowmc_t* lowmc);
#endif

/**
 * Clears the allocated LowMC parameters
 *
 * \param lowmc the LowMC parameters to be cleared
 */
void lowmc_clear(lowmc_t* lowmc);

bool lowmc_read_file(lowmc_t* lowmc, unsigned int m, unsigned int n, unsigned int r,
                     unsigned int k);

#endif
