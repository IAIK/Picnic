/*
 *  This file is part of the optimized implementation of the Picnic signature scheme.
 *  See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "lowmc_pars.h"

#include "macros.h"
#include "mzd_additional.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(WITH_CUSTOM_INSTANCES)
static bool prepare_masks(mask_t* mask, unsigned int n, unsigned int m) {
  mask->x0   = mzd_local_init(1, n);
  mask->x1   = mzd_local_init_ex(1, n, false);
  mask->x2   = mzd_local_init_ex(1, n, false);
  mask->mask = mzd_local_init(1, n);
  if (!mask->x0 || !mask->x1 || !mask->x2 || !mask->mask) {
    return false;
  }

  const unsigned int bound = n - 3 * m;
  // 1s for linear part
  for (unsigned int i = 0; i < bound; ++i) {
    mzd_local_write_bit(mask->mask, 0, i, 1);
  }
  // 1s for a
  for (unsigned int i = bound; i < n; i += 3) {
    mzd_local_write_bit(mask->x0, 0, i, 1);
  }
  // 1s for b
  mzd_shift_left(mask->x1, mask->x0, 1);
  // 1s for c
  mzd_shift_left(mask->x2, mask->x0, 2);
  return true;
}
#endif

#if defined(MUL_M4RI)
bool lowmc_init(lowmc_t* lowmc) {
  if (!lowmc) {
    return false;
  }

  if (lowmc->n - 3 * lowmc->m < 2 || lowmc->n != lowmc->k) {
    return false;
  }

  lowmc->k0_lookup = mzd_precompute_matrix_lookup(lowmc->k0_matrix);
#if defined(REDUCED_LINEAR_LAYER)
  lowmc->precomputed_non_linear_part_lookup =
      mzd_precompute_matrix_lookup(lowmc->precomputed_non_linear_part_matrix);
#endif
  for (unsigned int i = 0; i < lowmc->r; ++i) {
    lowmc->rounds[i].l_lookup = mzd_precompute_matrix_lookup(lowmc->rounds[i].l_matrix);
#if !defined(REDUCED_LINEAR_LAYER)
    lowmc->rounds[i].k_lookup = mzd_precompute_matrix_lookup(lowmc->rounds[i].k_matrix);
#endif
  }

  return true;
}
#endif

#if defined(WITH_CUSTOM_INSTANCES)
static mzd_local_t* read_mzd_t(FILE* file) {
  uint32_t nrows = 0;
  uint32_t ncols = 0;

  size_t ret = fread(&(nrows), sizeof(uint32_t), 1, file);
  ret += fread(&(ncols), sizeof(uint32_t), 1, file);
  if (ret != 2) {
    return NULL;
  }
  ret = 0;

  mzd_local_t* A = mzd_local_init_ex(nrows, ncols, false);
  if (!A) {
    return NULL;
  }

  for (unsigned int i = 0; i < A->nrows; i++) {
    ret += fread(ROW(A, i), A->rowstride * sizeof(word), 1, file);
  }

  if (ret != A->nrows) {
    mzd_local_free(A);
    return NULL;
  }

  return A;
}

const uint32_t supported_instance_type =
#if defined(REDUCED_LINEAR_LAYER)
  LOWMC_INSTANCE_RLL |
#endif
  0;

bool lowmc_read_file(lowmc_t* lowmc, const char* file_name) {
  if (!lowmc) {
    return false;
  }

  FILE* file = fopen(file_name, "r+");
  if (!file) {
    return false;
  }

  uint32_t instance_type = 0;
  size_t ret = fread(&instance_type, sizeof(instance_type), 1, file);
  if (ret != 1 || instance_type != supported_instance_type) {
    fclose(file);
    return false;
  }

  ret += fread(&lowmc->n, sizeof(lowmc->n), 1, file);
  ret += fread(&lowmc->k, sizeof(lowmc->k), 1, file);
  ret += fread(&lowmc->m, sizeof(lowmc->m), 1, file);
  ret += fread(&lowmc->r, sizeof(lowmc->r), 1, file);
  if (ret != 4 || lowmc->n != lowmc->k || lowmc->n < 3 * lowmc->m) {
    fclose(file);
    return false;
  }

  lowmc->needs_free = true;
  lowmc->k0_matrix  = read_mzd_t(file);
  if (!lowmc->k0_matrix) {
    goto error;
  }

  lowmc->rounds     = calloc(lowmc->r, sizeof(lowmc_round_t));
  if (!lowmc->rounds) {
    goto error;
  }

  for (uint32_t i = 0; i < lowmc->r; ++i) {
    lowmc->rounds[i].l_matrix = read_mzd_t(file);
#if !defined(REDUCED_LINEAR_LAYER)
    lowmc->rounds[i].k_matrix = read_mzd_t(file);
    lowmc->rounds[i].constant = read_mzd_t(file);
#endif
  }
#if defined(REDUCED_LINEAR_LAYER)
  lowmc->precomputed_non_linear_part_matrix = read_mzd_t(file);
  lowmc->precomputed_constant_linear        = read_mzd_t(file);
  lowmc->precomputed_constant_non_linear    = read_mzd_t(file);
#endif

  if (lowmc->m != 10 && !prepare_masks(&lowmc->mask, lowmc->n, lowmc->m)) {
    goto error;
  }

  fclose(file);
  return true;

error:
  fclose(file);
  lowmc_clear(lowmc);
  return false;
}
#endif

void lowmc_clear(lowmc_t* lowmc) {
  for (unsigned int i = 0; i < lowmc->r; ++i) {
#if defined(MUL_M4RI)
#if !defined(REDUCED_LINEAR_LAYER)
    mzd_local_free(lowmc->rounds[i].k_lookup);
#endif
    mzd_local_free(lowmc->rounds[i].l_lookup);
#endif
#if defined(WITH_CUSTOM_INSTANCES)
    if (lowmc->needs_free) {
#if !defined(REDUCED_LINEAR_LAYER)
      mzd_local_free((mzd_local_t*)lowmc->rounds[i].constant);
#endif
      mzd_local_free((mzd_local_t*)lowmc->rounds[i].l_matrix);
#if !defined(REDUCED_LINEAR_LAYER)
      mzd_local_free((mzd_local_t*)lowmc->rounds[i].k_matrix);
#endif
    }
#endif
  }
#if defined(REDUCED_LINEAR_LAYER) && defined(WITH_CUSTOM_INSTANCES)
  if (lowmc->needs_free) {
    mzd_local_free((mzd_local_t*)lowmc->precomputed_constant_non_linear);
    mzd_local_free((mzd_local_t*)lowmc->precomputed_constant_linear);
    mzd_local_free((mzd_local_t*)lowmc->precomputed_non_linear_part_matrix);
  }
#endif
#if defined(MUL_M4RI)
  mzd_local_free(lowmc->k0_lookup);
#if defined(REDUCED_LINEAR_LAYER)
  mzd_local_free(lowmc->precomputed_non_linear_part_lookup);
#endif
#endif
#if defined(WITH_CUSTOM_INSTANCES)
  if (lowmc->needs_free) {
    mzd_local_free((mzd_local_t*)lowmc->k0_matrix);
  }
#endif

#if defined(WITH_CUSTOM_INSTANCES)
  mzd_local_free(lowmc->mask.x0);
  mzd_local_free(lowmc->mask.x1);
  mzd_local_free(lowmc->mask.x2);
  mzd_local_free(lowmc->mask.mask);
  lowmc->needs_free = false;
#endif
}
