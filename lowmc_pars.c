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
static void prepare_masks(mask_t* mask, unsigned int n, unsigned int m) {
  mask->x0   = mzd_local_init(1, n);
  mask->x1   = mzd_local_init_ex(1, n, false);
  mask->x2   = mzd_local_init_ex(1, n, false);
  mask->mask = mzd_local_init(1, n);

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
static mzd_local_t* readMZD_TStructFromFile(FILE* file) {
  int ret   = 0;
  int nrows = 0;
  int ncols = 0;
  ret += fread(&(nrows), sizeof(uint32_t), 1, file);
  ret += fread(&(ncols), sizeof(uint32_t), 1, file);

  mzd_local_t* A = mzd_local_init_ex(nrows, ncols, false);
  for (unsigned int i = 0; i < A->nrows; i++) {
    ret += fread(ROW(A, i), A->rowstride * sizeof(word), 1, file);
  }

  return A;
}

bool lowmc_read_file(lowmc_t* lowmc, unsigned int m, unsigned int n, unsigned int r,
                     unsigned int k) {
  if (!lowmc) {
    return false;
  }

  char file_name[4 * 11] = {0};
  if (snprintf(file_name, sizeof(file_name), "%u-%u-%u-%u", m, n, r, k) < 0) {
    return false;
  }

  FILE* file = fopen(file_name, "r+");
  if (file) {
    int ret = 0;
    ret     = fread(&lowmc->m, sizeof(lowmc->m), 1, file);
    ret += fread(&lowmc->n, sizeof(lowmc->n), 1, file);
    ret += fread(&lowmc->r, sizeof(lowmc->r), 1, file);
    ret += fread(&lowmc->k, sizeof(lowmc->k), 1, file);

    if (lowmc->m != m || lowmc->n != n || lowmc->r != r || lowmc->k != k) {
      fclose(file);
      return false;
    }

    lowmc->needs_free = true;
    lowmc->k0_matrix  = readMZD_TStructFromFile(file);
    lowmc->rounds     = calloc(r, sizeof(lowmc_round_t));
    for (unsigned int i = 0; i < lowmc->r; ++i) {
#if !defined(REDUCED_LINEAR_LAYER)
      lowmc->rounds[i].k_matrix = readMZD_TStructFromFile(file);
#endif
      lowmc->rounds[i].l_matrix = readMZD_TStructFromFile(file);
#if !defined(REDUCED_LINEAR_LAYER)
      lowmc->rounds[i].constant = readMZD_TStructFromFile(file);
#endif
    }
#if defined(REDUCED_LINEAR_LAYER)
    lowmc->precomputed_non_linear_part_matrix = readMZD_TStructFromFile(file);
    lowmc->precomputed_constant_linear        = readMZD_TStructFromFile(file);
    lowmc->precomputed_constant_non_linear    = readMZD_TStructFromFile(file);
#endif

    fclose(file);
  }

  prepare_masks(&lowmc->mask, n, m);
  return true;
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
#endif
}
