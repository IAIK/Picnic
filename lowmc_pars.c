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

#if defined(WITH_LOWMC_128_128_20)
#include "lowmc_128_128_20.h"
#endif
#if defined(WITH_LOWMC_192_192_30)
#include "lowmc_192_192_30.h"
#endif
#if defined(WITH_LOWMC_256_256_38)
#include "lowmc_256_256_38.h"
#endif

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

bool lowmc_init(lowmc_t* lowmc, unsigned int m, unsigned int n, unsigned int r, unsigned int k) {
  if (!lowmc) {
    return false;
  }

  if (n - 3 * m < 2) {
    return false;
  }

  lowmc->m = m;
  lowmc->n = n;
  lowmc->r = r;
  lowmc->k = k;
#if defined(WITH_CUSTOM_INSTANCES)
  lowmc->needs_free = false;
#endif

  lowmc->rounds = calloc(sizeof(lowmc_round_t), r);

#define LOAD_OPT(N, K, R)                                                                          \
  lowmc->precomputed_non_linear_part_matrix =                                                      \
      lowmc_##N##_##K##_##R##_get_precomputed_round_key_matrix_non_linear_part();                  \
  lowmc->k0_matrix = lowmc_##N##_##K##_##R##_get_precomputed_round_key_matrix_linear_part();       \
  lowmc->precomputed_constant_linear =                                                             \
      lowmc_##N##_##K##_##R##_get_precomputed_constant_linear_part();                              \
  lowmc->precomputed_constant_non_linear =                                                         \
      lowmc_##N##_##K##_##R##_get_precomputed_constant_non_linear_part()

#define LOAD(N, K, R)                                                                              \
  lowmc->k0_matrix = lowmc_##N##_##K##_##R##_get_round_key(0);                                     \
  for (unsigned int i = 0; i < (R); ++i) {                                                         \
    lowmc->rounds[i].k_matrix = lowmc_##N##_##K##_##R##_get_round_key(i + 1);                      \
    lowmc->rounds[i].constant = lowmc_##N##_##K##_##R##_get_round_const(i);                        \
  }

#define LOAD_FROM_FIXED_IMPL(N, K, R, PREC)                                                        \
  for (unsigned int i = 0; i < (R); ++i) {                                                         \
    lowmc->rounds[i].l_matrix = lowmc_##N##_##K##_##R##_get_linear_layer(i);                       \
  }                                                                                                \
  LOAD##PREC(N, K, R);

#if defined(REDUCED_LINEAR_LAYER)
#define LOAD_FROM_FIXED(N, K, R) LOAD_FROM_FIXED_IMPL(N, K, R, _OPT)
#else
#define LOAD_FROM_FIXED(N, K, R) LOAD_FROM_FIXED_IMPL(N, K, R, )
#endif

#if defined(WITH_LOWMC_128_128_20)
  if (n == 128 && k == 128 && r == 20) {
    LOAD_FROM_FIXED(128, 128, 20);
    goto precomp;
  }
#endif
#if defined(WITH_LOWMC_192_192_30)
  if (n == 192 && k == 192 && r == 30) {
    LOAD_FROM_FIXED(192, 192, 30);
    goto precomp;
  }
#endif
#if defined(WITH_LOWMC_256_256_38)
  if (n == 256 && k == 256 && r == 38) {
    LOAD_FROM_FIXED(256, 256, 38);
    goto precomp;
  }
#endif

  lowmc_clear(lowmc);
  return false;

precomp:
#if defined(MUL_M4RI)
  lowmc->k0_lookup = mzd_precompute_matrix_lookup(lowmc->k0_matrix);
#if defined(REDUCED_LINEAR_LAYER)
  lowmc->precomputed_non_linear_part_lookup =
      mzd_precompute_matrix_lookup(lowmc->precomputed_non_linear_part_matrix);
#endif
  for (unsigned int i = 0; i < r; ++i) {
    lowmc->rounds[i].l_lookup = mzd_precompute_matrix_lookup(lowmc->rounds[i].l_matrix);
#if !defined(REDUCED_LINEAR_LAYER)
    lowmc->rounds[i].k_lookup = mzd_precompute_matrix_lookup(lowmc->rounds[i].k_matrix);
#endif
  }
#endif

#if defined(WITH_CUSTOM_INSTANCES)
  prepare_masks(&lowmc->mask, n, m);
#endif
  return true;
}

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
#if !defineD(REDUCED_LINEAR_LAYER)
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
  free(lowmc->rounds);

#if defined(WITH_CUSTOM_INSTANCES)
  mzd_local_free(lowmc->mask.x0);
  mzd_local_free(lowmc->mask.x1);
  mzd_local_free(lowmc->mask.x2);
  mzd_local_free(lowmc->mask.mask);
#endif
}
