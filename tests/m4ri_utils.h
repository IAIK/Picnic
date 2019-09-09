/*
 *  This file is part of the optimized implementation of the Picnic signature scheme.
 *  See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */
#ifndef PICNIC_TEST_M4RI_UTILS_H
#define PICNIC_TEST_M4RI_UTILS_H

#include <m4ri/m4ri.h>

static inline mzd_local_t* mzd_convert_128(const mzd_t* v) {
  mzd_local_t* r = mzd_local_init(v->nrows, 128);

  for (rci_t row = 0; row < v->nrows; ++row) {
    memcpy(&BLOCK(r, row >> 1)->w64[(row & 0x1) << 1], v->rows[row], 2 * sizeof(word));
  }
  return r;
}

static inline mzd_local_t* mzd_convert(const mzd_t* v) {
  if (v->ncols == 128) {
    return mzd_convert_128(v);
  }

  const unsigned int num_uints = (v->ncols + 63) / 64;
  const unsigned int num_blocks = (num_uints + 3) / 4;

  mzd_local_t* r = mzd_local_init(v->nrows, v->ncols);

  for (rci_t i = 0; i < v->nrows; ++i) {
    for (unsigned int j = 0; j < num_blocks; ++j) {
      const unsigned int s = ((j == num_blocks - 1) && (num_uints % 4)) ? (num_uints % 4) : 4;
      memcpy(BLOCK(r, j + i * num_blocks)->w64, &v->rows[i][j * 4], s * sizeof(word));
    }
  }

  return r;
}

#endif
