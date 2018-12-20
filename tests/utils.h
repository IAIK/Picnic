/*
 *  This file is part of the optimized implementation of the Picnic signature scheme.
 *  See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */

#include <string.h>

static inline bool mzd_local_equal(mzd_local_t const* first, mzd_local_t const* second, unsigned int rows, unsigned cols) {
  const unsigned int num_uints = (cols + 63) / 64;
  const unsigned int num_blocks = (num_uints + 3) / 4;

  for (unsigned int i = 0; i < rows; ++i) {
    for (unsigned int j = 0; j < num_blocks; ++j) {
      const unsigned int s = ((j == num_blocks - 1) && (num_uints % 4)) ? (num_uints % 4) : 4;
      if (memcmp(CONST_BLOCK(first, j + i * num_blocks)->w64,
                 CONST_BLOCK(second, j + i * num_blocks)->w64, s * sizeof(word)) != 0) {
        return false;
      }
    }
  }

  return true;
}

#if defined(WITH_M4RI_UTILS)
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
