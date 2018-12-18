/*
 *  This file is part of the optimized implementation of the Picnic signature scheme.
 *  See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */

#include <openssl/rand.h>
#include <string.h>

static inline void mzd_randomize_ssl(mzd_local_t* val) {
  for (unsigned int i = 0; i < val->nrows; ++i) {
    RAND_bytes((unsigned char*)ROW(val, i), val->width * sizeof(word));
  }
}

static inline bool mzd_local_equal(mzd_local_t const* first, mzd_local_t const* second) {
  const unsigned int rows  = first->nrows;
  const unsigned int width = first->width;

  for (unsigned int r = 0; r < rows; ++r) {
    if (memcmp(ASSUME_ALIGNED(CONST_ROW(first, r), 32), ASSUME_ALIGNED(CONST_ROW(second, r), 32),
               sizeof(word) * width) != 0) {
      return false;
    }
  }

  return true;
}

#if defined(WITH_M4RI_UTILS)
static inline mzd_local_t* mzd_convert(const mzd_t* v) {
  mzd_local_t* r = mzd_local_init(v->nrows, v->ncols);

  for (rci_t i = 0; i < v->nrows; ++i) {
    memcpy(ROW(r, i), v->rows[i], v->width * sizeof(word));
  }

  return r;
}
#endif
