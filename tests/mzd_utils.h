/*
 *  This file is part of the optimized implementation of the Picnic signature scheme.
 *  See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */
#ifndef PICNIC_TEST_MZD_UTILS_H
#define PICNIC_TEST_MZD_UTILS_H

#include "mzd_additional.h"

#include <string.h>
#include <stdbool.h>

static inline bool mzd_local_equal(mzd_local_t const* first, mzd_local_t const* second,
                                   unsigned int rows, unsigned cols) {
  const unsigned int num_uints  = (cols + 63) / 64;
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

#endif
