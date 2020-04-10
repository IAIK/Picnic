/*
 *  This file is part of the optimized implementation of the Picnic signature scheme.
 *  See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */
#ifndef PICNIC_TEST_UTILS_H
#define PICNIC_TEST_UTILS_H

#include "picnic.h"
#include "mzd_additional.h"

#include <errno.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#if defined(_WIN32) || defined(_WIN64)
#define strcasecmp(s1, s2) _stricmp((s1), (s2))
#endif

static inline picnic_params_t argument_to_params(const char* arg, bool support_m1) {
  for (unsigned int param = Picnic_L1_FS; param < PARAMETER_SET_MAX_INDEX; ++param) {
    if (!strcasecmp(arg, picnic_get_param_name(param))) {
      return param;
    }
  }

  const long idx = strtol(arg, NULL, 10);
  if ((errno == ERANGE && (idx == LONG_MAX || idx == LONG_MIN)) || (errno != 0 && idx == 0) ||
      idx < 1 || (size_t)idx >= PARAMETER_SET_MAX_INDEX) {
    return PARAMETER_SET_INVALID;
  }
  if (!support_m1 && (size_t)idx > Picnic3_L5) {
    return PARAMETER_SET_INVALID;
  }

  return idx;
}

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

#endif
