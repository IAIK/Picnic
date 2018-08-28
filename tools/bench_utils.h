/*
 *  This file is part of the optimized implementation of the Picnic signature scheme.
 *  See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */

#ifndef BENCH_UTILS_H
#define BENCH_UTILS_H

#include <stdbool.h>
#include <stdint.h>

#include "picnic.h"

typedef struct {
  picnic_params_t params;
  const char* lowmc_file;
  uint32_t iter;
  bool lowmc;
} bench_options_t;

bool parse_args(bench_options_t* options, int argc, char** argv);

#endif
