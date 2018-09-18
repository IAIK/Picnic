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

#include "bench_timing.h"
#include "bench_utils.h"
#include "io.h"
#include "lowmc.h"
#include "picnic_impl.h"
#include "randomness.h"

#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

static void print_timings(const uint64_t* timing, unsigned int iter) {
  for (unsigned int i = 0; i < iter; i++) {
    printf("%" PRIu64 "\n", timing[i]);
  }
}

static void bench_lowmc(const bench_options_t* options) {
  uint64_t* timings = calloc(options->iter, sizeof(uint64_t));

  timing_context_t ctx;
  if (!timing_init(&ctx)) {
    printf("Failed to initialize timing functionality.\n");
    return;
  }

  lowmc_t custom_lowmc = { 0 };
  lowmc_implementation_f custom_impl = NULL;
#if defined(WITH_CUSTOM_INSTANCES)
  if (options->lowmc_file) {
    if (!lowmc_read_file(&custom_lowmc, options->lowmc_file)) {
      printf("Failed to read LowMC instance.\n");
      return;
    }

    custom_impl = lowmc_get_implementation(&custom_lowmc);
  }
#endif

  const picnic_instance_t* pp = !custom_impl ? picnic_instance_get(options->params) : NULL;
  if (!custom_impl && !pp) {
    printf("Failed to initialize LowMC instance.\n");
    return;
  }

  const lowmc_t* lowmc                    = pp ? pp->lowmc : &custom_lowmc;
  const lowmc_implementation_f lowmc_impl = pp ? pp->lowmc_impl : custom_impl;

  mzd_local_t* sk = mzd_local_init(1, lowmc->k);
  mzd_local_t* pt = mzd_local_init(1, lowmc->n);

  const size_t input_size = (lowmc->k + 7) >> 3;
  const size_t output_size = (lowmc->n + 7) >> 3;

  uint8_t* rand = malloc(input_size + output_size);
  rand_bytes(rand, input_size + output_size);
  mzd_from_char_array(sk, rand, input_size);
  mzd_from_char_array(pt, rand + input_size, output_size);
  free(rand);

  for (unsigned int i = 0; i != options->iter; ++i) {
    const uint64_t start_time = timing_read(&ctx);

    mzd_local_t* ct = lowmc_impl(lowmc, sk, pt);
    mzd_local_free(pt);
    pt = ct;

    const uint64_t end_time = timing_read(&ctx);
    timings[i]              = end_time - start_time;
  }

  mzd_local_free(pt);
  mzd_local_free(sk);

  timing_close(&ctx);
  print_timings(timings, options->iter);

#if defined(WITH_CUSTOM_INSTANCES)
  if (custom_impl) {
    lowmc_clear(&custom_lowmc);
  }
#endif

  free(timings);
}

int main(int argc, char** argv) {
  bench_options_t opts = {0};
  int ret              = parse_args(&opts, argc, argv) ? 0 : -1;

  if (!ret) {
    bench_lowmc(&opts);
  }

  return ret;
}
