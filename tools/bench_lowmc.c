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
#include "../io.h"
#include "../lowmc.h"
#include "../picnic_instances.h"
#include "../randomness.h"

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
  const picnic_instance_t* pp = picnic_instance_get(options->params);
  if (!pp) {
    printf("Failed to initialize LowMC instance.\n");
    return;
  }

  timing_context_t ctx;
  if (!timing_init(&ctx)) {
    printf("Failed to initialize timing functionality.\n");
    return;
  }

  uint64_t* timings = calloc(options->iter, sizeof(uint64_t));

  const lowmc_parameters_t* lowmc         = &pp->lowmc;
  const lowmc_implementation_f lowmc_impl = pp->impl_lowmc;

  mzd_local_t* sk = mzd_local_init(1, lowmc->k);
  mzd_local_t* pt = mzd_local_init(1, lowmc->n);
  mzd_local_t* ct = mzd_local_init(1, lowmc->n);

  const size_t input_size = (lowmc->k + 7) >> 3;
  const size_t output_size = (lowmc->n + 7) >> 3;

  uint8_t* rand = malloc(input_size + output_size);
  rand_bits(rand, (input_size + output_size) * 8);
  mzd_from_char_array(sk, rand, input_size);
  mzd_from_char_array(pt, rand + input_size, output_size);
  free(rand);

  for (unsigned int i = 0; i != options->iter; ++i) {
    const uint64_t start_time = timing_read(&ctx);
    lowmc_impl(sk, pt, ct);
    timings[i] = timing_read(&ctx) - start_time;

    mzd_local_t* tmp = pt;
    pt = ct;
    ct = tmp;
  }

  mzd_local_free(ct);
  mzd_local_free(pt);
  mzd_local_free(sk);

  timing_close(&ctx);
  print_timings(timings, options->iter);

  free(timings);
}

int main(int argc, char** argv) {
  bench_options_t opts = {PARAMETER_SET_INVALID, 0};
  int ret              = parse_args(&opts, argc, argv) ? 0 : -1;

  if (!ret) {
    bench_lowmc(&opts);
  }

  return ret;
}
