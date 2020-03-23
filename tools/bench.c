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
#include "picnic.h"

#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct {
  uint64_t keygen, sign, verify, size, max_size;
} timing_and_size_t;

static void print_timings(timing_and_size_t* timings, unsigned int iter) {
  for (unsigned int i = 0; i < iter; i++) {
    printf("%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64 "\n", timings[i].keygen,
           timings[i].sign, timings[i].verify, timings[i].size, timings[i].max_size);
  }
}

static void bench_sign_and_verify(const bench_options_t* options) {
  static const uint8_t m[] = {1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15, 16,
                              17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32};

  const size_t max_signature_size = picnic_signature_size(options->params);
  if (!max_signature_size) {
    printf("Failed to create Picnic instance.\n");
    return;
  }

  timing_context_t ctx;
  if (!timing_init(&ctx)) {
    printf("Failed to initialize timing functionality.\n");
    return;
  }

  timing_and_size_t* timings = calloc(options->iter, sizeof(timing_and_size_t));
  uint8_t sig[PICNIC_MAX_SIGNATURE_SIZE];

  for (unsigned int i = 0; i != options->iter; ++i) {
    timing_and_size_t* timing = &timings[i];
    timing->max_size          = max_signature_size;

    uint64_t start_time = timing_read(&ctx);
    picnic_privatekey_t private_key;
    picnic_publickey_t public_key;

    if (picnic_keygen(options->params, &public_key, &private_key)) {
      printf("picnic_keygen: failed.\n");
      break;
    }

    uint64_t tmp_time = timing_read(&ctx);
    timing->keygen    = tmp_time - start_time;
    start_time        = timing_read(&ctx);

    size_t siglen = max_signature_size;
    if (!picnic_sign(&private_key, m, sizeof(m), sig, &siglen)) {
      tmp_time     = timing_read(&ctx);
      timing->sign = tmp_time - start_time;
      timing->size = siglen;
      start_time   = timing_read(&ctx);

      if (picnic_verify(&public_key, m, sizeof(m), sig, siglen)) {
        printf("picnic_verify: failed\n");
      }
      tmp_time       = timing_read(&ctx);
      timing->verify = tmp_time - start_time;
    } else {
      printf("picnic_sign: failed\n");
    }
  }

  timing_close(&ctx);
  print_timings(timings, options->iter);

  free(timings);
}

int main(int argc, char** argv) {
  bench_options_t opts = {PARAMETER_SET_INVALID, 0};
  int ret              = parse_args(&opts, argc, argv) ? 0 : -1;

  if (!ret) {
    bench_sign_and_verify(&opts);
  }

  return ret;
}
