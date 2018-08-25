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

#include "picnic.h"
#include "timing.h"
#include "bench_timing.h"
#include "bench_utils.h"

#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef VERBOSE
static void print_timings(timing_and_size_t* timings, unsigned int iter) {
  const unsigned int numt = sizeof(*timings) / sizeof(timings->data[0]);

  for (unsigned int i = 0; i < iter; i++) {
    for (unsigned int j = 0; j < numt; j++) {
      printf("%" PRIu64, timings[i].data[j]);
      if (j < numt - 1)
        printf(",");
    }
    printf("\n");
  }
}
#else
#ifdef WITH_DETAILED_TIMING
static void print_timings(timing_and_size_t* timings, unsigned int iter) {
  for (unsigned int i = 0; i != iter; ++i, ++timings) {
    printf("Setup:\n");
    printf("LowMC setup               %9" PRIu64 "\n", timings->gen.lowmc_init);
    printf("LowMC key generation      %9" PRIu64 "\n", timings->gen.keygen);
    printf("Public key computation    %9" PRIu64 "\n", timings->gen.pubkey);
    printf("\n");
    printf("Prove:\n");
    printf("MPC randomess generation  %9" PRIu64 "\n", timings->sign.rand);
    printf("MPC secret sharing        %9" PRIu64 "\n", timings->sign.secret_sharing);
    printf("MPC LowMC encryption      %9" PRIu64 "\n", timings->sign.lowmc_enc);
    printf("Hashing views             %9" PRIu64 "\n", timings->sign.views);
    printf("Generating challenge      %9" PRIu64 "\n", timings->sign.challenge);
    printf("Overall hash time         %6" PRIu64 "\n", timings->sign.hash);
    printf("\n");
    printf("Verify:\n");
    printf("Recomputing challenge     %9" PRIu64 "\n", timings->verify.challenge);
    printf("Verifying output shares   %9" PRIu64 "\n", timings->verify.output_shares);
    printf("Comparing output views    %9" PRIu64 "\n", timings->verify.output_views);
    printf("Verifying views           %9" PRIu64 "\n", timings->verify.verify);
    printf("Overall hash time         %9" PRIu64 "\n", timings->verify.hash);
    printf("\n");
  }
}
#else
static void print_timings(timing_and_size_t* timings, unsigned int iter) {
  for (unsigned int i = 0; i != iter; ++i, ++timings) {
    printf("Sign                      %9" PRIu64 "\n", timings->sign);
    printf("Verify                    %9" PRIu64 "\n", timings->verify);
    printf("\n");
  }
}
#endif
#endif

static void bench_sign_and_verify(const bench_options_t* options) {
  static const uint8_t m[] = {1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15, 16,
                              17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32};

  timing_and_size_t* timings_fis = calloc(options->iter, sizeof(timing_and_size_t));

  const size_t max_signature_size = picnic_signature_size(options->params);
  if (!max_signature_size) {
    printf("Failed to create Picnic instance.\n");
    return;
  }

  uint8_t sig[PICNIC_MAX_SIGNATURE_SIZE];

  timing_context_t ctx;
  if (!timing_init(&ctx)) {
    printf("Failed to initialize timing functionality.\n");
    return;
  }

  for (unsigned int i = 0; i != options->iter; ++i) {
#ifndef WITH_DETAILED_TIMING
    timing_and_size_t* timing_and_size;
    uint64_t start_time = timing_read(&ctx);
#endif
    timing_and_size           = &timings_fis[i];
    timing_and_size->max_size = max_signature_size;

    picnic_privatekey_t private_key;
    picnic_publickey_t public_key;

    if (picnic_keygen(options->params, &public_key, &private_key)) {
      printf("Failed to create key.\n");
      break;
    }

#ifndef WITH_DETAILED_TIMING
    uint64_t tmp_time       = timing_read(&ctx);
    timing_and_size->keygen = tmp_time - start_time;
    start_time              = timing_read(&ctx);
#endif
    size_t siglen = max_signature_size;
    if (!picnic_sign(&private_key, m, sizeof(m), sig, &siglen)) {
#ifndef WITH_DETAILED_TIMING
      tmp_time              = timing_read(&ctx);
      timing_and_size->sign = tmp_time - start_time;
      timing_and_size->size = siglen;
      start_time            = timing_read(&ctx);
#endif

      if (picnic_verify(&public_key, m, sizeof(m), sig, siglen)) {
        printf("picnic_verify: failed\n");
      }
#ifndef WITH_DETAILED_TIMING
      tmp_time                = timing_read(&ctx);
      timing_and_size->verify = tmp_time - start_time;
#endif
    } else {
      printf("picnic_sign: failed\n");
    }
  }

#ifdef VERBOSE
  printf("Picnic signature:\n\n");
#endif
  timing_close(&ctx);
  print_timings(timings_fis, options->iter);

  free(timings_fis);
}

int main(int argc, char** argv) {
  bench_options_t opts = { 0 };
  int ret = parse_args(&opts, argc, argv) ? 0 : -1;

  if (!ret) {
    bench_sign_and_verify(&opts);
  }

  return ret;
}
