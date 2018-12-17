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

#include "../lowmc.h"
#include "../picnic_impl.h"

#include <m4ri/m4ri.h>
#include <stdint.h>

#include "utils.h"

static int lowmc_enc_str(const picnic_params_t param, const char* key, const char* plaintext,
                         const char* expected) {
  const picnic_instance_t* pp = picnic_instance_get(param);
  if (!pp) {
    return -1;
  }
  const lowmc_t* lowmc = pp->lowmc;

  mzd_t* sk = mzd_from_str(1, lowmc->k, key);
  mzd_t* pt = mzd_from_str(1, lowmc->n, plaintext);
  mzd_t* ct = mzd_from_str(1, lowmc->n, expected);

  mzd_local_t* skl = mzd_convert(sk);
  mzd_local_t* ptl = mzd_convert(pt);
  mzd_local_t* ctl = mzd_convert(ct);

  int ret          = 0;
  mzd_local_t* ctr = pp->impls.lowmc(skl, ptl);
  if (!ctr) {
    ret = 1;
    goto end;
  }

  if (!mzd_local_equal(ctr, ctl)) {
    ret = 2;
  }

end:
  mzd_local_free(ctr);
  mzd_local_free(ctl);
  mzd_local_free(ptl);
  mzd_local_free(skl);
  mzd_free(ct);
  mzd_free(pt);
  mzd_free(sk);

  return ret;
}

static const char key_L1_1[] = "0000000000000000000000000000000000000000000000000000000000000000000"
                               "0000100100011010001010110011110001001101010111100110111101111";
static const char plaintext_L1_1[] = "0000000000000000000000000000000000000000000000000000000000000"
                                     "0001111111011011100101110101001100001110110010101000011001000"
                                     "010000";
static const char expected_L1_1[] = "01111000101001001001011000101001011010110110010000111100000011"
                                    "10101010100100000011100110110000111000011110110111011000001100"
                                    "0000";

static const char key_L1_2[] = "0000000000000000000000000000000000000000000000000000000000000000111"
                               "1111011011100101110101001100001110110010101000011001000010000";
static const char plaintext_L1_2[] = "0000000000000000000000000000000000000000000000000000000000000"
                                     "0000000000100100011010001010110011110001001101010111100110111"
                                     "101111";
static const char expected_L1_2[] = "10110010110111000100101000011101111111111100001010100110111000"
                                    "11011011010110001001100001111010111000010011000110001110101000"
                                    "1001";

static const char key_256_256_10_38[] =
    "0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"
    "0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"
    "00000000000000000000000000000000000000000000000000000000000000000001";
static const char expected_256_256_10_38[] =
    "1010001100111111101000101100011100101000011100101100111101110110111001111011101010001001100010"
    "1101110110101000100101101001100010010010000110101001110000000100010000011111000001111110010001"
    "01100100011111110011111011000110100001110000001111001010110111001000";

static const char key_192_192_10_30[] = "0000000000000000000000000000000000000000000000000000000000"
                                        "0000000000000000000000000000000000000000000000000000000000"
                                        "0000000000000000000000000000000000000000000000000000000000"
                                        "000000000000000001";
static const char expected_192_192_10_30[] = "10111000100100100001101010010110000100011010000111101"
                                             "00001011011100111100010101000010100001010001110000100"
                                             "00011011001000000010110000010010010110011001110110111"
                                             "011101100111000011111001000111110";

static const struct {
  picnic_params_t param;
  const char* key;
  const char* plaintext;
  const char* expected;
} str_tests[] = {{Picnic_L1_FS, key_L1_1, plaintext_L1_1, expected_L1_1},
                 {Picnic_L1_FS, key_L1_2, plaintext_L1_2, expected_L1_2},
                 {Picnic_L3_FS, key_192_192_10_30, key_192_192_10_30, expected_192_192_10_30},
                 {Picnic_L5_FS, key_256_256_10_38, key_256_256_10_38, expected_256_256_10_38}};

static const size_t num_str_tests = sizeof(str_tests) / sizeof(str_tests[0]);

int main(void) {
  int ret = 0;
  for (size_t s = 0; s < num_str_tests; ++s) {
    const int t = lowmc_enc_str(str_tests[s].param, str_tests[s].key, str_tests[s].plaintext,
                                str_tests[s].expected);
    if (t) {
      printf("ERR: lowmc_enc_str %zu FAILED (%d)\n", s, t);
      ret = -1;
    }
  }

  return ret;
}
