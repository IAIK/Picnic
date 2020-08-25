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

#include <stdio.h>
#include <stdint.h>

#include "../kdf_shake.h"

static int test_shake128_x4(void) {
  const uint8_t data1[8] = {0x01, 0x23, 0x45, 0x67, 0x89, 0xab, 0xcd, 0xef};
  const uint8_t data2[8] = {0xfe, 0xdc, 0xba, 0x98, 0x76, 0x54, 0x32, 0x10};
  const uint8_t* data_ptrs[4] = {data1, data2, data2, data1};

  uint8_t output1[32]         = {0};
  uint8_t output2[32]         = {0};
  uint8_t output3[32]         = {0};
  uint8_t output4[32]         = {0};
  uint8_t* output_ptrs[4]     = {output1, output2, output3, output4};

  hash_context_x4 ctx;
  hash_init_x4(&ctx, 32);
  hash_update_x4(&ctx, data_ptrs, 8);
  hash_final_x4(&ctx);
  hash_squeeze_x4(&ctx, output_ptrs, 32);

  uint8_t out1[32] = {0};
  uint8_t out2[32] = {0};
  hash_context ctx1;
  hash_init(&ctx1, 32);
  hash_update(&ctx1, data1, 8);
  hash_final(&ctx1);
  hash_squeeze(&ctx1, out1, 32);

  hash_context ctx2;
  hash_init(&ctx2, 32);
  hash_update(&ctx2, data2, 8);
  hash_final(&ctx2);
  hash_squeeze(&ctx2, out2, 32);

  const uint8_t expected_output1[32] = {0xf1, 0xf7, 0x2b, 0xe3, 0x41, 0x61, 0xe0, 0x8a,
                                        0x80, 0xff, 0x73, 0x7f, 0x8c, 0xac, 0xa2, 0x30,
                                        0x25, 0x72, 0xb7, 0x9a, 0x13, 0xc5, 0x39, 0xff,
                                        0x3a, 0x28, 0xd1, 0x2f, 0x93, 0x48, 0xdc, 0x15};
  const uint8_t expected_output2[32] = {0xa,  0xc4, 0x6d, 0x96, 0x5a, 0x19, 0x2a, 0x9b,
                                        0x5b, 0x25, 0x18, 0x7d, 0x38, 0xf9, 0xa1, 0xc3,
                                        0x9,  0xe7, 0x82, 0x4c, 0x13, 0xbb, 0xf2, 0x83,
                                        0xa1, 0x17, 0xac, 0xf4, 0x55, 0x3e, 0xff, 0x14};

  if (memcmp(out1, expected_output1, 32) != 0) return -1;
  if (memcmp(out2, expected_output2, 32) != 0)
    return -1;

  if (memcmp(output1, expected_output1, 32) != 0)
    return -1;
  if (memcmp(output2, expected_output2, 32) != 0)
    return -1;
  if (memcmp(output3, expected_output2, 32) != 0)
    return -1;
  if (memcmp(output4, expected_output1, 32) != 0)
    return -1;

  return 0;
}

static int test_shake256_x4(void) {
  const uint8_t data1[8] = {0x01, 0x23, 0x45, 0x67, 0x89, 0xab, 0xcd, 0xef};
  const uint8_t data2[8] = {0xfe, 0xdc, 0xba, 0x98, 0x76, 0x54, 0x32, 0x10};
  const uint8_t* data_ptrs[4] = {data1, data2, data2, data1};

  uint8_t output1[32]         = {0};
  uint8_t output2[32]         = {0};
  uint8_t output3[32]         = {0};
  uint8_t output4[32]         = {0};
  uint8_t* output_ptrs[4]     = {output1, output2, output3, output4};

  hash_context_x4 ctx;
  hash_init_x4(&ctx, 64);
  hash_update_x4(&ctx, data_ptrs, 8);
  hash_final_x4(&ctx);
  hash_squeeze_x4(&ctx, output_ptrs, 32);

  uint8_t out1[32] = {0};
  uint8_t out2[32] = {0};
  hash_context ctx1;
  hash_init(&ctx1, 64);
  hash_update(&ctx1, data1, 8);
  hash_final(&ctx1);
  hash_squeeze(&ctx1, out1, 32);

  hash_context ctx2;
  hash_init(&ctx2, 64);
  hash_update(&ctx2, data2, 8);
  hash_final(&ctx2);
  hash_squeeze(&ctx2, out2, 32);

  const uint8_t expected_output1[32] = {0xd5, 0x9c, 0x54, 0x97, 0x73, 0xd2, 0xb3, 0x91,
                                        0x57, 0x49, 0x76, 0xb9, 0x59, 0xeb, 0xf6, 0xe0,
                                        0xcf, 0x23, 0xdd, 0xf3, 0x01, 0xb4, 0x55, 0xac,
                                        0x60, 0x0d, 0xb8, 0xe5, 0xa3, 0xd3, 0xab, 0xf1};
  const uint8_t expected_output2[32] = {0x00, 0xd2, 0xd4, 0x1b, 0xe6, 0xc7, 0xc2, 0xf4,
                                        0x7f, 0xf9, 0x4c, 0x7d, 0x1c, 0x50, 0x70, 0x8c,
                                        0xe3, 0x8d, 0xc7, 0x36, 0x4c, 0xbd, 0xa3, 0xb7,
                                        0xfb, 0x66, 0xfc, 0x05, 0x25, 0xfe, 0xef, 0xc8};

  if (memcmp(out1, expected_output1, 32) != 0)
    return -1;
  if (memcmp(out2, expected_output2, 32) != 0)
    return -1;

  if (memcmp(output1, expected_output1, 32) != 0)
    return -1;
  if (memcmp(output2, expected_output2, 32) != 0)
    return -1;
  if (memcmp(output3, expected_output2, 32) != 0)
    return -1;
  if (memcmp(output4, expected_output1, 32) != 0)
    return -1;

  return 0;
}

int main(void) {
  int ret = 0;

  printf("testing SHAKE128x4 ...");
  int r = test_shake128_x4();
  if (r == -2) {
    printf("SKIPPED\n");
  } else if (r) {
    printf("FAILED\n");
    ret = -1;
  } else {
    printf("OK\n");
  }

  printf("testing SHAKE256x4 ...");
  r = test_shake256_x4();
  if (r == -2) {
    printf("SKIPPED\n");
  } else if (r) {
    printf("FAILED\n");
    ret = -1;
  } else {
    printf("OK\n");
  }

  return ret;
}
