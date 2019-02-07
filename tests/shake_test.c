/*
 *  This file is part of the optimized implementation of the Picnic signature scheme.
 *  See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */

#include <stdio.h>
#include <stdint.h>
#include "KeccakHash.h"
#include "KeccakHashtimes4.h"

//#define BENCH_SHAKE
#ifdef BENCH_SHAKE
#include "tools/bench_timing.h"
#endif

int main(void) {
  const uint8_t data1[8] = {0x01, 0x23, 0x45, 0x67, 0x89, 0xab, 0xcd, 0xef};
  const uint8_t data2[8] = {0xfe, 0xdc, 0xba, 0x98, 0x76, 0x54, 0x32, 0x10};

  const uint8_t* data_ptrs[4] = {data1, data2, data2, data1};
  uint8_t output1[32] = {0,};
  uint8_t output2[32] = {0,};
  uint8_t output3[32] = {0,};
  uint8_t output4[32] = {0,};
  uint8_t* output_ptrs[4] = {output1, output2, output3, output4};

  Keccak_HashInstancetimes4 ctx;
  if(Keccak_HashInitializetimes4_SHAKE256(&ctx))
    perror("Initialization 4x failed.");
  if(Keccak_HashUpdatetimes4(&ctx, data_ptrs, 64))
    perror("Update 4x failed");
  if(Keccak_HashFinaltimes4(&ctx, NULL))
    perror("Final 4x failed");
  if(Keccak_HashSqueezetimes4(&ctx, output_ptrs, 256))
    perror("Squeeze 4x failed");

  uint8_t out1[32] = {0,};
  uint8_t out2[32] = {0,};
  Keccak_HashInstance ctx1;
  if(Keccak_HashInitialize_SHAKE256(&ctx1))
    perror("Init 1 failed");
  if(Keccak_HashUpdate(&ctx1, data1, 64))
    perror("Update 1 failed");
  if(Keccak_HashFinal(&ctx1, NULL))
    perror("Final 1 failed");
  if(Keccak_HashSqueeze(&ctx1, out1, 256))
    perror("Squeeze 1 failed");

  Keccak_HashInstance ctx2;
  if(Keccak_HashInitialize_SHAKE256(&ctx2))
    perror("Init 2 failed");
  if(Keccak_HashUpdate(&ctx2, data2, 64))
    perror("Update 2 failed");
  if(Keccak_HashFinal(&ctx2, NULL))
    perror("Final 2 failed");
  if(Keccak_HashSqueeze(&ctx2, out2, 256))
    perror("Squeeze 2 failed");


  const uint8_t expected_output[32] = {0xc3, 0x1f, 0xcb, 0xe0, 0x76, 0x48, 0x87, 0x58,
                                       0x73, 0x0d, 0xa7, 0xe5, 0xf8, 0x96, 0x4a, 0xfe,
                                       0xfa, 0x65, 0x9a, 0x6f, 0x52, 0x6b, 0x6f, 0x9c,
                                       0xa6, 0x7e, 0x59, 0x9f, 0x31, 0x8a, 0x7e, 0xa1};

  printf("Output1: ");
  for(uint32_t i = 0; i < 32; i++){
    printf("%02X", output1[i]);
  }
  printf("\n");
  printf("Output2: ");
  for(uint32_t i = 0; i < 32; i++){
    printf("%02X", output2[i]);
  }
  printf("\n");
  printf("Output3: ");
  for(uint32_t i = 0; i < 32; i++){
    printf("%02X", output3[i]);
  }
  printf("\n");
  printf("Output4: ");
  for(uint32_t i = 0; i < 32; i++){
    printf("%02X", output4[i]);
  }
  printf("\n");


  printf("Out1:    ");
  for(uint32_t i = 0; i < 32; i++){
    printf("%02X", out1[i]);
  }
  printf("\n");
  printf("Out2:    ");
  for(uint32_t i = 0; i < 32; i++){
    printf("%02X", out2[i]);
  }
  printf("\n");

#ifdef BENCH_SHAKE
#define BENCH_ITERS (1000000)
  timing_context_t timer;
  timing_init(&timer);

  uint64_t start1 = timing_read(&timer);
  for(unsigned int counter = 0; counter < BENCH_ITERS; counter ++) {
    Keccak_HashInstancetimes4 ctx_bench_x4;
    Keccak_HashInitializetimes4_SHAKE256(&ctx_bench_x4),
    Keccak_HashUpdatetimes4(&ctx_bench_x4, data_ptrs, 64);
    Keccak_HashFinaltimes4(&ctx_bench_x4, NULL);
    Keccak_HashSqueezetimes4(&ctx_bench_x4, output_ptrs, 256);
  }
  uint64_t end1 = timing_read(&timer);

  uint64_t start2 = timing_read(&timer);
  for(unsigned int counter = 0; counter < BENCH_ITERS * 4; counter ++) {
    Keccak_HashInstance ctx_bench_x1;
    Keccak_HashInitialize_SHAKE256(&ctx_bench_x1),
    Keccak_HashUpdate(&ctx_bench_x1, data1, 64);
    Keccak_HashFinal(&ctx_bench_x1, NULL);
    Keccak_HashSqueeze(&ctx_bench_x1, out1, 256);
  }
  uint64_t end2 = timing_read(&timer);
  printf("Timings:\n 4x: %zu\n", end1-start1);
  printf(" 1x: %zu\n", end2-start2);
#endif

}
