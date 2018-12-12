/*
 *  This file is part of the optimized implementation of the Picnic signature scheme.
 *  See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */

#include "api.h"

#include <stdio.h>
#include <string.h>

#if defined(__WIN32__)
#define SIZET_FMT "%Iu"
#define LL_FMT "%I64u"
#else
#define SIZET_FMT "%zu"
#define LL_FMT "%llu"
#endif

int main(void) {
  unsigned char pk[CRYPTO_PUBLICKEYBYTES]          = {0};
  unsigned char sk[CRYPTO_SECRETKEYBYTES]          = {0};
  const unsigned char message[50]                  = {0};
  unsigned char omessage[sizeof(message)]          = {0};
  unsigned char sm[sizeof(message) + CRYPTO_BYTES] = {0};

  int ret = crypto_sign_keypair(pk, sk);
  if (ret != 0) {
    printf("Failed to generate key pair\n");
    return -1;
  }

  unsigned long long smlen = sizeof(sm);
  ret                      = crypto_sign(sm, &smlen, message, sizeof(message), sk);
  if (ret != 0) {
    printf("Failed to sign\n");
    return -1;
  }

  unsigned long long mlen = sizeof(omessage);
  ret                     = crypto_sign_open(omessage, &mlen, sm, smlen, pk);
  if (ret != 0) {
    printf("Failed to verify (ret = %d)\n", ret);
    return -1;
  }

  if (mlen != sizeof(message)) {
    printf("length of message after verify incorrect, got " LL_FMT ", expected " SIZET_FMT "\n",
           mlen, sizeof(message));
    return -1;
  }
  if (memcmp(message, omessage, sizeof(message)) != 0) {
    printf("message mismatch after verification\n");
    return -1;
  }

  // test special case where message and signature overlap
  memcpy(sm, message, sizeof(message));

  smlen = sizeof(sm);
  ret   = crypto_sign(sm, &smlen, sm, sizeof(message), sk);
  if (ret != 0) {
    printf("Failed to sign\n");
    return -1;
  }

  mlen = smlen;
  ret  = crypto_sign_open(sm, &mlen, sm, smlen, pk);
  if (ret != 0) {
    printf("Failed to verify (ret = %d)\n", ret);
    return -1;
  }

  if (mlen != sizeof(message)) {
    printf("length of message after verify incorrect, got " LL_FMT ", expected " SIZET_FMT "\n",
           mlen, sizeof(message));
    return -1;
  }
  if (memcmp(message, sm, sizeof(message)) != 0) {
    printf("message mismatch after verification\n");
    return -1;
  }

  printf("Sign/Verify test passed\n");

  return 0;
}
