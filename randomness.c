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

#include "randomness.h"

#if defined(HAVE_RANDOMBYTES) || defined(SUPERCOP)
extern void randombytes(unsigned char* x, unsigned long long xlen);

int rand_bytes(uint8_t* dst, size_t len) {
  randombytes(dst, len);
  return 1;
}
#else

#if defined(__linux__)
#if defined(HAVE_LINUX_RANDOM_H) && defined(HAVE_GETRANDOM)
#include <linux/random.h>

int rand_bytes(uint8_t* dst, size_t len) {
  if (getrandom(dst, len, GRND_NONBLOCK) != len) {
    return 0;
  }
  return 1;
}
#else
#include <stdio.h>

int rand_bytes(uint8_t* dst, size_t len) {
  FILE* urandom = fopen("/dev/urandom", "r");
  int ret       = 1;
  if (urandom) {
    ret = fread(dst, 1, len, urandom) == len ? 1 : 0;
    fclose(urandom);
  }
  return ret;
}
#endif
#elif defined(_WIN32)
#include <windows.h>

int rand_bytes(uint8_t* dst, size_t len) {
  if (!BCRYPT_SUCCESS(BCryptGenRandom(NULL, dst, len, BCRYPT_USE_SYSTEM_PREFERRED_RNG))) {
    return 0;
  }
  return 1;
}
#else
#error "Unsupported OS! Please implement rand_bytes."
#endif
#endif
