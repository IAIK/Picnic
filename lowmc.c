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

#include "io.h"
#include "lowmc.h"
#include "mzd_additional.h"
#if defined(WITH_KKW)
#include "picnic2_impl.h"
#include "picnic2_simulate_mul.h"
#endif
#if defined(WITH_OPT)
#include "simd.h"
#endif

#if !defined(_MSC_VER)
#include <stdalign.h>
#endif
#include <string.h>

// TODO: these sboxes are pretty inefficient, look for better implementation maybe
/**
 * S-box for m = 42
 */
static void sbox_layer_42(mzd_local_t* in) {
  uint8_t tmp[16];
  mzd_to_char_array(tmp, in, 16);
  for (uint32_t i = 0; i < 42; i++) {
    uint8_t a = getBit(tmp, 3 * i + 2);
    uint8_t b = getBit(tmp, 3 * i + 1);
    uint8_t c = getBit(tmp, 3 * i + 0);

    uint8_t d = a ^ (b & c);
    uint8_t e = a ^ b ^ (a & c);
    uint8_t f = a ^ b ^ c ^ (a & b);

    setBit(tmp, 3 * i + 2, d);
    setBit(tmp, 3 * i + 1, e);
    setBit(tmp, 3 * i + 0, f);
  }
  mzd_from_char_array(in, tmp, 16);
}
/**
 * S-box for m = 64
 */
static void sbox_layer_64(mzd_local_t* in) {
  uint8_t tmp[24];
  mzd_to_char_array(tmp, in, 24);
  for (uint32_t i = 0; i < 64; i++) {
    uint8_t a = getBit(tmp, 3 * i + 2);
    uint8_t b = getBit(tmp, 3 * i + 1);
    uint8_t c = getBit(tmp, 3 * i + 0);

    uint8_t d = a ^ (b & c);
    uint8_t e = a ^ b ^ (a & c);
    uint8_t f = a ^ b ^ c ^ (a & b);

    setBit(tmp, 3 * i + 2, d);
    setBit(tmp, 3 * i + 1, e);
    setBit(tmp, 3 * i + 0, f);
  }
  mzd_from_char_array(in, tmp, 24);
}
/**
 * S-box for m = 85
 */
static void sbox_layer_85(mzd_local_t* in) {
  uint8_t tmp[32];
  mzd_to_char_array(tmp, in, 32);
  for (uint32_t i = 0; i < 85; i++) {
    uint8_t a = getBit(tmp, 3 * i + 2);
    uint8_t b = getBit(tmp, 3 * i + 1);
    uint8_t c = getBit(tmp, 3 * i + 0);

    uint8_t d = a ^ (b & c);
    uint8_t e = a ^ b ^ (a & c);
    uint8_t f = a ^ b ^ c ^ (a & b);

    setBit(tmp, 3 * i + 2, d);
    setBit(tmp, 3 * i + 1, e);
    setBit(tmp, 3 * i + 0, f);
  }
  mzd_from_char_array(in, tmp, 32);
}


#if defined(WITH_LOWMC_126_126_4)
#include "lowmc_126_126_4.h"
#endif
#if defined(WITH_LOWMC_192_192_4)
#include "lowmc_192_192_4.h"
#endif
#if defined(WITH_LOWMC_255_255_4)
#include "lowmc_255_255_4.h"
#endif

#if !defined(NO_UINT64_FALLBACK)
// uint64 based implementation
#include "lowmc_fns_uint64_L1.h"
#define LOWMC lowmc_uint64_128
#include "lowmc.c.i"

#include "lowmc_fns_uint64_L3.h"
#undef LOWMC
#define LOWMC lowmc_uint64_192
#include "lowmc.c.i"

#include "lowmc_fns_uint64_L5.h"
#undef LOWMC
#define LOWMC lowmc_uint64_256
#include "lowmc.c.i"
#endif

#if defined(WITH_OPT)
#if defined(WITH_SSE2) || defined(WITH_NEON)
#if defined(WITH_SSE2)
#define FN_ATTR ATTR_TARGET_SSE2
#endif

// L1 using SSE2/NEON
#include "lowmc_fns_s128_L1.h"
#undef LOWMC
#define LOWMC lowmc_s128_128
#include "lowmc.c.i"

// L3 using SSE2/NEON
#include "lowmc_fns_s128_L3.h"
#undef LOWMC
#define LOWMC lowmc_s128_192
#include "lowmc.c.i"

// L5 using SSE2/NEON
#include "lowmc_fns_s128_L5.h"
#undef LOWMC
#define LOWMC lowmc_s128_256
#include "lowmc.c.i"

#undef FN_ATTR
#endif

#if defined(WITH_AVX2)
#define FN_ATTR ATTR_TARGET_AVX2

// L1 using AVX2
#include "lowmc_fns_s256_L1.h"
#undef LOWMC
#define LOWMC lowmc_s256_128
#include "lowmc.c.i"

// L3 using AVX2
#include "lowmc_fns_s256_L3.h"
#undef LOWMC
#define LOWMC lowmc_s256_192
#include "lowmc.c.i"

// L5 using AVX2
#include "lowmc_fns_s256_L5.h"
#undef LOWMC
#define LOWMC lowmc_s256_256
#include "lowmc.c.i"

#undef FN_ATTR

#endif
#endif

lowmc_implementation_f lowmc_get_implementation(const lowmc_t* lowmc) {
  ASSUME(lowmc->m == 10);
  ASSUME(lowmc->n == 128 || lowmc->n == 192 || lowmc->n == 256);

#if defined(WITH_OPT)
#if defined(WITH_AVX2)
  if (CPU_SUPPORTS_AVX2) {
    if (lowmc->m == 10) {
      switch (lowmc->n) {
#if defined(WITH_LOWMC_126_126_4)
      case 128:
        return lowmc_s256_128_10;
#endif
#if defined(WITH_LOWMC_192_192_4)
      case 192:
        return lowmc_s256_192_10;
#endif
#if defined(WITH_LOWMC_255_255_4)
      case 256:
        return lowmc_s256_256_10;
#endif
      }
    }
  }
#endif
#if defined(WITH_SSE2) || defined(WITH_NEON)
  if (CPU_SUPPORTS_SSE2 || CPU_SUPPORTS_NEON) {
    if (lowmc->m == 10) {
      switch (lowmc->n) {
#if defined(WITH_LOWMC_126_126_4)
      case 128:
        return lowmc_s128_128_10;
#endif
#if defined(WITH_LOWMC_192_192_4)
      case 192:
        return lowmc_s128_192_10;
#endif
#if defined(WITH_LOWMC_255_255_4)
      case 256:
        return lowmc_s128_256_10;
#endif
      }
    }
  }
#endif
#endif

#if !defined(NO_UINT64_FALLBACK)
  if (lowmc->m == 10) {
    switch (lowmc->n) {
#if defined(WITH_LOWMC_126_126_4)
    case 128:
      return lowmc_uint64_128_10;
#endif
#if defined(WITH_LOWMC_192_192_4)
    case 192:
      return lowmc_uint64_192_10;
#endif
#if defined(WITH_LOWMC_255_255_4)
    case 256:
      return lowmc_uint64_256_10;
#endif
    }
  }

#endif

  return NULL;
}

#if defined(WITH_ZKBPP)
lowmc_store_implementation_f lowmc_store_get_implementation(const lowmc_t* lowmc) {
  ASSUME(lowmc->m == 10);
  ASSUME(lowmc->n == 128 || lowmc->n == 192 || lowmc->n == 256);

#if defined(WITH_OPT)
#if defined(WITH_AVX2)
  if (CPU_SUPPORTS_AVX2) {
    if (lowmc->m == 10) {
      switch (lowmc->n) {
#if defined(WITH_LOWMC_126_126_4)
      case 128:
        return lowmc_s256_128_store_10;
#endif
#if defined(WITH_LOWMC_192_192_4)
      case 192:
        return lowmc_s256_192_store_10;
#endif
#if defined(WITH_LOWMC_255_255_4)
      case 256:
        return lowmc_s256_256_store_10;
#endif
      }
    }
  }
#endif
#if defined(WITH_SSE2) || defined(WITH_NEON)
  if (CPU_SUPPORTS_SSE2 || CPU_SUPPORTS_NEON) {
    if (lowmc->m == 10) {
      switch (lowmc->n) {
#if defined(WITH_LOWMC_126_126_4)
      case 128:
        return lowmc_s128_128_store_10;
#endif
#if defined(WITH_LOWMC_192_192_4)
      case 192:
        return lowmc_s128_192_store_10;
#endif
#if defined(WITH_LOWMC_255_255_4)
      case 256:
        return lowmc_s128_256_store_10;
#endif
      }
    }
  }
#endif
#endif

#if !defined(NO_UINT64_FALLBACK)
  if (lowmc->m == 10) {
    switch (lowmc->n) {
#if defined(WITH_LOWMC_126_126_4)
    case 128:
      return lowmc_uint64_128_store_10;
#endif
#if defined(WITH_LOWMC_192_192_4)
    case 192:
      return lowmc_uint64_192_store_10;
#endif
#if defined(WITH_LOWMC_255_255_4)
    case 256:
      return lowmc_uint64_256_store_10;
#endif
    }
  }

#endif

  return NULL;
}
#endif

#if defined(WITH_KKW)
lowmc_compute_aux_implementation_f lowmc_compute_aux_get_implementation(const lowmc_t* lowmc) {
  ASSUME(lowmc->m == 10);
  ASSUME(lowmc->n == 128 || lowmc->n == 192 || lowmc->n == 256);

#if defined(WITH_OPT)
#if defined(WITH_AVX2)
  if (CPU_SUPPORTS_AVX2) {
    if (lowmc->m == 10) {
      switch (lowmc->n) {
#if defined(WITH_LOWMC_126_126_4)
      case 128:
        return lowmc_s256_128_compute_aux_10;
#endif
#if defined(WITH_LOWMC_192_192_4)
      case 192:
        return lowmc_s256_192_compute_aux_10;
#endif
#if defined(WITH_LOWMC_255_255_4)
      case 256:
        return lowmc_s256_256_compute_aux_10;
#endif
      }
    }
  }
#endif
#if defined(WITH_SSE2) || defined(WITH_NEON)
  if (CPU_SUPPORTS_SSE2 || CPU_SUPPORTS_NEON) {
    if (lowmc->m == 10) {
      switch (lowmc->n) {
#if defined(WITH_LOWMC_126_126_4)
      case 128:
        return lowmc_s128_128_compute_aux_10;
#endif
#if defined(WITH_LOWMC_192_192_4)
      case 192:
        return lowmc_s128_192_compute_aux_10;
#endif
#if defined(WITH_LOWMC_255_255_4)
      case 256:
        return lowmc_s128_256_compute_aux_10;
#endif
      }
    }
  }
#endif
#endif

#if !defined(NO_UINT64_FALLBACK)
  if (lowmc->m == 10) {
    switch (lowmc->n) {
#if defined(WITH_LOWMC_126_126_4)
    case 128:
      return lowmc_uint64_128_compute_aux_10;
#endif
#if defined(WITH_LOWMC_192_192_4)
    case 192:
      return lowmc_uint64_192_compute_aux_10;
#endif
#if defined(WITH_LOWMC_255_255_4)
    case 256:
      return lowmc_uint64_256_compute_aux_10;
#endif
    }
  }
#endif

  return NULL;
}
#endif
