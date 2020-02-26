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

#include "mpc_lowmc.h"
#include "mzd_additional.h"

#if !defined(_MSC_VER)
#include <stdalign.h>
#endif
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#if defined(WITH_OPT)
#include "simd.h"
#endif

#define MPC_LOOP_CONST(function, result, first, second, sc)                                        \
  do {                                                                                             \
    for (unsigned int e = 0; e < (sc); ++e) {                                                      \
      function((result)[e], (first)[e], (second));                                                 \
    }                                                                                              \
  } while (0)

#define MPC_LOOP_SHARED(function, result, first, second, sc)                                       \
  do {                                                                                             \
    for (unsigned int o = 0; o < (sc); ++o) {                                                      \
      function((result)[o], (first)[o], (second)[o]);                                              \
    }                                                                                              \
  } while (0)

#define MPC_LOOP_SHARED_1(function, result, first, sc)                                             \
  do {                                                                                             \
    for (unsigned int o = 0; o < (sc); ++o) {                                                      \
      function((result)[o], (first)[o]);                                                           \
    }                                                                                              \
  } while (0)

#define MPC_LOOP_CONST_C_0(function, result, first, second, sc)                                    \
  function((result)[0], (first)[0], (second))

#define MPC_LOOP_CONST_C_ch(function, result, first, second, sc, c)                                \
  do {                                                                                             \
    if (!(c)) {                                                                                    \
      MPC_LOOP_CONST_C_0(function, result, first, second, sc);                                     \
    } else if ((c) == (sc)) {                                                                      \
      function((result)[(sc)-1], first[(sc)-1], (second));                                         \
    }                                                                                              \
  } while (0)

static void mpc_and_uint64_128(mzd_local_t* res, const mzd_local_t* first,
                               const mzd_local_t* second, const mzd_local_t* r, view_t* view,
                               unsigned viewshift) {
  mzd_local_t tmp;

  for (unsigned m = 0; m < SC_PROOF; ++m) {
    const unsigned j = (m + 1) % SC_PROOF;

    // f[m] & s[m]
    mzd_and_uint64_128(&res[m], &first[m], &second[m]);

    // f[m + 1] & s[m]
    mzd_and_uint64_128(&tmp, &first[j], &second[m]);
    mzd_xor_uint64_128(&res[m], &res[m], &tmp);

    // f[m] & s[m + 1]
    mzd_and_uint64_128(&tmp, &first[m], &second[j]);
    mzd_xor_uint64_128(&res[m], &res[m], &tmp);

    // ... ^ r[m] ^ r[m + 1]
    mzd_xor_uint64_128(&res[m], &res[m], &r[m]);
    mzd_xor_uint64_128(&res[m], &res[m], &r[j]);

    if (viewshift) {
      mzd_shift_right_uint64_128(&tmp, &res[m], viewshift);
      mzd_xor_uint64_128(&view->s[m], &view->s[m], &tmp);
    } else {
      // on first call (viewshift == 0), view->t[0..2] == 0
      memcpy(&view->s[m], &res[m], sizeof(mzd_local_t));
    }
  }

#if 0
  mpc_shift_right(buffer, res, viewshift, SC_PROOF);
  mpc_xor(view->s, view->s, buffer, SC_PROOF);
#endif
}

static void mpc_and_verify_uint64_128(mzd_local_t* res, const mzd_local_t* first,
                                      const mzd_local_t* second, const mzd_local_t* r,
                                      view_t* view, const mzd_local_t* mask, unsigned viewshift) {
  mzd_local_t tmp;

  for (unsigned m = 0; m < (SC_VERIFY - 1); ++m) {
    const unsigned j = m + 1;

    mzd_and_uint64_128(&res[m], &first[m], &second[m]);

    mzd_and_uint64_128(&tmp, &first[j], &second[m]);
    mzd_xor_uint64_128(&res[m], &res[m], &tmp);

    mzd_and_uint64_128(&tmp, &first[m], &second[j]);
    mzd_xor_uint64_128(&res[m], &res[m], &tmp);

    mzd_xor_uint64_128(&res[m], &res[m], &r[m]);
    mzd_xor_uint64_128(&res[m], &res[m], &r[j]);

    if (viewshift || m) {
      mzd_shift_right_uint64_128(&tmp, &res[m], viewshift);
      mzd_xor_uint64_128(&view->s[m], &view->s[m], &tmp);
    } else {
      // on first call (viewshift == 0), view->s[0] == 0
      memcpy(&view->s[m], &res[m], sizeof(mzd_local_t));
    }
  }

#if 0
  for (unsigned m = 0; m < (SC_VERIFY - 1); ++m) {
    mzd_shift_right(b, res[m], viewshift);
    mzd_xor(view->s[m], view->s[m], b);
  }
#endif

  mzd_shift_left_uint64_128(&res[SC_VERIFY - 1], &view->s[SC_VERIFY - 1], viewshift);
  mzd_and_uint64_128(&res[SC_VERIFY - 1], &res[SC_VERIFY - 1], mask);
}

#define bitsliced_step_1(sc, AND, ROL, MASK_A, MASK_B, MASK_C)                                     \
  mzd_local_t x2m[sc];                                                                             \
  mzd_local_t r0m[sc], r1m[sc], r2m[sc];                                                           \
  mzd_local_t x0s[sc], x1s[sc], r0s[sc], r1s[sc];                                                  \
                                                                                                   \
  for (unsigned int m = 0; m < (sc); ++m) {                                                        \
    (AND)(&x0s[m], &in[m], MASK_A);                                                                \
    (AND)(&x1s[m], &in[m], MASK_B);                                                                \
    (AND)(&x2m[m], &in[m], MASK_C);                                                                \
                                                                                                   \
    (ROL)(&x0s[m], &x0s[m], 2);                                                                    \
    (ROL)(&x1s[m], &x1s[m], 1);                                                                    \
                                                                                                   \
    (AND)(&r0m[m], &rvec->s[m], MASK_A);                                                           \
    (AND)(&r1m[m], &rvec->s[m], MASK_B);                                                           \
    (AND)(&r2m[m], &rvec->s[m], MASK_C);                                                           \
                                                                                                   \
    (ROL)(&r0s[m], &r0m[m], 2);                                                                    \
    (ROL)(&r1s[m], &r1m[m], 1);                                                                    \
  }

#define bitsliced_step_2(sc, XOR, ROR)                                                             \
  for (unsigned int m = 0; m < sc; ++m) {                                                          \
    (XOR)(&r2m[m], &r2m[m], &x0s[m]);                                                              \
    (XOR)(&x0s[m], &x0s[m], &x1s[m]);                                                              \
    (XOR)(&r1m[m], &r1m[m], &x0s[m]);                                                              \
    (XOR)(&r0m[m], &r0m[m], &x0s[m]);                                                              \
    (XOR)(&r0m[m], &r0m[m], &x2m[m]);                                                              \
                                                                                                   \
    (ROR)(&x0s[m], &r2m[m], 2);                                                                    \
    (ROR)(&x1s[m], &r1m[m], 1);                                                                    \
                                                                                                   \
    (XOR)(&x0s[m], &x0s[m], &r0m[m]);                                                              \
    (XOR)(&out[m], &x0s[m], &x1s[m]);                                                              \
  }

static const mzd_local_t mask_126_126_42_a[1] = {
    {{UINT64_C(0x4924924924924924), UINT64_C(0x2492492492492492), UINT64_C(0x0), UINT64_C(0x0)}}};
static const mzd_local_t mask_126_126_42_b[1] = {
    {{UINT64_C(0x9249249249249248), UINT64_C(0x4924924924924924), UINT64_C(0x0), UINT64_C(0x0)}}};
static const mzd_local_t mask_126_126_42_c[1] = {
    {{UINT64_C(0x2492492492492490), UINT64_C(0x9249249249249249), UINT64_C(0x0), UINT64_C(0x0)}}};

static const mzd_local_t mask_192_192_64_a[1] = {
    {{UINT64_C(0x9249249249249249), UINT64_C(0x4924924924924924), UINT64_C(0x2492492492492492),
      UINT64_C(0x0)}}};
static const mzd_local_t mask_192_192_64_b[1] = {
    {{UINT64_C(0x2492492492492492), UINT64_C(0x9249249249249249), UINT64_C(0x4924924924924924),
      UINT64_C(0x0)}}};
static const mzd_local_t mask_192_192_64_c[1] = {
    {{UINT64_C(0x4924924924924924), UINT64_C(0x2492492492492492), UINT64_C(0x9249249249249249),
      UINT64_C(0x0)}}};

static const mzd_local_t mask_255_255_83_a[1] = {
    {{UINT64_C(0x2492492492492492), UINT64_C(0x9249249249249249), UINT64_C(0x4924924924924924),
      UINT64_C(0x2492492492492492)}}};
static const mzd_local_t mask_255_255_83_b[1] = {
    {{UINT64_C(0x4924924924924924), UINT64_C(0x2492492492492492), UINT64_C(0x9249249249249249),
      UINT64_C(0x4924924924924924)}}};
static const mzd_local_t mask_255_255_83_c[1] = {
    {{UINT64_C(0x9249249249249248), UINT64_C(0x4924924924924924), UINT64_C(0x2492492492492492),
      UINT64_C(0x9249249249249249)}}};


static void mpc_sbox_uint64_42(mzd_local_t* out, const mzd_local_t* in, view_t* view, const rvec_t* rvec) {
  bitsliced_step_1(SC_PROOF, mzd_and_uint64_128, mzd_shift_left_uint64_128, mask_126_126_42_a, mask_126_126_42_b, mask_126_126_42_c);

  // mpc_clear(view->s, SC_PROOF);
  // a & b
  mpc_and_uint64_128(r0m, x0s, x1s, r2m, view, 0);
  // b & c
  mpc_and_uint64_128(r2m, x1s, x2m, r1s, view, 1);
  // c & a
  mpc_and_uint64_128(r1m, x0s, x2m, r0s, view, 2);

  bitsliced_step_2(SC_PROOF, mzd_xor_uint64_128, mzd_shift_right_uint64_128);
}

static void mpc_sbox_verify_uint64_42(mzd_local_t* out, const mzd_local_t* in, view_t* view, const rvec_t* rvec) {
  bitsliced_step_1(SC_VERIFY, mzd_and_uint64_128, mzd_shift_left_uint64_128, mask_126_126_42_a, mask_126_126_42_b, mask_126_126_42_c);

  // mzd_local_clear(view->s[0]);
  // a & b
  mpc_and_verify_uint64_128(r0m, x0s, x1s, r2m, view, mask_126_126_42_c, 0);
  // b & c
  mpc_and_verify_uint64_128(r2m, x1s, x2m, r1s, view, mask_126_126_42_c, 1);
  // c & a
  mpc_and_verify_uint64_128(r1m, x0s, x2m, r0s, view, mask_126_126_42_c, 2);

  bitsliced_step_2(SC_VERIFY, mzd_xor_uint64_128, mzd_shift_right_uint64_128);
}

static void mpc_and_uint64(uint64_t* res, uint64_t const* first, uint64_t const* second,
                           uint64_t const* r, view_t* view, unsigned viewshift) {
  for (unsigned m = 0; m < SC_PROOF; ++m) {
    const unsigned j = (m + 1) % SC_PROOF;
    uint64_t tmp1    = second[m] ^ second[j];
    uint64_t tmp2    = first[j] & second[m];
    tmp1             = tmp1 & first[m];
    tmp1             = tmp1 ^ tmp2;
    tmp2             = r[m] ^ r[j];
    res[m] = tmp1 = tmp1 ^ tmp2;
    if (viewshift) {
      tmp1       = tmp1 >> viewshift;
      view->t[m] = view->t[m] ^ tmp1;
    } else {
      // on first call (viewshift == 0), view->t[0..2] == 0
      view->t[m] = tmp1;
    }
  }
}

static void mpc_and_verify_uint64(uint64_t* res, uint64_t const* first, uint64_t const* second,
                                  uint64_t const* r, view_t* view, uint64_t const mask,
                                  unsigned viewshift) {
  for (unsigned m = 0; m < (SC_VERIFY - 1); ++m) {
    const unsigned j = (m + 1);
    uint64_t tmp1    = second[m] ^ second[j];
    uint64_t tmp2    = first[j] & second[m];
    tmp1             = tmp1 & first[m];
    tmp1             = tmp1 ^ tmp2;
    tmp2             = r[m] ^ r[j];
    res[m] = tmp1 = tmp1 ^ tmp2;
    if (viewshift || m) {
      tmp1       = tmp1 >> viewshift;
      view->t[m] = view->t[m] ^ tmp1;
    } else {
      // on first call (viewshift == 0), view->t[0] == 0
      view->t[m] = tmp1;
    }
  }

  const uint64_t rsc = view->t[SC_VERIFY - 1] << viewshift;
  res[SC_VERIFY - 1] = rsc & mask;
}


static void mpc_sbox_layer_bitsliced_uint64_10(uint64_t* in, view_t* view, uint64_t const* rvec) {
  //bitsliced_step_1_uint64_10(SC_PROOF);

  //mpc_and_uint64(r0m, x0s, x1s, r2m, view, 0);
  //mpc_and_uint64(r2m, x1s, x2m, r1s, view, 1);
  //mpc_and_uint64(r1m, x0s, x2m, r0s, view, 2);

  //bitsliced_step_2_uint64_10(SC_PROOF - 1);
}

static void mpc_sbox_layer_bitsliced_verify_uint64_10(uint64_t* in, view_t* view,
                                                      uint64_t const* rvec) {
  //bitsliced_step_1_uint64_10(SC_VERIFY);

  //mpc_and_verify_uint64(r0m, x0s, x1s, r2m, view, MASK_X2I, 0);
  //mpc_and_verify_uint64(r2m, x1s, x2m, r1s, view, MASK_X2I, 1);
  //mpc_and_verify_uint64(r1m, x0s, x2m, r0s, view, MASK_X2I, 2);

  //bitsliced_step_2_uint64_10(SC_VERIFY);
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

#define SBOX_uint64(sbox, y, x, views, r, n, shares, shares2)                                      \
  do {                                                                                             \
    uint64_t in[shares];                                                                           \
    for (unsigned int count = 0; count < shares; ++count) {                                        \
      in[count] = CONST_BLOCK(x[count], 0)->w64[(n) / (sizeof(word) * 8) - 1];                     \
    }                                                                                              \
    sbox(in, views, r->t);                                                                            \
    for (unsigned int count = 0; count < shares2; ++count) {                                       \
      memcpy(BLOCK(y[count], 0)->w64, CONST_BLOCK(x[count], 0)->w64,                               \
             ((n) / (sizeof(word) * 8) - 1) * sizeof(word));                                       \
      BLOCK(y[count], 0)->w64[(n) / (sizeof(word) * 8) - 1] = in[count];                           \
    }                                                                                              \
  } while (0)

#define SBOX_uint64_126_126_42(sbox, y, x, views, rvec, n, shares, shares2)                        \
  {                                                                                                \
    mzd_local_t tmp[shares];                                                                       \
    for (unsigned int count = 0; count < shares; ++count) {                                        \
      memcpy(tmp[count].w64, CONST_BLOCK(x[count], 0)->w64, sizeof(mzd_local_t));                  \
    }                                                                                              \
    sbox(tmp, tmp, views, rvec);                                                                   \
    for (unsigned int count = 0; count < shares; ++count) {                                        \
      memcpy(BLOCK(y[count], 0)->w64, tmp[count].w64, sizeof(mzd_local_t));                        \
    }                                                                                              \
  }                                                                                                \
  while (0)

#if !defined(NO_UINT64_FALLBACK)
#define SBOX SBOX_uint64_126_126_42
#define SBOX_SIGN mpc_sbox_uint64_42
#define SBOX_VERIFY mpc_sbox_verify_uint64_42

// uint64 based implementation
#include "lowmc_fns_uint64_L1.h"
#define SIGN mpc_lowmc_prove_uint64_126_126
#define VERIFY mpc_lowmc_verify_uint64_126_126
#include "mpc_lowmc.c.i"

#undef SBOX
#undef SBOX_SIGN
#undef SBOX_VERIFY
#define SBOX SBOX_uint64
#define SBOX_SIGN mpc_sbox_layer_bitsliced_uint64_10
#define SBOX_VERIFY mpc_sbox_layer_bitsliced_verify_uint64_10

#include "lowmc_fns_uint64_L3.h"
#define SIGN mpc_lowmc_call_uint64_192
#define VERIFY mpc_lowmc_call_verify_uint64_192
#include "mpc_lowmc.c.i"

#include "lowmc_fns_uint64_L5.h"
#define SIGN mpc_lowmc_call_uint64_256
#define VERIFY mpc_lowmc_call_verify_uint64_256
#include "mpc_lowmc.c.i"
#endif

#if defined(WITH_OPT)
#if defined(WITH_SSE2) || defined(WITH_NEON)
#if defined(WITH_SSE2)
#define FN_ATTR ATTR_TARGET_SSE2
#endif

#undef SBOX
#undef SBOX_SIGN
#undef SBOX_VERIFY
#define SBOX SBOX_uint64_126_126_42
#define SBOX_SIGN mpc_sbox_uint64_42
#define SBOX_VERIFY mpc_sbox_verify_uint64_42

// L1 using SSE2/NEON
#include "lowmc_fns_s128_L1.h"
#define SIGN mpc_lowmc_prove_s128_126_126
#define VERIFY mpc_lowmc_verify_s128_126_126
#include "mpc_lowmc.c.i"

#undef SBOX
#undef SBOX_SIGN
#undef SBOX_VERIFY
#define SBOX SBOX_uint64
#define SBOX_SIGN mpc_sbox_layer_bitsliced_uint64_10
#define SBOX_VERIFY mpc_sbox_layer_bitsliced_verify_uint64_10

// L3 using SSE2/NEON
#include "lowmc_fns_s128_L3.h"
#define SIGN mpc_lowmc_call_s128_192
#define VERIFY mpc_lowmc_call_verify_s128_192
#include "mpc_lowmc.c.i"

// L5 using SSE2/NEON
#include "lowmc_fns_s128_L5.h"
#define SIGN mpc_lowmc_call_s128_256
#define VERIFY mpc_lowmc_call_verify_s128_256
#include "mpc_lowmc.c.i"

#undef FN_ATTR
#endif

#if defined(WITH_AVX2)
#define FN_ATTR ATTR_TARGET_AVX2

#undef SBOX
#undef SBOX_SIGN
#undef SBOX_VERIFY
#define SBOX SBOX_uint64_126_126_42
#define SBOX_SIGN mpc_sbox_uint64_42
#define SBOX_VERIFY mpc_sbox_verify_uint64_42

// L1 using AVX2
#include "lowmc_fns_s256_L1.h"
#define SIGN mpc_lowmc_prove_s256_126_126
#define VERIFY mpc_lowmc_verify_s256_126_126
#include "mpc_lowmc.c.i"

#undef SBOX
#undef SBOX_SIGN
#undef SBOX_VERIFY
#define SBOX SBOX_uint64
#define SBOX_SIGN mpc_sbox_layer_bitsliced_uint64_10
#define SBOX_VERIFY mpc_sbox_layer_bitsliced_verify_uint64_10


// L3 using AVX2
#include "lowmc_fns_s256_L3.h"
#define SIGN mpc_lowmc_call_s256_192
#define VERIFY mpc_lowmc_call_verify_s256_192
#include "mpc_lowmc.c.i"

// L5 using AVX2
#include "lowmc_fns_s256_L5.h"
#define SIGN mpc_lowmc_call_s256_256
#define VERIFY mpc_lowmc_call_verify_s256_256
#include "mpc_lowmc.c.i"

#undef FN_ATTR

#endif
#endif

zkbpp_lowmc_implementation_f get_zkbpp_lowmc_implementation(const lowmc_t* lowmc) {
  // ASSUME(lowmc->m == 10);
  // ASSUME(lowmc->n == 128 || lowmc->n == 192 || lowmc->n == 256);

#if defined(WITH_OPT)
#if defined(WITH_AVX2)
  if (CPU_SUPPORTS_AVX2) {
    if (lowmc->n == 126 && lowmc->m == 42) {
#if defined(WITH_LOWMC_126_126_4)
      return mpc_lowmc_prove_s256_126_126_42;
#endif
    }

    /*
    if (lowmc->m == 10) {
      switch (lowmc->n) {
#if defined(WITH_LOWMC_126_126_4)
      case 128:
        return mpc_lowmc_call_s256_128_10;
#endif
#if defined(WITH_LOWMC_192_192_4)
      case 192:
        return mpc_lowmc_call_s256_192_10;
#endif
#if defined(WITH_LOWMC_255_255_4)
      case 256:
        return mpc_lowmc_call_s256_256_10;
#endif
      }
    }
    */
  }
#endif
#if defined(WITH_SSE2) || defined(WITH_NEON)
  if (CPU_SUPPORTS_SSE2 || CPU_SUPPORTS_NEON) {
    if (lowmc->n == 126 && lowmc->m == 42) {
#if defined(WITH_LOWMC_126_126_4)
      return mpc_lowmc_prove_s128_126_126_42;
#endif
    }

    /*
    if (lowmc->m == 10) {
      switch (lowmc->n) {
#if defined(WITH_LOWMC_126_126_4)
      case 128:
        return mpc_lowmc_call_s128_128_10;
#endif
#if defined(WITH_LOWMC_192_192_4)
      case 192:
        return mpc_lowmc_call_s128_192_10;
#endif
#if defined(WITH_LOWMC_255_255_4)
      case 256:
        return mpc_lowmc_call_s128_256_10;
#endif
      }
    }
    */
  }
#endif
#endif


#if !defined(NO_UINT64_FALLBACK)
  if (lowmc->n == 126 && lowmc->m == 42) {
#if defined(WITH_LOWMC_126_126_4)
    return mpc_lowmc_prove_uint64_126_126_42;
#endif
  }

  /*
  if (lowmc->m == 10) {
    switch (lowmc->n) {
#if defined(WITH_LOWMC_126_126_4)
    case 128:
      return mpc_lowmc_call_uint64_128_10;
#endif
#if defined(WITH_LOWMC_192_192_4)
    case 192:
      return mpc_lowmc_call_uint64_192_10;
#endif
#if defined(WITH_LOWMC_255_255_4)
    case 256:
      return mpc_lowmc_call_uint64_256_10;
#endif
    }
  }
  */

#endif

  return NULL;
}

zkbpp_lowmc_verify_implementation_f get_zkbpp_lowmc_verify_implementation(const lowmc_t* lowmc) {
  // ASSUME(lowmc->m == 10);
  // ASSUME(lowmc->n == 128 || lowmc->n == 192 || lowmc->n == 256);

#if defined(WITH_OPT)
#if defined(WITH_AVX2)
  if (CPU_SUPPORTS_AVX2) {
    if (lowmc->n == 126 && lowmc->m == 42) {
#if defined(WITH_LOWMC_126_126_4)
      return mpc_lowmc_verify_s128_126_126_42;
#endif
    }

/*
    if (lowmc->m == 10) {
      switch (lowmc->n) {
#if defined(WITH_LOWMC_126_126_4)
      case 128:
        return mpc_lowmc_call_verify_s256_128_10;
#endif
#if defined(WITH_LOWMC_192_192_4)
      case 192:
        return mpc_lowmc_call_verify_s256_192_10;
#endif
#if defined(WITH_LOWMC_255_255_4)
      case 256:
        return mpc_lowmc_call_verify_s256_256_10;
#endif
      }
    }
    */
  }
#endif
#if defined(WITH_SSE2) || defined(WITH_NEON)
  if (CPU_SUPPORTS_SSE2 || CPU_SUPPORTS_NEON) {
    if (lowmc->n == 126 && lowmc->m == 42) {
#if defined(WITH_LOWMC_126_126_4)
      return mpc_lowmc_verify_s128_126_126_42;
#endif
    }

    /*
    if (lowmc->m == 10) {
      switch (lowmc->n) {
#if defined(WITH_LOWMC_126_126_4)
      case 128:
        return mpc_lowmc_call_verify_s128_128_10;
#endif
#if defined(WITH_LOWMC_192_192_4)
      case 192:
        return mpc_lowmc_call_verify_s128_192_10;
#endif
#if defined(WITH_LOWMC_255_255_4)
      case 256:
        return mpc_lowmc_call_verify_s128_256_10;
#endif
      }
    }
    */
  }
#endif
#endif

#if !defined(NO_UINT64_FALLBACK)
  if (lowmc->n == 126 && lowmc->m == 42) {
    return mpc_lowmc_verify_uint64_126_126_42;
  }

  /*
  if (lowmc->m == 10) {
    switch (lowmc->n) {
#if defined(WITH_LOWMC_126_126_4)
    case 128:
      return mpc_lowmc_call_verify_uint64_128_10;
#endif
#if defined(WITH_LOWMC_192_192_4)
    case 192:
      return mpc_lowmc_call_verify_uint64_192_10;
#endif
#if defined(WITH_LOWMC_255_255_4)
    case 256:
      return mpc_lowmc_call_verify_uint64_256_10;
#endif
    }
  }
  */
#endif

  return NULL;
}

#if !defined(NO_UINT64_FALLBACK)
static void mzd_share_uint64_128(mzd_local_t* r, const mzd_local_t* v1, const mzd_local_t* v2,
                                 const mzd_local_t* v3) {
  mzd_xor_uint64_128(r, v1, v2);
  mzd_xor_uint64_128(r, r, v3);
}

static void mzd_share_uint64_192(mzd_local_t* r, const mzd_local_t* v1, const mzd_local_t* v2,
                                 const mzd_local_t* v3) {
  mzd_xor_uint64_192(r, v1, v2);
  mzd_xor_uint64_192(r, r, v3);
}

static void mzd_share_uint64_256(mzd_local_t* r, const mzd_local_t* v1, const mzd_local_t* v2,
                                 const mzd_local_t* v3) {
  mzd_xor_uint64_256(r, v1, v2);
  mzd_xor_uint64_256(r, r, v3);
}
#endif

#if defined(WITH_OPT)
#if defined(WITH_SSE2) || defined(WITH_NEON)
#if defined(WITH_SSE2)
#define FN_ATTR ATTR_TARGET_SSE2
#else
#define FN_ATTR
#endif

FN_ATTR
static void mzd_share_s128_128(mzd_local_t* r, const mzd_local_t* v1, const mzd_local_t* v2,
                               const mzd_local_t* v3) {
  mzd_xor_s128_128(r, v1, v2);
  mzd_xor_s128_128(r, r, v3);
}

FN_ATTR
static void mzd_share_s128_256(mzd_local_t* r, const mzd_local_t* v1, const mzd_local_t* v2,
                               const mzd_local_t* v3) {
  mzd_xor_s128_256(r, v1, v2);
  mzd_xor_s128_256(r, r, v3);
}

#undef FN_ATTR
#endif

#if defined(WITH_AVX2)
ATTR_TARGET_AVX2
static void mzd_share_s256_128(mzd_local_t* r, const mzd_local_t* v1, const mzd_local_t* v2,
                               const mzd_local_t* v3) {
  mzd_xor_s256_128(r, v1, v2);
  mzd_xor_s256_128(r, r, v3);
}

ATTR_TARGET_AVX2
static void mzd_share_s256_256(mzd_local_t* r, const mzd_local_t* v1, const mzd_local_t* v2,
                               const mzd_local_t* v3) {
  mzd_xor_s256_256(r, v1, v2);
  mzd_xor_s256_256(r, r, v3);
}
#endif
#endif

zkbpp_share_implementation_f get_zkbpp_share_implentation(const lowmc_t* lowmc) {
#if defined(WITH_OPT)
#if defined(WITH_AVX2)
  if (CPU_SUPPORTS_AVX2) {
    if (lowmc->n <= 128) {
      return mzd_share_s256_128;
    } else {
      return mzd_share_s256_256;
    }
  }
#endif
#if defined(WITH_SSE2) || defined(WITH_NEON)
  if (CPU_SUPPORTS_SSE2 || CPU_SUPPORTS_NEON) {
    if (lowmc->n <= 128) {
      return mzd_share_s128_128;
    } else {
      return mzd_share_s128_256;
    }
  }
#endif
#endif

#if !defined(NO_UINT64_FALLBACK)
  if (lowmc->n <= 128) {
    return mzd_share_uint64_128;
  } else if (lowmc->n <= 192) {
    return mzd_share_uint64_192;
  } else {
    return mzd_share_uint64_256;
  }
#endif
}
