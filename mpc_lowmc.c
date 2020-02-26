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
      mzd_copy_uint64_128(&view->s[m], &res[m]);
    }
  }
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
      mzd_copy_uint64_128(&view->s[m], &res[m]);
    }
  }

  mzd_shift_left_uint64_128(&res[SC_VERIFY - 1], &view->s[SC_VERIFY - 1], viewshift);
  mzd_and_uint64_128(&res[SC_VERIFY - 1], &res[SC_VERIFY - 1], mask);
}

static void mpc_and_uint64_192(mzd_local_t* res, const mzd_local_t* first,
                               const mzd_local_t* second, const mzd_local_t* r, view_t* view,
                               unsigned viewshift) {
  mzd_local_t tmp;

  for (unsigned m = 0; m < SC_PROOF; ++m) {
    const unsigned j = (m + 1) % SC_PROOF;

    // f[m] & s[m]
    mzd_and_uint64_192(&res[m], &first[m], &second[m]);

    // f[m + 1] & s[m]
    mzd_and_uint64_192(&tmp, &first[j], &second[m]);
    mzd_xor_uint64_192(&res[m], &res[m], &tmp);

    // f[m] & s[m + 1]
    mzd_and_uint64_192(&tmp, &first[m], &second[j]);
    mzd_xor_uint64_192(&res[m], &res[m], &tmp);

    // ... ^ r[m] ^ r[m + 1]
    mzd_xor_uint64_192(&res[m], &res[m], &r[m]);
    mzd_xor_uint64_192(&res[m], &res[m], &r[j]);

    if (viewshift) {
      mzd_shift_right_uint64_192(&tmp, &res[m], viewshift);
      mzd_xor_uint64_192(&view->s[m], &view->s[m], &tmp);
    } else {
      // on first call (viewshift == 0), view->t[0..2] == 0
      mzd_copy_uint64_192(&view->s[m], &res[m]);
    }
  }
}

static void mpc_and_verify_uint64_192(mzd_local_t* res, const mzd_local_t* first,
                                      const mzd_local_t* second, const mzd_local_t* r,
                                      view_t* view, const mzd_local_t* mask, unsigned viewshift) {
  mzd_local_t tmp;

  for (unsigned m = 0; m < (SC_VERIFY - 1); ++m) {
    const unsigned j = m + 1;

    mzd_and_uint64_192(&res[m], &first[m], &second[m]);

    mzd_and_uint64_192(&tmp, &first[j], &second[m]);
    mzd_xor_uint64_192(&res[m], &res[m], &tmp);

    mzd_and_uint64_192(&tmp, &first[m], &second[j]);
    mzd_xor_uint64_192(&res[m], &res[m], &tmp);

    mzd_xor_uint64_192(&res[m], &res[m], &r[m]);
    mzd_xor_uint64_192(&res[m], &res[m], &r[j]);

    if (viewshift || m) {
      mzd_shift_right_uint64_192(&tmp, &res[m], viewshift);
      mzd_xor_uint64_192(&view->s[m], &view->s[m], &tmp);
    } else {
      // on first call (viewshift == 0), view->s[0] == 0
      mzd_copy_uint64_192(&view->s[m], &res[m]);
    }
  }

  mzd_shift_left_uint64_192(&res[SC_VERIFY - 1], &view->s[SC_VERIFY - 1], viewshift);
  mzd_and_uint64_192(&res[SC_VERIFY - 1], &res[SC_VERIFY - 1], mask);
}

static void mpc_and_uint64_256(mzd_local_t* res, const mzd_local_t* first,
                               const mzd_local_t* second, const mzd_local_t* r, view_t* view,
                               unsigned viewshift) {
  mzd_local_t tmp;

  for (unsigned m = 0; m < SC_PROOF; ++m) {
    const unsigned j = (m + 1) % SC_PROOF;

    // f[m] & s[m]
    mzd_and_uint64_256(&res[m], &first[m], &second[m]);

    // f[m + 1] & s[m]
    mzd_and_uint64_256(&tmp, &first[j], &second[m]);
    mzd_xor_uint64_256(&res[m], &res[m], &tmp);

    // f[m] & s[m + 1]
    mzd_and_uint64_256(&tmp, &first[m], &second[j]);
    mzd_xor_uint64_256(&res[m], &res[m], &tmp);

    // ... ^ r[m] ^ r[m + 1]
    mzd_xor_uint64_256(&res[m], &res[m], &r[m]);
    mzd_xor_uint64_256(&res[m], &res[m], &r[j]);

    if (viewshift) {
      mzd_shift_right_uint64_256(&tmp, &res[m], viewshift);
      mzd_xor_uint64_256(&view->s[m], &view->s[m], &tmp);
    } else {
      // on first call (viewshift == 0), view->t[0..2] == 0
      mzd_copy_uint64_256(&view->s[m], &res[m]);
    }
  }
}

static void mpc_and_verify_uint64_256(mzd_local_t* res, const mzd_local_t* first,
                                      const mzd_local_t* second, const mzd_local_t* r,
                                      view_t* view, const mzd_local_t* mask, unsigned viewshift) {
  mzd_local_t tmp;

  for (unsigned m = 0; m < (SC_VERIFY - 1); ++m) {
    const unsigned j = m + 1;

    mzd_and_uint64_256(&res[m], &first[m], &second[m]);

    mzd_and_uint64_256(&tmp, &first[j], &second[m]);
    mzd_xor_uint64_256(&res[m], &res[m], &tmp);

    mzd_and_uint64_256(&tmp, &first[m], &second[j]);
    mzd_xor_uint64_256(&res[m], &res[m], &tmp);

    mzd_xor_uint64_256(&res[m], &res[m], &r[m]);
    mzd_xor_uint64_256(&res[m], &res[m], &r[j]);

    if (viewshift || m) {
      mzd_shift_right_uint64_256(&tmp, &res[m], viewshift);
      mzd_xor_uint64_256(&view->s[m], &view->s[m], &tmp);
    } else {
      // on first call (viewshift == 0), view->s[0] == 0
      mzd_copy_uint64_256(&view->s[m], &res[m]);
    }
  }

  mzd_shift_left_uint64_256(&res[SC_VERIFY - 1], &view->s[SC_VERIFY - 1], viewshift);
  mzd_and_uint64_256(&res[SC_VERIFY - 1], &res[SC_VERIFY - 1], mask);
}

static void mpc_and_s128_128(mzd_local_t* res, const mzd_local_t* first,
                               const mzd_local_t* second, const mzd_local_t* r, view_t* view,
                               unsigned viewshift) {
  mzd_local_t tmp;

  for (unsigned m = 0; m < SC_PROOF; ++m) {
    const unsigned j = (m + 1) % SC_PROOF;

    // f[m] & s[m]
    mzd_and_s128_128(&res[m], &first[m], &second[m]);

    // f[m + 1] & s[m]
    mzd_and_s128_128(&tmp, &first[j], &second[m]);
    mzd_xor_s128_128(&res[m], &res[m], &tmp);

    // f[m] & s[m + 1]
    mzd_and_s128_128(&tmp, &first[m], &second[j]);
    mzd_xor_s128_128(&res[m], &res[m], &tmp);

    // ... ^ r[m] ^ r[m + 1]
    mzd_xor_s128_128(&res[m], &res[m], &r[m]);
    mzd_xor_s128_128(&res[m], &res[m], &r[j]);

    if (viewshift) {
      mzd_shift_right_s128_128(&tmp, &res[m], viewshift);
      mzd_xor_s128_128(&view->s[m], &view->s[m], &tmp);
    } else {
      // on first call (viewshift == 0), view->t[0..2] == 0
      mzd_copy_s128_128(&view->s[m], &res[m]);
    }
  }
}

static void mpc_and_verify_s128_128(mzd_local_t* res, const mzd_local_t* first,
                                      const mzd_local_t* second, const mzd_local_t* r,
                                      view_t* view, const mzd_local_t* mask, unsigned viewshift) {
  mzd_local_t tmp;

  for (unsigned m = 0; m < (SC_VERIFY - 1); ++m) {
    const unsigned j = m + 1;

    mzd_and_s128_128(&res[m], &first[m], &second[m]);

    mzd_and_s128_128(&tmp, &first[j], &second[m]);
    mzd_xor_s128_128(&res[m], &res[m], &tmp);

    mzd_and_s128_128(&tmp, &first[m], &second[j]);
    mzd_xor_s128_128(&res[m], &res[m], &tmp);

    mzd_xor_s128_128(&res[m], &res[m], &r[m]);
    mzd_xor_s128_128(&res[m], &res[m], &r[j]);

    if (viewshift || m) {
      mzd_shift_right_s128_128(&tmp, &res[m], viewshift);
      mzd_xor_s128_128(&view->s[m], &view->s[m], &tmp);
    } else {
      // on first call (viewshift == 0), view->s[0] == 0
      mzd_copy_s128_128(&view->s[m], &res[m]);
    }
  }

  mzd_shift_left_s128_128(&res[SC_VERIFY - 1], &view->s[SC_VERIFY - 1], viewshift);
  mzd_and_s128_128(&res[SC_VERIFY - 1], &res[SC_VERIFY - 1], mask);
}

static void mpc_and_s128_256(mzd_local_t* res, const mzd_local_t* first,
                               const mzd_local_t* second, const mzd_local_t* r, view_t* view,
                               unsigned viewshift) {
  mzd_local_t tmp;

  for (unsigned m = 0; m < SC_PROOF; ++m) {
    const unsigned j = (m + 1) % SC_PROOF;

    // f[m] & s[m]
    mzd_and_s128_256(&res[m], &first[m], &second[m]);

    // f[m + 1] & s[m]
    mzd_and_s128_256(&tmp, &first[j], &second[m]);
    mzd_xor_s128_256(&res[m], &res[m], &tmp);

    // f[m] & s[m + 1]
    mzd_and_s128_256(&tmp, &first[m], &second[j]);
    mzd_xor_s128_256(&res[m], &res[m], &tmp);

    // ... ^ r[m] ^ r[m + 1]
    mzd_xor_s128_256(&res[m], &res[m], &r[m]);
    mzd_xor_s128_256(&res[m], &res[m], &r[j]);

    if (viewshift) {
      mzd_shift_right_s128_256(&tmp, &res[m], viewshift);
      mzd_xor_s128_256(&view->s[m], &view->s[m], &tmp);
    } else {
      // on first call (viewshift == 0), view->t[0..2] == 0
      mzd_copy_s128_256(&view->s[m], &res[m]);
    }
  }
}

static void mpc_and_verify_s128_256(mzd_local_t* res, const mzd_local_t* first,
                                      const mzd_local_t* second, const mzd_local_t* r,
                                      view_t* view, const mzd_local_t* mask, unsigned viewshift) {
  mzd_local_t tmp;

  for (unsigned m = 0; m < (SC_VERIFY - 1); ++m) {
    const unsigned j = m + 1;

    mzd_and_s128_256(&res[m], &first[m], &second[m]);

    mzd_and_s128_256(&tmp, &first[j], &second[m]);
    mzd_xor_s128_256(&res[m], &res[m], &tmp);

    mzd_and_s128_256(&tmp, &first[m], &second[j]);
    mzd_xor_s128_256(&res[m], &res[m], &tmp);

    mzd_xor_s128_256(&res[m], &res[m], &r[m]);
    mzd_xor_s128_256(&res[m], &res[m], &r[j]);

    if (viewshift || m) {
      mzd_shift_right_s128_256(&tmp, &res[m], viewshift);
      mzd_xor_s128_256(&view->s[m], &view->s[m], &tmp);
    } else {
      // on first call (viewshift == 0), view->s[0] == 0
      mzd_copy_s128_256(&view->s[m], &res[m]);
    }
  }

  mzd_shift_left_s128_256(&res[SC_VERIFY - 1], &view->s[SC_VERIFY - 1], viewshift);
  mzd_and_s128_256(&res[SC_VERIFY - 1], &res[SC_VERIFY - 1], mask);
}

static void mpc_and_s256_128(mzd_local_t* res, const mzd_local_t* first,
                               const mzd_local_t* second, const mzd_local_t* r, view_t* view,
                               unsigned viewshift) {
  mzd_local_t tmp;

  for (unsigned m = 0; m < SC_PROOF; ++m) {
    const unsigned j = (m + 1) % SC_PROOF;

    // f[m] & s[m]
    mzd_and_s256_128(&res[m], &first[m], &second[m]);

    // f[m + 1] & s[m]
    mzd_and_s256_128(&tmp, &first[j], &second[m]);
    mzd_xor_s256_128(&res[m], &res[m], &tmp);

    // f[m] & s[m + 1]
    mzd_and_s256_128(&tmp, &first[m], &second[j]);
    mzd_xor_s256_128(&res[m], &res[m], &tmp);

    // ... ^ r[m] ^ r[m + 1]
    mzd_xor_s256_128(&res[m], &res[m], &r[m]);
    mzd_xor_s256_128(&res[m], &res[m], &r[j]);

    if (viewshift) {
      mzd_shift_right_s256_128(&tmp, &res[m], viewshift);
      mzd_xor_s256_128(&view->s[m], &view->s[m], &tmp);
    } else {
      // on first call (viewshift == 0), view->t[0..2] == 0
      mzd_copy_s256_128(&view->s[m], &res[m]);
    }
  }
}

static void mpc_and_verify_s256_128(mzd_local_t* res, const mzd_local_t* first,
                                      const mzd_local_t* second, const mzd_local_t* r,
                                      view_t* view, const mzd_local_t* mask, unsigned viewshift) {
  mzd_local_t tmp;

  for (unsigned m = 0; m < (SC_VERIFY - 1); ++m) {
    const unsigned j = m + 1;

    mzd_and_s256_128(&res[m], &first[m], &second[m]);

    mzd_and_s256_128(&tmp, &first[j], &second[m]);
    mzd_xor_s256_128(&res[m], &res[m], &tmp);

    mzd_and_s256_128(&tmp, &first[m], &second[j]);
    mzd_xor_s256_128(&res[m], &res[m], &tmp);

    mzd_xor_s256_128(&res[m], &res[m], &r[m]);
    mzd_xor_s256_128(&res[m], &res[m], &r[j]);

    if (viewshift || m) {
      mzd_shift_right_s256_128(&tmp, &res[m], viewshift);
      mzd_xor_s256_128(&view->s[m], &view->s[m], &tmp);
    } else {
      // on first call (viewshift == 0), view->s[0] == 0
      mzd_copy_s256_128(&view->s[m], &res[m]);
    }
  }

  mzd_shift_left_s256_128(&res[SC_VERIFY - 1], &view->s[SC_VERIFY - 1], viewshift);
  mzd_and_s256_128(&res[SC_VERIFY - 1], &res[SC_VERIFY - 1], mask);
}

static void mpc_and_s256_256(mzd_local_t* res, const mzd_local_t* first,
                               const mzd_local_t* second, const mzd_local_t* r, view_t* view,
                               unsigned viewshift) {
  mzd_local_t tmp;

  for (unsigned m = 0; m < SC_PROOF; ++m) {
    const unsigned j = (m + 1) % SC_PROOF;

    // f[m] & s[m]
    mzd_and_s256_256(&res[m], &first[m], &second[m]);

    // f[m + 1] & s[m]
    mzd_and_s256_256(&tmp, &first[j], &second[m]);
    mzd_xor_s256_256(&res[m], &res[m], &tmp);

    // f[m] & s[m + 1]
    mzd_and_s256_256(&tmp, &first[m], &second[j]);
    mzd_xor_s256_256(&res[m], &res[m], &tmp);

    // ... ^ r[m] ^ r[m + 1]
    mzd_xor_s256_256(&res[m], &res[m], &r[m]);
    mzd_xor_s256_256(&res[m], &res[m], &r[j]);

    if (viewshift) {
      mzd_shift_right_s256_256(&tmp, &res[m], viewshift);
      mzd_xor_s256_256(&view->s[m], &view->s[m], &tmp);
    } else {
      // on first call (viewshift == 0), view->t[0..2] == 0
      mzd_copy_s256_256(&view->s[m], &res[m]);
    }
  }
}

static void mpc_and_verify_s256_256(mzd_local_t* res, const mzd_local_t* first,
                                      const mzd_local_t* second, const mzd_local_t* r,
                                      view_t* view, const mzd_local_t* mask, unsigned viewshift) {
  mzd_local_t tmp;

  for (unsigned m = 0; m < (SC_VERIFY - 1); ++m) {
    const unsigned j = m + 1;

    mzd_and_s256_256(&res[m], &first[m], &second[m]);

    mzd_and_s256_256(&tmp, &first[j], &second[m]);
    mzd_xor_s256_256(&res[m], &res[m], &tmp);

    mzd_and_s256_256(&tmp, &first[m], &second[j]);
    mzd_xor_s256_256(&res[m], &res[m], &tmp);

    mzd_xor_s256_256(&res[m], &res[m], &r[m]);
    mzd_xor_s256_256(&res[m], &res[m], &r[j]);

    if (viewshift || m) {
      mzd_shift_right_s256_256(&tmp, &res[m], viewshift);
      mzd_xor_s256_256(&view->s[m], &view->s[m], &tmp);
    } else {
      // on first call (viewshift == 0), view->s[0] == 0
      mzd_copy_s256_256(&view->s[m], &res[m]);
    }
  }

  mzd_shift_left_s256_256(&res[SC_VERIFY - 1], &view->s[SC_VERIFY - 1], viewshift);
  mzd_and_s256_256(&res[SC_VERIFY - 1], &res[SC_VERIFY - 1], mask);
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

static void mpc_sbox_prove_uint64_lowmc_126_126_4(mzd_local_t* out, const mzd_local_t* in, view_t* view, const rvec_t* rvec) {
  bitsliced_step_1(SC_PROOF, mzd_and_uint64_128, mzd_shift_left_uint64_128, mask_126_126_42_a, mask_126_126_42_b, mask_126_126_42_c);

  // a & b
  mpc_and_uint64_128(r0m, x0s, x1s, r2m, view, 0);
  // b & c
  mpc_and_uint64_128(r2m, x1s, x2m, r1s, view, 1);
  // c & a
  mpc_and_uint64_128(r1m, x0s, x2m, r0s, view, 2);

  bitsliced_step_2(SC_PROOF, mzd_xor_uint64_128, mzd_shift_right_uint64_128);
}

static void mpc_sbox_verify_uint64_lowmc_126_126_4(mzd_local_t* out, const mzd_local_t* in, view_t* view, const rvec_t* rvec) {
  bitsliced_step_1(SC_VERIFY, mzd_and_uint64_128, mzd_shift_left_uint64_128, mask_126_126_42_a, mask_126_126_42_b, mask_126_126_42_c);

  // a & b
  mpc_and_verify_uint64_128(r0m, x0s, x1s, r2m, view, mask_126_126_42_c, 0);
  // b & c
  mpc_and_verify_uint64_128(r2m, x1s, x2m, r1s, view, mask_126_126_42_c, 1);
  // c & a
  mpc_and_verify_uint64_128(r1m, x0s, x2m, r0s, view, mask_126_126_42_c, 2);

  bitsliced_step_2(SC_VERIFY, mzd_xor_uint64_128, mzd_shift_right_uint64_128);
}

static void mpc_sbox_prove_uint64_lowmc_192_192_4(mzd_local_t* out, const mzd_local_t* in, view_t* view, const rvec_t* rvec) {
  bitsliced_step_1(SC_PROOF, mzd_and_uint64_192, mzd_rotate_left_uint64_192, mask_192_192_64_a, mask_192_192_64_b, mask_192_192_64_c);

  // a & b
  mpc_and_uint64_192(r0m, x0s, x1s, r2m, view, 0);
  // b & c
  mpc_and_uint64_192(r2m, x1s, x2m, r1s, view, 1);
  // c & a
  mpc_and_uint64_192(r1m, x0s, x2m, r0s, view, 2);

  bitsliced_step_2(SC_PROOF, mzd_xor_uint64_192, mzd_rotate_right_uint64_192);
}

static void mpc_sbox_verify_uint64_lowmc_192_192_4(mzd_local_t* out, const mzd_local_t* in, view_t* view, const rvec_t* rvec) {
  bitsliced_step_1(SC_VERIFY, mzd_and_uint64_192, mzd_rotate_left_uint64_192, mask_192_192_64_a, mask_192_192_64_b, mask_192_192_64_c);

  // a & b
  mpc_and_verify_uint64_192(r0m, x0s, x1s, r2m, view, mask_192_192_64_c, 0);
  // b & c
  mpc_and_verify_uint64_192(r2m, x1s, x2m, r1s, view, mask_192_192_64_c, 1);
  // c & a
  mpc_and_verify_uint64_192(r1m, x0s, x2m, r0s, view, mask_192_192_64_c, 2);

  bitsliced_step_2(SC_VERIFY, mzd_xor_uint64_192, mzd_rotate_right_uint64_192);
}

static void mpc_sbox_prove_uint64_lowmc_255_255_4(mzd_local_t* out, const mzd_local_t* in, view_t* view, const rvec_t* rvec) {
  bitsliced_step_1(SC_PROOF, mzd_and_uint64_256, mzd_rotate_left_uint64_256, mask_255_255_85_a, mask_255_255_85_b, mask_255_255_85_c);

  // a & b
  mpc_and_uint64_256(r0m, x0s, x1s, r2m, view, 0);
  // b & c
  mpc_and_uint64_256(r2m, x1s, x2m, r1s, view, 1);
  // c & a
  mpc_and_uint64_256(r1m, x0s, x2m, r0s, view, 2);

  bitsliced_step_2(SC_PROOF, mzd_xor_uint64_256, mzd_rotate_right_uint64_256);
}

static void mpc_sbox_verify_uint64_lowmc_255_255_4(mzd_local_t* out, const mzd_local_t* in, view_t* view, const rvec_t* rvec) {
  bitsliced_step_1(SC_VERIFY, mzd_and_uint64_256, mzd_shift_left_uint64_256, mask_255_255_85_a, mask_255_255_85_b, mask_255_255_85_c);

  // a & b
  mpc_and_verify_uint64_256(r0m, x0s, x1s, r2m, view, mask_255_255_85_c, 0);
  // b & c
  mpc_and_verify_uint64_256(r2m, x1s, x2m, r1s, view, mask_255_255_85_c, 1);
  // c & a
  mpc_and_verify_uint64_256(r1m, x0s, x2m, r0s, view, mask_255_255_85_c, 2);

  bitsliced_step_2(SC_VERIFY, mzd_xor_uint64_256, mzd_shift_right_uint64_256);
}

static void mpc_sbox_prove_s128_lowmc_126_126_4(mzd_local_t* out, const mzd_local_t* in, view_t* view, const rvec_t* rvec) {
  bitsliced_step_1(SC_PROOF, mzd_and_s128_128, mzd_shift_left_s128_128, mask_126_126_42_a, mask_126_126_42_b, mask_126_126_42_c);

  // a & b
  mpc_and_s128_128(r0m, x0s, x1s, r2m, view, 0);
  // b & c
  mpc_and_s128_128(r2m, x1s, x2m, r1s, view, 1);
  // c & a
  mpc_and_s128_128(r1m, x0s, x2m, r0s, view, 2);

  bitsliced_step_2(SC_PROOF, mzd_xor_s128_128, mzd_shift_right_s128_128);
}

static void mpc_sbox_verify_s128_lowmc_126_126_4(mzd_local_t* out, const mzd_local_t* in, view_t* view, const rvec_t* rvec) {
  bitsliced_step_1(SC_VERIFY, mzd_and_s128_128, mzd_shift_left_s128_128, mask_126_126_42_a, mask_126_126_42_b, mask_126_126_42_c);

  // a & b
  mpc_and_verify_s128_128(r0m, x0s, x1s, r2m, view, mask_126_126_42_c, 0);
  // b & c
  mpc_and_verify_s128_128(r2m, x1s, x2m, r1s, view, mask_126_126_42_c, 1);
  // c & a
  mpc_and_verify_s128_128(r1m, x0s, x2m, r0s, view, mask_126_126_42_c, 2);

  bitsliced_step_2(SC_VERIFY, mzd_xor_s128_128, mzd_shift_right_s128_128);
}

static void mpc_sbox_prove_s128_lowmc_192_192_4(mzd_local_t* out, const mzd_local_t* in, view_t* view, const rvec_t* rvec) {
  bitsliced_step_1(SC_PROOF, mzd_and_s128_256, mzd_shift_left_s128_256, mask_192_192_64_a, mask_192_192_64_b, mask_192_192_64_c);

  // a & b
  mpc_and_s128_256(r0m, x0s, x1s, r2m, view, 0);
  // b & c
  mpc_and_s128_256(r2m, x1s, x2m, r1s, view, 1);
  // c & a
  mpc_and_s128_256(r1m, x0s, x2m, r0s, view, 2);

  bitsliced_step_2(SC_PROOF, mzd_xor_s128_256, mzd_shift_right_s128_256);
}

static void mpc_sbox_verify_s128_lowmc_192_192_4(mzd_local_t* out, const mzd_local_t* in, view_t* view, const rvec_t* rvec) {
  bitsliced_step_1(SC_VERIFY, mzd_and_s128_256, mzd_shift_left_s128_256, mask_192_192_64_a, mask_192_192_64_b, mask_192_192_64_c);

  // a & b
  mpc_and_verify_s128_256(r0m, x0s, x1s, r2m, view, mask_192_192_64_c, 0);
  // b & c
  mpc_and_verify_s128_256(r2m, x1s, x2m, r1s, view, mask_192_192_64_c, 1);
  // c & a
  mpc_and_verify_s128_256(r1m, x0s, x2m, r0s, view, mask_192_192_64_c, 2);

  bitsliced_step_2(SC_VERIFY, mzd_xor_s128_256, mzd_shift_right_s128_256);
}

static void mpc_sbox_prove_s128_lowmc_255_255_4(mzd_local_t* out, const mzd_local_t* in, view_t* view, const rvec_t* rvec) {
  bitsliced_step_1(SC_PROOF, mzd_and_s128_256, mzd_rotate_left_s128_256, mask_255_255_85_a, mask_255_255_85_b, mask_255_255_85_c);

  // a & b
  mpc_and_s128_256(r0m, x0s, x1s, r2m, view, 0);
  // b & c
  mpc_and_s128_256(r2m, x1s, x2m, r1s, view, 1);
  // c & a
  mpc_and_s128_256(r1m, x0s, x2m, r0s, view, 2);

  bitsliced_step_2(SC_PROOF, mzd_xor_s128_256, mzd_rotate_right_s128_256);
}

static void mpc_sbox_verify_s128_lowmc_255_255_4(mzd_local_t* out, const mzd_local_t* in, view_t* view, const rvec_t* rvec) {
  bitsliced_step_1(SC_VERIFY, mzd_and_s128_256, mzd_rotate_left_s128_256, mask_255_255_85_a, mask_255_255_85_b, mask_255_255_85_c);

  // a & b
  mpc_and_verify_s128_256(r0m, x0s, x1s, r2m, view, mask_255_255_85_c, 0);
  // b & c
  mpc_and_verify_s128_256(r2m, x1s, x2m, r1s, view, mask_255_255_85_c, 1);
  // c & a
  mpc_and_verify_s128_256(r1m, x0s, x2m, r0s, view, mask_255_255_85_c, 2);

  bitsliced_step_2(SC_VERIFY, mzd_xor_s128_256, mzd_rotate_right_s128_256);
}

static void mpc_sbox_prove_s256_lowmc_126_126_4(mzd_local_t* out, const mzd_local_t* in, view_t* view, const rvec_t* rvec) {
  bitsliced_step_1(SC_PROOF, mzd_and_s256_128, mzd_shift_left_s256_128, mask_126_126_42_a, mask_126_126_42_b, mask_126_126_42_c);

  // a & b
  mpc_and_s256_128(r0m, x0s, x1s, r2m, view, 0);
  // b & c
  mpc_and_s256_128(r2m, x1s, x2m, r1s, view, 1);
  // c & a
  mpc_and_s256_128(r1m, x0s, x2m, r0s, view, 2);

  bitsliced_step_2(SC_PROOF, mzd_xor_s256_128, mzd_shift_right_s256_128);
}

static void mpc_sbox_verify_s256_lowmc_126_126_4(mzd_local_t* out, const mzd_local_t* in, view_t* view, const rvec_t* rvec) {
  bitsliced_step_1(SC_VERIFY, mzd_and_s256_128, mzd_shift_left_s256_128, mask_126_126_42_a, mask_126_126_42_b, mask_126_126_42_c);

  // a & b
  mpc_and_verify_s256_128(r0m, x0s, x1s, r2m, view, mask_126_126_42_c, 0);
  // b & c
  mpc_and_verify_s256_128(r2m, x1s, x2m, r1s, view, mask_126_126_42_c, 1);
  // c & a
  mpc_and_verify_s256_128(r1m, x0s, x2m, r0s, view, mask_126_126_42_c, 2);

  bitsliced_step_2(SC_VERIFY, mzd_xor_s256_128, mzd_shift_right_s256_128);
}

static void mpc_sbox_prove_s256_lowmc_192_192_4(mzd_local_t* out, const mzd_local_t* in, view_t* view, const rvec_t* rvec) {
  bitsliced_step_1(SC_PROOF, mzd_and_s256_256, mzd_shift_left_s256_256, mask_192_192_64_a, mask_192_192_64_b, mask_192_192_64_c);

  // a & b
  mpc_and_s256_256(r0m, x0s, x1s, r2m, view, 0);
  // b & c
  mpc_and_s256_256(r2m, x1s, x2m, r1s, view, 1);
  // c & a
  mpc_and_s256_256(r1m, x0s, x2m, r0s, view, 2);

  bitsliced_step_2(SC_PROOF, mzd_xor_s256_256, mzd_shift_right_s256_256);
}

static void mpc_sbox_verify_s256_lowmc_192_192_4(mzd_local_t* out, const mzd_local_t* in, view_t* view, const rvec_t* rvec) {
  bitsliced_step_1(SC_VERIFY, mzd_and_s256_256, mzd_shift_left_s256_256, mask_192_192_64_a, mask_192_192_64_b, mask_192_192_64_c);

  // a & b
  mpc_and_verify_s256_256(r0m, x0s, x1s, r2m, view, mask_192_192_64_c, 0);
  // b & c
  mpc_and_verify_s256_256(r2m, x1s, x2m, r1s, view, mask_192_192_64_c, 1);
  // c & a
  mpc_and_verify_s256_256(r1m, x0s, x2m, r0s, view, mask_192_192_64_c, 2);

  bitsliced_step_2(SC_VERIFY, mzd_xor_s256_256, mzd_shift_right_s256_256);
}

static void mpc_sbox_prove_s256_lowmc_255_255_4(mzd_local_t* out, const mzd_local_t* in, view_t* view, const rvec_t* rvec) {
  bitsliced_step_1(SC_PROOF, mzd_and_s256_256, mzd_rotate_left_s256_256, mask_255_255_85_a, mask_255_255_85_b, mask_255_255_85_c);

  // a & b
  mpc_and_s256_256(r0m, x0s, x1s, r2m, view, 0);
  // b & c
  mpc_and_s256_256(r2m, x1s, x2m, r1s, view, 1);
  // c & a
  mpc_and_s256_256(r1m, x0s, x2m, r0s, view, 2);

  bitsliced_step_2(SC_PROOF, mzd_xor_s256_256, mzd_rotate_right_s256_256);
}

static void mpc_sbox_verify_s256_lowmc_255_255_4(mzd_local_t* out, const mzd_local_t* in, view_t* view, const rvec_t* rvec) {
  bitsliced_step_1(SC_VERIFY, mzd_and_s256_256, mzd_rotate_left_s256_256, mask_255_255_85_a, mask_255_255_85_b, mask_255_255_85_c);

  // a & b
  mpc_and_verify_s256_256(r0m, x0s, x1s, r2m, view, mask_255_255_85_c, 0);
  // b & c
  mpc_and_verify_s256_256(r2m, x1s, x2m, r1s, view, mask_255_255_85_c, 1);
  // c & a
  mpc_and_verify_s256_256(r1m, x0s, x2m, r0s, view, mask_255_255_85_c, 2);

  bitsliced_step_2(SC_VERIFY, mzd_xor_s256_256, mzd_rotate_right_s256_256);
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

/* TODO: get rid of the copies */
#define SBOX_mzd(sbox, y, x, views, rvec, n, shares, shares2)                                      \
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

#define SBOX SBOX_mzd

#if !defined(NO_UINT64_FALLBACK)
#define IMPL uint64

// uint64 based implementation
#include "lowmc_fns_uint64_L1.h"
#include "mpc_lowmc.c.i"

#include "lowmc_fns_uint64_L3.h"
#include "mpc_lowmc.c.i"

#include "lowmc_fns_uint64_L5.h"
#include "mpc_lowmc.c.i"
#endif

#if defined(WITH_OPT)
#if defined(WITH_SSE2) || defined(WITH_NEON)
#if defined(WITH_SSE2)
#define FN_ATTR ATTR_TARGET_SSE2
#endif

#undef IMPL
#define IMPL s128

// L1 using SSE2/NEON
#include "lowmc_fns_s128_L1.h"
#include "mpc_lowmc.c.i"

// L3 using SSE2/NEON
#include "lowmc_fns_s128_L3.h"
#include "mpc_lowmc.c.i"

// L5 using SSE2/NEON
#include "lowmc_fns_s128_L5.h"
#include "mpc_lowmc.c.i"

#undef FN_ATTR
#endif

#if defined(WITH_AVX2)
#define FN_ATTR ATTR_TARGET_AVX2

#undef IMPL
#define IMPL s256

// L1 using AVX2
#include "lowmc_fns_s256_L1.h"
#include "mpc_lowmc.c.i"

// L3 using AVX2
#include "lowmc_fns_s256_L3.h"
#include "mpc_lowmc.c.i"

// L5 using AVX2
#include "lowmc_fns_s256_L5.h"
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
      return mpc_lowmc_prove_s256_lowmc_126_126_4;
#endif
    }
    if (lowmc->n == 192 && lowmc->m == 64) {
#if defined(WITH_LOWMC_192_192_4)
      return mpc_lowmc_prove_s256_lowmc_192_192_4;
#endif
    }
    if (lowmc->n == 255 && lowmc->m == 85) {
#if defined(WITH_LOWMC_255_255_4)
      return mpc_lowmc_prove_s256_lowmc_255_255_4;
#endif
    }
  }
#endif
#if defined(WITH_SSE2) || defined(WITH_NEON)
  if (CPU_SUPPORTS_SSE2 || CPU_SUPPORTS_NEON) {
    if (lowmc->n == 126 && lowmc->m == 42) {
#if defined(WITH_LOWMC_126_126_4)
      return mpc_lowmc_prove_s128_lowmc_126_126_4;
#endif
    }
    if (lowmc->n == 192 && lowmc->m == 64) {
#if defined(WITH_LOWMC_192_192_4)
      return mpc_lowmc_prove_s128_lowmc_192_192_4;
#endif
    }
    if (lowmc->n == 255 && lowmc->m == 85) {
#if defined(WITH_LOWMC_255_255_4)
      return mpc_lowmc_prove_s128_lowmc_255_255_4;
#endif
    }
  }
#endif
#endif


#if !defined(NO_UINT64_FALLBACK)
  if (lowmc->n == 126 && lowmc->m == 42) {
#if defined(WITH_LOWMC_126_126_4)
    return mpc_lowmc_prove_uint64_lowmc_126_126_4;
#endif
  }
  if (lowmc->n == 192 && lowmc->m == 64) {
#if defined(WITH_LOWMC_192_192_4)
    return mpc_lowmc_prove_uint64_lowmc_192_192_4;
#endif
  }
  if (lowmc->n == 255 && lowmc->m == 85) {
#if defined(WITH_LOWMC_255_255_4)
    return mpc_lowmc_prove_uint64_lowmc_255_255_4;
#endif
  }
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
      return mpc_lowmc_verify_s256_lowmc_126_126_4;
#endif
    }
    if (lowmc->n == 192 && lowmc->m == 64) {
#if defined(WITH_LOWMC_192_192_4)
      return mpc_lowmc_verify_s256_lowmc_192_192_4;
#endif
    }
    if (lowmc->n == 255 && lowmc->m == 85) {
#if defined(WITH_LOWMC_255_255_4)
      return mpc_lowmc_verify_s256_lowmc_255_255_4;
#endif
    }
  }
#endif
#if defined(WITH_SSE2) || defined(WITH_NEON)
  if (CPU_SUPPORTS_SSE2 || CPU_SUPPORTS_NEON) {
    if (lowmc->n == 126 && lowmc->m == 42) {
#if defined(WITH_LOWMC_126_126_4)
      return mpc_lowmc_verify_s128_lowmc_126_126_4;
#endif
    }
    if (lowmc->n == 192 && lowmc->m == 64) {
#if defined(WITH_LOWMC_192_192_4)
      return mpc_lowmc_verify_s128_lowmc_192_192_4;
#endif
    }
    if (lowmc->n == 255 && lowmc->m == 85) {
#if defined(WITH_LOWMC_255_255_4)
      return mpc_lowmc_verify_s128_lowmc_255_255_4;
#endif
    }
  }
#endif
#endif


#if !defined(NO_UINT64_FALLBACK)
  if (lowmc->n == 126 && lowmc->m == 42) {
#if defined(WITH_LOWMC_126_126_4)
    return mpc_lowmc_verify_uint64_lowmc_126_126_4;
#endif
  }
  if (lowmc->n == 192 && lowmc->m == 64) {
#if defined(WITH_LOWMC_192_192_4)
    return mpc_lowmc_verify_uint64_lowmc_192_192_4;
#endif
  }
  if (lowmc->n == 255 && lowmc->m == 85) {
#if defined(WITH_LOWMC_255_255_4)
    return mpc_lowmc_verify_uint64_lowmc_255_255_4;
#endif
  }
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
