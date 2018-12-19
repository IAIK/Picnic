/*
 *  This file is part of the optimized implementation of the Picnic signature scheme.
 *  See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */

#include "../mzd_additional.h"

#ifdef WITH_OPT
#include "../simd.h"
#if !defined(_MSC_VER)
#include <stdalign.h>
#endif
#endif

#include <m4ri/m4ri.h>

#include "utils.h"

static int test_mzd_local_equal(void) {
  int ret = 0;

  for (unsigned int i = 0; i < 10; ++i) {
    mzd_local_t* a = mzd_local_init(1, (i + 1) * 64);
    mzd_randomize_ssl(a);
    mzd_local_t* b = mzd_local_init(1, (i + 1) * 64);
    mzd_local_copy(b, a);

    if (mzd_local_equal(a, b)) {
      printf("equal: ok [%u]\n", (i + 1) * 64);
    } else {
      printf("equal: fail [%u]\n", (i + 1) * 64);
      ret = -1;
    }

    mzd_xor_uint64(b, b, a);
    if (!mzd_local_equal(a, b)) {
      printf("equal: ok [%u]\n", (i + 1) * 64);
    } else {
      printf("equal: fail [%u]\n", (i + 1) * 64);
      ret = -1;
    }

    mzd_local_free(a);
    mzd_local_free(b);
  }

  return ret;
}

typedef void (*mul_fn)(mzd_local_t*, const mzd_local_t*, const mzd_local_t*);

static int test_mzd_mul_f(const char* n, unsigned int rows, unsigned int cols, mul_fn f,
                          bool is_addmul) {
  int ret = 0;

  mzd_t* A = mzd_init(rows, cols);
  mzd_t* v = mzd_init(1, rows);
  mzd_t* c = mzd_init(1, cols);

  mzd_randomize(A);
  mzd_randomize(v);
  mzd_randomize(c);

  mzd_local_t* Al = mzd_convert(A);
  mzd_local_t* vl = mzd_convert(v);
  mzd_local_t* c2 = mzd_convert(c);

  for (unsigned int k = 0; k < 3; ++k) {
    mzd_t* r = is_addmul ? mzd_addmul_naive(c, v, A) : mzd_mul_naive(c, v, A);
    f(c2, vl, Al);

    mzd_local_t* rc = mzd_convert(r);

    if (!mzd_local_equal(rc, c2)) {
      printf("%s: fail [%u x %u]\n", n, rows, cols);
      ret = -1;
    } else {
      printf("%s: ok [%u x %u]\n", n, rows, cols);
    }

    mzd_local_free(rc);
  }

  mzd_local_free(c2);
  mzd_local_free(vl);
  mzd_local_free(Al);

  mzd_free(c);
  mzd_free(v);
  mzd_free(A);

  return ret;
}

#if defined(MUL_M4RI)
static int test_mzd_mul_l_f(const char* n, unsigned int rows, unsigned int cols, mul_fn f,
                            bool is_addmul) {
  int ret = 0;

  mzd_t* A = mzd_init(rows, cols);
  mzd_t* v = mzd_init(1, rows);
  mzd_t* c = mzd_init(1, cols);

  mzd_randomize(A);
  mzd_randomize(v);
  mzd_randomize(c);

  mzd_local_t* Al  = mzd_convert(A);
  mzd_local_t* All = mzd_precompute_matrix_lookup(Al);
  mzd_local_t* vl  = mzd_convert(v);
  mzd_local_t* c2  = mzd_convert(c);

  for (unsigned int k = 0; k < 3; ++k) {
    mzd_t* r = is_addmul ? mzd_addmul_naive(c, v, A) : mzd_mul_naive(c, v, A);
    f(c2, vl, All);

    mzd_local_t* rc = mzd_convert(r);

    if (!mzd_local_equal(rc, c2)) {
      printf("%s: fail [%u x %u]\n", n, rows, cols);
      ret = -1;
    } else {
      printf("%s: ok [%u x %u]\n", n, rows, cols);
    }

    mzd_local_free(rc);
  }

  mzd_local_free(c2);
  mzd_local_free(vl);
  mzd_local_free(All);
  mzd_local_free(Al);

  mzd_free(c);
  mzd_free(v);
  mzd_free(A);

  return ret;
}
#endif

static int test_mzd_mul_uint64_128(void) {
  return test_mzd_mul_f("mul uint64 128", 128, 128, mzd_mul_v_uint64, false);
}

static int test_mzd_mul_uint64_192(void) {
  return test_mzd_mul_f("mul uint64 192", 192, 192, mzd_mul_v_uint64, false);
}

static int test_mzd_mul_uint64_256(void) {
  return test_mzd_mul_f("mul uint64 256", 256, 256, mzd_mul_v_uint64, false);
}

static int test_mzd_addmul_uint64_128(void) {
  return test_mzd_mul_f("addmul uint64 128", 128, 128, mzd_addmul_v_uint64, true);
}

static int test_mzd_addmul_uint64_192(void) {
  return test_mzd_mul_f("addmul uint64 192", 192, 192, mzd_addmul_v_uint64, true);
}

static int test_mzd_addmul_uint64_256(void) {
  return test_mzd_mul_f("addmul uint64 256", 256, 256, mzd_addmul_v_uint64, true);
}

#if defined(MUL_M4RI)
static int test_mzd_mull_uint64_128(void) {
  return test_mzd_mul_l_f("mull uint64 128", 128, 128, mzd_mul_vl_uint64, false);
}

static int test_mzd_mull_uint64_192(void) {
  return test_mzd_mul_l_f("mull uint64 192", 192, 192, mzd_mul_vl_uint64, false);
}

static int test_mzd_mull_uint64_256(void) {
  return test_mzd_mul_l_f("mull uint64 256", 256, 256, mzd_mul_vl_uint64, false);
}

static int test_mzd_addmull_uint64_128(void) {
  return test_mzd_mul_l_f("addmull uint64 128", 128, 128, mzd_addmul_vl_uint64, true);
}

static int test_mzd_addmull_uint64_192(void) {
  return test_mzd_mul_l_f("addmull uint64 192", 192, 192, mzd_addmul_vl_uint64, true);
}

static int test_mzd_addmull_uint64_256(void) {
  return test_mzd_mul_l_f("addmull uint64 256", 256, 256, mzd_addmul_vl_uint64, true);
}
#endif

#if defined(WITH_AVX2)
static int test_mzd_mul_s256(void) {
  return test_mzd_mul_f("mul s256", 192, 192, mzd_mul_v_s256, false);
}

static int test_mzd_mul_s256_128(void) {
  return test_mzd_mul_f("mul s256 128", 128, 128, mzd_mul_v_s256_128, false);
}

static int test_mzd_mul_s256_192(void) {
  return test_mzd_mul_f("mul s256 192", 192, 192, mzd_mul_v_s256_192, false);
}

static int test_mzd_mul_s256_256(void) {
  return test_mzd_mul_f("mul s256 256", 256, 256, mzd_mul_v_s256_256, false);
}

static int test_mzd_addmul_s256_128(void) {
  return test_mzd_mul_f("addmul s256 128", 128, 128, mzd_addmul_v_s256_128, true);
}

static int test_mzd_addmul_s256_192(void) {
  return test_mzd_mul_f("addmul s256 192", 192, 192, mzd_addmul_v_s256_192, true);
}

static int test_mzd_addmul_s256_256(void) {
  return test_mzd_mul_f("addmul s256 256", 256, 256, mzd_addmul_v_s256_256, true);
}

#if defined(MUL_M4RI)
static int test_mzd_mull_s256(void) {
  return test_mzd_mul_l_f("mull s256", 192, 192, mzd_mul_vl_s256, false);
}

static int test_mzd_mull_s256_128(void) {
  return test_mzd_mul_l_f("mull s256 128", 128, 128, mzd_mul_vl_s256_128, false);
}

static int test_mzd_mull_s256_192(void) {
  return test_mzd_mul_l_f("mull s256 192", 192, 192, mzd_mul_vl_s256_192, false);
}

static int test_mzd_mull_s256_256(void) {
  return test_mzd_mul_l_f("mull s256 256", 256, 256, mzd_mul_vl_s256_256, false);
}

static int test_mzd_addmull_s256(void) {
  return test_mzd_mul_l_f("addmull s256", 192, 192, mzd_addmul_vl_s256, true);
}

static int test_mzd_addmull_s256_128(void) {
  return test_mzd_mul_l_f("addmull s256 128", 128, 128, mzd_addmul_vl_s256_128, true);
}

static int test_mzd_addmull_s256_192(void) {
  return test_mzd_mul_l_f("addmull s256 192", 192, 192, mzd_addmul_vl_s256_192, true);
}

static int test_mzd_addmull_s256_256(void) {
  return test_mzd_mul_l_f("addmull s256 256", 256, 256, mzd_addmul_vl_s256_256, true);
}
#endif
#endif

#if defined(WITH_SSE2) || defined(WITH_NEON)
static int test_mzd_mul_s128(void) {
  return test_mzd_mul_f("mul s128", 192, 192, mzd_mul_v_s128, false);
}

static int test_mzd_mul_s128_128(void) {
  return test_mzd_mul_f("mul s128 128", 128, 128, mzd_mul_v_s128_128, false);
}

static int test_mzd_mul_s128_192(void) {
  return test_mzd_mul_f("mul s128 192", 192, 192, mzd_mul_v_s128_192, false);
}

static int test_mzd_mul_s128_256(void) {
  return test_mzd_mul_f("mul s128 256", 256, 256, mzd_mul_v_s128_256, false);
}

static int test_mzd_addmul_s128_128(void) {
  return test_mzd_mul_f("addmul s128 128", 128, 128, mzd_addmul_v_s128_128, true);
}

static int test_mzd_addmul_s128_192(void) {
  return test_mzd_mul_f("addmul s128 192", 192, 192, mzd_addmul_v_s128_192, true);
}

static int test_mzd_addmul_s128_256(void) {
  return test_mzd_mul_f("addmul s128 256", 256, 256, mzd_addmul_v_s128_256, true);
}

#if defined(MUL_M4RI)
static int test_mzd_mull_s128(void) {
  return test_mzd_mul_l_f("mull s128", 192, 192, mzd_mul_vl_s128, false);
}

static int test_mzd_mull_s128_128(void) {
  return test_mzd_mul_l_f("mull s128 128", 128, 128, mzd_mul_vl_s128_128, false);
}

static int test_mzd_mull_s128_192(void) {
  return test_mzd_mul_l_f("mull s128 192", 192, 192, mzd_mul_vl_s128_192, false);
}

static int test_mzd_mull_s128_256(void) {
  return test_mzd_mul_l_f("mull s128 256", 256, 256, mzd_mul_vl_s128_256, false);
}

static int test_mzd_addmull_s128(void) {
  return test_mzd_mul_l_f("addmull s128", 192, 192, mzd_addmul_vl_s128, true);
}

static int test_mzd_addmull_s128_128(void) {
  return test_mzd_mul_l_f("addmull s128 128", 128, 128, mzd_addmul_vl_s128_128, true);
}

static int test_mzd_addmull_s128_192(void) {
  return test_mzd_mul_l_f("addmull s128 192", 192, 192, mzd_addmul_vl_s128_192, true);
}

static int test_mzd_addmull_s128_256(void) {
  return test_mzd_mul_l_f("addmull s128 256", 256, 256, mzd_addmul_vl_s128_256, true);
}
#endif
#endif

int main(void) {
  int ret = 0;

  ret |= test_mzd_local_equal();
  ret |= test_mzd_mul_uint64_128();
  ret |= test_mzd_mul_uint64_192();
  ret |= test_mzd_mul_uint64_256();
  ret |= test_mzd_addmul_uint64_128();
  ret |= test_mzd_addmul_uint64_192();
  ret |= test_mzd_addmul_uint64_256();
#if defined(MUL_M4RI)
  ret |= test_mzd_mull_uint64_128();
  ret |= test_mzd_mull_uint64_192();
  ret |= test_mzd_mull_uint64_256();
  ret |= test_mzd_addmull_uint64_128();
  ret |= test_mzd_addmull_uint64_192();
  ret |= test_mzd_addmull_uint64_256();
#endif
#ifdef WITH_AVX2
  if (CPU_SUPPORTS_AVX2) {
    ret |= test_mzd_mul_s256();
    ret |= test_mzd_mul_s256_128();
    ret |= test_mzd_mul_s256_192();
    ret |= test_mzd_mul_s256_256();
    ret |= test_mzd_addmul_s256_128();
    ret |= test_mzd_addmul_s256_192();
    ret |= test_mzd_addmul_s256_256();
#if defined(MUL_M4RI)
    ret |= test_mzd_mull_s256();
    ret |= test_mzd_mull_s256_128();
    ret |= test_mzd_mull_s256_192();
    ret |= test_mzd_mull_s256_256();
    ret |= test_mzd_addmull_s256();
    ret |= test_mzd_addmull_s256_128();
    ret |= test_mzd_addmull_s256_192();
    ret |= test_mzd_addmull_s256_256();
#endif
  }
#endif
#if defined(WITH_SSE2) || defined(WITH_NEON)
  if (CPU_SUPPORTS_SSE2 || CPU_SUPPORTS_NEON) {
    ret |= test_mzd_mul_s128();
    ret |= test_mzd_mul_s128_128();
    ret |= test_mzd_mul_s128_192();
    ret |= test_mzd_mul_s128_256();
    ret |= test_mzd_addmul_s128_128();
    ret |= test_mzd_addmul_s128_192();
    ret |= test_mzd_addmul_s128_256();
#if defined(MUL_M4RI)
    ret |= test_mzd_mull_s128();
    ret |= test_mzd_mull_s128_128();
    ret |= test_mzd_mull_s128_192();
    ret |= test_mzd_mull_s128_256();
    ret |= test_mzd_addmull_s128();
    ret |= test_mzd_addmull_s128_128();
    ret |= test_mzd_addmull_s128_192();
    ret |= test_mzd_addmull_s128_256();
#endif
  }
#endif
  return ret;
}
