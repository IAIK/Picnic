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

  for (unsigned int i = 1; i < 4; ++i) {
    const unsigned int cols = (i + 1) * 64;

    mzd_t* ma = mzd_init(1, cols);
    mzd_randomize(ma);

    mzd_local_t* a = mzd_convert(ma);
    mzd_free(ma);

    mzd_local_t* b = mzd_local_init(1, cols);
    memcpy(BLOCK(b, 0)->w64, CONST_BLOCK(a, 0)->w64, cols / 8);

    if (mzd_local_equal(a, b, 1, cols)) {
      printf("equal: ok [%u]\n", cols);
    } else {
      printf("equal: fail [%u]\n", cols);
      ret = -1;
    }

    switch (i) {
      case 1:
        mzd_xor_uint64_128(b, b, a);
        break;
      case 2:
        mzd_xor_uint64_192(b, b, a);
        break;
      case 3:
        mzd_xor_uint64_256(b, b, a);
        break;
    }

    if (!mzd_local_equal(a, b, 1, cols)) {
      printf("equal: ok [%u]\n", cols);
    } else {
      printf("equal: fail [%u]\n", cols);
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

    if (!mzd_local_equal(rc, c2, 1, cols)) {
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
  mzd_local_t* All = mzd_precompute_matrix_lookup(Al, rows, cols);
  mzd_local_t* vl  = mzd_convert(v);
  mzd_local_t* c2  = mzd_convert(c);

  for (unsigned int k = 0; k < 3; ++k) {
    mzd_t* r = is_addmul ? mzd_addmul_naive(c, v, A) : mzd_mul_naive(c, v, A);
    f(c2, vl, All);

    mzd_local_t* rc = mzd_convert(r);

    if (!mzd_local_equal(rc, c2, 1, cols)) {
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
  return test_mzd_mul_f("mul uint64 128", 128, 128, mzd_mul_v_uint64_128, false);
}

static int test_mzd_mul_uint64_192(void) {
  return test_mzd_mul_f("mul uint64 192", 192, 192, mzd_mul_v_uint64_192, false);
}

static int test_mzd_mul_uint64_256(void) {
  return test_mzd_mul_f("mul uint64 256", 256, 256, mzd_mul_v_uint64_256, false);
}

static int test_mzd_mul_uint64_128_576(void) {
  return test_mzd_mul_f("mul uint64 128 576", 128, 576, mzd_mul_v_uint64_128_576, false);
}

static int test_mzd_mul_uint64_128_640(void) {
  return test_mzd_mul_f("mul uint64 128 640", 128, 640, mzd_mul_v_uint64_128_640, false);
}

static int test_mzd_mul_uint64_192_896(void) {
  return test_mzd_mul_f("mul uint64 192 896", 192, 896, mzd_mul_v_uint64_192_896, false);
}

static int test_mzd_mul_uint64_192_960(void) {
  return test_mzd_mul_f("mul uint64 192 960", 192, 960, mzd_mul_v_uint64_192_960, false);
}

static int test_mzd_mul_uint64_256_1152(void) {
  return test_mzd_mul_f("mul uint64 256 1152", 256, 1152, mzd_mul_v_uint64_256_1152, false);
}

static int test_mzd_mul_uint64_256_1216(void) {
  return test_mzd_mul_f("mul uint64 256 1216", 256, 1216, mzd_mul_v_uint64_256_1216, false);
}

static int test_mzd_addmul_uint64_128(void) {
  return test_mzd_mul_f("addmul uint64 128", 128, 128, mzd_addmul_v_uint64_128, true);
}

static int test_mzd_addmul_uint64_192(void) {
  return test_mzd_mul_f("addmul uint64 192", 192, 192, mzd_addmul_v_uint64_192, true);
}

static int test_mzd_addmul_uint64_256(void) {
  return test_mzd_mul_f("addmul uint64 256", 256, 256, mzd_addmul_v_uint64_256, true);
}

#if defined(WITH_AVX2)
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
#endif

#if defined(WITH_SSE2) || defined(WITH_NEON)
static int test_mzd_mul_s128_128(void) {
  return test_mzd_mul_f("mul s128 128", 128, 128, mzd_mul_v_s128_128, false);
}

static int test_mzd_mul_s128_192(void) {
  return test_mzd_mul_f("mul s128 192", 192, 192, mzd_mul_v_s128_192, false);
}

static int test_mzd_mul_s128_256(void) {
  return test_mzd_mul_f("mul s128 256", 256, 256, mzd_mul_v_s128_256, false);
}

static int test_mzd_mul_s128_128_640(void) {
  return test_mzd_mul_f("mul s128 128 640", 128, 640, mzd_mul_v_s128_128_640, false);
}

static int test_mzd_mul_s128_192_896(void) {
  return test_mzd_mul_f("mul s128 192 896", 192, 896, mzd_mul_v_s128_192_896, false);
}

static int test_mzd_mul_s128_192_1024(void) {
  return test_mzd_mul_f("mul s128 192 1024", 192, 1024, mzd_mul_v_s128_192_1024, false);
}

static int test_mzd_mul_s128_256_1152(void) {
  return test_mzd_mul_f("mul s128 256 1152", 256, 1152, mzd_mul_v_s128_256_1152, false);
}

static int test_mzd_mul_s128_256_1280(void) {
  return test_mzd_mul_f("mul s128 256 1280", 256, 1280, mzd_mul_v_s128_256_1280, false);
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
#endif

int main(void) {
  int ret = 0;

  ret |= test_mzd_local_equal();
  ret |= test_mzd_mul_uint64_128();
  ret |= test_mzd_mul_uint64_192();
  ret |= test_mzd_mul_uint64_256();
  ret |= test_mzd_mul_uint64_128_576();
  ret |= test_mzd_mul_uint64_128_640();
  ret |= test_mzd_mul_uint64_192_896();
  ret |= test_mzd_mul_uint64_192_960();
  ret |= test_mzd_mul_uint64_256_1152();
  ret |= test_mzd_mul_uint64_256_1216();
  ret |= test_mzd_addmul_uint64_128();
  ret |= test_mzd_addmul_uint64_192();
  ret |= test_mzd_addmul_uint64_256();
#ifdef WITH_AVX2
  if (CPU_SUPPORTS_AVX2) {
    ret |= test_mzd_mul_s256_128();
    ret |= test_mzd_mul_s256_192();
    ret |= test_mzd_mul_s256_256();
    ret |= test_mzd_addmul_s256_128();
    ret |= test_mzd_addmul_s256_192();
    ret |= test_mzd_addmul_s256_256();
  }
#endif
#if defined(WITH_SSE2) || defined(WITH_NEON)
  if (CPU_SUPPORTS_SSE2 || CPU_SUPPORTS_NEON) {
    ret |= test_mzd_mul_s128_128();
    ret |= test_mzd_mul_s128_128_640();
    ret |= test_mzd_mul_s128_192();
    ret |= test_mzd_mul_s128_192_896();
    ret |= test_mzd_mul_s128_192_1024();
    ret |= test_mzd_mul_s128_256();
    ret |= test_mzd_mul_s128_256_1152();
    ret |= test_mzd_mul_s128_256_1280();
    ret |= test_mzd_addmul_s128_128();
    ret |= test_mzd_addmul_s128_192();
    ret |= test_mzd_addmul_s128_256();
  }
#endif
  return ret;
}
