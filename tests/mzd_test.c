#include "../mzd_additional.h"

#ifdef WITH_OPT
#include "../simd.h"
#if !defined(_MSC_VER)
#include <stdalign.h>
#endif
#endif

#include <m4ri/m4ri.h>

#include "utils.h"

static void test_mzd_local_equal(void) {
  for (unsigned int i = 0; i < 10; ++i) {
    mzd_local_t* a = mzd_local_init(1, (i + 1) * 64);
    mzd_randomize_ssl(a);
    mzd_local_t* b = mzd_local_copy(NULL, a);

    if (mzd_local_equal(a, b)) {
      printf("equal: ok [%u]\n", (i + 1) * 64);
    }

    b = mzd_xor(b, b, a);
    if (mzd_local_equal(a, b)) {
      printf("equal: ok [%u]\n", (i + 1) * 64);
    }

    mzd_local_free(a);
    mzd_local_free(b);
  }
}

#ifdef WITH_AVX2
static int test_mzd_mul_avx(void) {
  int ret = 0;

  unsigned int size = 192;
  mzd_t* A    = mzd_init(size, size);
  mzd_t* v    = mzd_init(1, size);
  mzd_t* c    = mzd_init(1, size);

  mzd_randomize(A);
  mzd_randomize(v);
  mzd_randomize(c);

  mzd_local_t* Al = mzd_convert(A);
  mzd_local_t* vl = mzd_convert(v);
  mzd_local_t* c2 = mzd_convert(c);

  for (unsigned int k = 0; k < 3; ++k) {
    mzd_t* r  = mzd_mul_naive(c, v, A);
    mzd_local_t* rl = mzd_mul_v_avx(c2, vl, Al);

    mzd_local_t* rc = mzd_convert(r);

    if (!mzd_local_equal(rc, rl)) {
      printf("mul avx: fail [%u x %u]\n", size, size);
      ret = -1;
    } else {
      printf("mul avx: ok [%u x %u]\n", size, size);
    }

    mzd_local_free(rc);
  }

  mzd_free(A);
  mzd_free(v);
  mzd_free(c);

  mzd_local_free(c2);
  mzd_local_free(Al);
  mzd_local_free(vl);

  return ret;
}

static int test_mzd_mul_avx_128(void) {
  int ret = 0;

  unsigned int size = 128;
  mzd_t* A    = mzd_init(size, size);
  mzd_t* v    = mzd_init(1, size);
  mzd_t* c    = mzd_init(1, size);

  mzd_randomize(A);
  mzd_randomize(v);
  mzd_randomize(c);

  mzd_local_t* Al = mzd_convert(A);
  mzd_local_t* vl = mzd_convert(v);
  mzd_local_t* c2 = mzd_convert(c);

  for (unsigned int k = 0; k < 3; ++k) {
    mzd_t* r  = mzd_mul_naive(c, v, A);
    mzd_local_t* rl = mzd_mul_v_avx_128(c2, vl, Al);

    mzd_local_t* rc = mzd_convert(r);

    if (!mzd_local_equal(rc, rl)) {
      printf("mul avx 128: fail [%u x %u]\n", size, size);
      ret = -1;
    } else {
      printf("mul avx 128: ok [%u x %u]\n", size, size);
    }

    mzd_local_free(rc);
  }

  mzd_free(A);
  mzd_free(v);
  mzd_free(c);

  mzd_local_free(c2);
  mzd_local_free(Al);
  mzd_local_free(vl);

  return ret;
}
#endif

#ifdef WITH_NEON
static int test_mzd_mul_vl_neon_192(void) {
  int ret = 0;

  unsigned int size = 192;
  mzd_t* A    = mzd_init(size, size);
  mzd_t* v    = mzd_init(1, size);
  mzd_t* c    = mzd_init(1, size);

  mzd_randomize(A);
  mzd_randomize(v);
  mzd_randomize(c);

  mzd_local_t* Al = mzd_convert(A);
  mzd_local_t* vl = mzd_convert(v);
  mzd_local_t* c2 = mzd_convert(c);

  for (unsigned int k = 0; k < 3; ++k) {
    mzd_t* r  = mzd_mul_naive(c, v, A);
    mzd_local_t* rl = mzd_mul_v_neon(c2, vl, Al);

    mzd_local_t* rc = mzd_convert(r);

    if (!mzd_local_equal(rc, rl)) {
      printf("mul neon 192: fail [%u x %u]\n", size, size);
      ret = -1;
    } else {
      printf("mul neon 192: ok [%u x %u]\n", size, size);
    }

    mzd_local_free(rc);
  }

  mzd_free(A);
  mzd_free(v);
  mzd_free(c);

  mzd_local_free(c2);
  mzd_local_free(Al);
  mzd_local_free(vl);

  return ret;
}

static int test_mzd_mul_vl_neon_256(void) {
  int ret = 0;

  unsigned int size = 256;
  mzd_t* A    = mzd_init(size, size);
  mzd_t* v    = mzd_init(1, size);
  mzd_t* c    = mzd_init(1, size);

  mzd_randomize(A);
  mzd_randomize(v);
  mzd_randomize(c);

  mzd_local_t* Al = mzd_convert(A);
  mzd_local_t* vl = mzd_convert(v);
  mzd_local_t* c2 = mzd_convert(c);
  mzd_local_t* All = mzd_precompute_matrix_lookup(Al);

  for (unsigned int k = 0; k < 3; ++k) {
    mzd_t* r  = mzd_mul_naive(c, v, A);
    mzd_local_t* rl = mzd_mul_v_neon(c2, vl, Al);

    mzd_local_t* rc = mzd_convert(r);

    if (!mzd_local_equal(rc, rl)) {
      printf("mul neon 256: fail [%u x %u]\n", size, size);
      ret = -1;
    } else {
      printf("mul neon 256: ok [%u x %u]\n", size, size);
    }

    mzd_local_free(rc);
  }

  mzd_free(A);
  mzd_free(v);
  mzd_free(c);

  mzd_local_free(All);
  mzd_local_free(c2);
  mzd_local_free(Al);
  mzd_local_free(vl);

  return ret;
}

static int test_mzd_addmul_vl_neon_192(void) {
  int ret = 0;

  unsigned int size = 192;
  mzd_t* A    = mzd_init(size, size);
  mzd_t* v    = mzd_init(1, size);
  mzd_t* c    = mzd_init(1, size);

  mzd_randomize(A);
  mzd_randomize(v);
  mzd_randomize(c);

  mzd_local_t* Al = mzd_convert(A);
  mzd_local_t* vl = mzd_convert(v);
  mzd_local_t* c2 = mzd_convert(c);

  for (unsigned int k = 0; k < 3; ++k) {
    mzd_t* r  = mzd_addmul_naive(c, v, A);
    mzd_local_t* rl = mzd_addmul_v_neon(c2, vl, Al);

    mzd_local_t* rc = mzd_convert(r);

    if (!mzd_local_equal(rc, rl)) {
      printf("addmul neon 192: fail [%u x %u]\n", size, size);
      ret = -1;
    } else {
      printf("addmul neon 192: ok [%u x %u]\n", size, size);
    }

    mzd_local_free(rc);
  }

  mzd_free(A);
  mzd_free(v);
  mzd_free(c);

  mzd_local_free(c2);
  mzd_local_free(Al);
  mzd_local_free(vl);

  return ret;
}

static int test_mzd_addmul_vl_neon_256(void) {
  int ret = 0;

  unsigned int size = 256;
  mzd_t* A    = mzd_init(size, size);
  mzd_t* v    = mzd_init(1, size);
  mzd_t* c    = mzd_init(1, size);

  mzd_randomize(A);
  mzd_randomize(v);
  mzd_randomize(c);

  mzd_local_t* Al = mzd_convert(A);
  mzd_local_t* vl = mzd_convert(v);
  mzd_local_t* c2 = mzd_convert(c);

  for (unsigned int k = 0; k < 3; ++k) {
    mzd_t* r  = mzd_addmul_naive(c, v, A);
    mzd_local_t* rl = mzd_addmul_v_neon(c2, vl, Al);

    mzd_local_t* rc = mzd_convert(r);

    if (!mzd_local_equal(rc, rl)) {
      printf("addmul neon 256: fail [%u x %u]\n", size, size);
      ret = -1;
    } else {
      printf("addmul neon 256: ok [%u x %u]\n", size, size);
    }

    mzd_local_free(rc);
  }

  mzd_free(A);
  mzd_free(v);
  mzd_free(c);

  mzd_local_free(c2);
  mzd_local_free(Al);
  mzd_local_free(vl);

  return ret;
}
#endif

static void test_mzd_mul(void) {
  for (unsigned int i = 1; i <= 10; ++i) {
    for (unsigned int j = 1; j <= 10; ++j) {
      mzd_t* A = mzd_init(i * 64, j * 64);
      mzd_t* v = mzd_init(1, i * 64);
      mzd_t* c = mzd_init(1, j * 64);

      mzd_randomize(A);
      mzd_randomize(v);
      mzd_randomize(c);

      mzd_local_t* Al  = mzd_convert(A);
      mzd_local_t* All = mzd_precompute_matrix_lookup(Al);
      mzd_local_t* vl  = mzd_convert(v);
      mzd_local_t* cl  = mzd_convert(c);
      mzd_local_t* cll = mzd_convert(c);

      mzd_t* At = mzd_transpose(NULL, A);
      mzd_t* vt = mzd_transpose(NULL, v);
      mzd_t* c2 = mzd_copy(NULL, c);
      mzd_t* c3 = mzd_transpose(NULL, c);

      for (unsigned int k = 0; k < 3; ++k) {
        mzd_local_t* r  = mzd_mul_v(cl, vl, Al);
        mzd_local_t* rl = mzd_mul_vl(cll, vl, All);
        mzd_t* r2 = mzd_mul(c2, v, A, __M4RI_STRASSEN_MUL_CUTOFF);
        mzd_t* r3 = mzd_mul(c3, At, vt, __M4RI_STRASSEN_MUL_CUTOFF);

        if (!mzd_local_equal(r, rl)) {
          printf("mul: fail [%u x %u]\n", i * 64, j * 64);
        }

        mzd_local_t* rc = mzd_convert(r2);
        if (!mzd_local_equal(r, rc)) {
          printf("mul: fail [%u x %u]\n", i * 64, j * 64);
        }
        mzd_local_free(rc);

        mzd_t* r4 = mzd_transpose(NULL, r3);
        if (mzd_cmp(r4, r2) != 0) {
          printf("mul: fail [%u x %u]\n", i * 64, j * 64);
        }
        mzd_free(r4);
      }

      mzd_free(At);
      mzd_free(A);
      mzd_free(v);
      mzd_free(c);

      mzd_free(vt);
      mzd_free(c2);
      mzd_free(c3);

      mzd_local_free(All);
      mzd_local_free(Al);
      mzd_local_free(cll);
      mzd_local_free(cl);
      mzd_local_free(vl);
    }
  }
}

static void test_mzd_shift(void) {
#ifdef WITH_OPT
#ifdef WITH_SSE2
  if (CPU_SUPPORTS_SSE2) {
    mzd_local_t* v = mzd_local_init(1, 128);
    mzd_local_t* w = mzd_local_copy(NULL, v);
    mzd_local_t* r = mzd_local_copy(NULL, v);
    __m128i* wr    = __builtin_assume_aligned(FIRST_ROW(w), 16);

    for (unsigned int i = 0; i < 32; ++i) {
      mzd_randomize_ssl(v);
      mzd_local_copy(w, v);

      mzd_shift_left(r, v, i);
      *wr = mm128_shift_left(*wr, i);

      if (!mzd_local_equal(r, w)) {
        printf("lshift fail\n");
      }
    }

    for (unsigned int i = 0; i < 32; ++i) {
      mzd_randomize_ssl(v);
      mzd_local_copy(w, v);

      mzd_shift_right(r, v, i);
      *wr = mm128_shift_right(*wr, i);

      if (!mzd_local_equal(r, w)) {
        printf("rshift fail\n");
      }
    }

    mzd_local_free(w);
    mzd_local_free(v);
    mzd_local_free(r);
  }
#endif
#ifdef WITH_AVX2
  if (CPU_SUPPORTS_AVX2) {
    mzd_local_t* v = mzd_local_init(1, 256);
    mzd_local_t* w = mzd_local_copy(NULL, v);
    mzd_local_t* r = mzd_local_copy(NULL, v);
    __m256i* wr    = __builtin_assume_aligned(FIRST_ROW(w), 32);

    for (unsigned int i = 0; i < 32; ++i) {
      mzd_randomize_ssl(v);
      mzd_local_copy(w, v);

      mzd_shift_left(r, v, i);
      *wr = mm256_shift_left(*wr, i);

      if (!mzd_local_equal(r, w)) {
        printf("lshift fail\n");
      }
    }

    for (unsigned int i = 0; i < 32; ++i) {
      mzd_randomize_ssl(v);
      mzd_local_copy(w, v);

      mzd_shift_right(r, v, i);
      *wr = mm256_shift_right(*wr, i);

      if (!mzd_local_equal(r, w)) {
        printf("rshift fail\n");
      }
    }

    mzd_local_free(w);
    mzd_local_free(v);
    mzd_local_free(r);
  }
#endif
#ifdef WITH_NEON
  if (CPU_SUPPORTS_NEON) {
    mzd_local_t* v = mzd_local_init(1, 256);
    mzd_local_t* w = mzd_local_copy(NULL, v);
    mzd_local_t* r = mzd_local_copy(NULL, v);
    uint32x4_t* wr = __builtin_assume_aligned(FIRST_ROW(w), alignof(uint32x4_t));

    for (unsigned int i = 0; i < 32; ++i) {
      mzd_randomize_ssl(v);
      mzd_local_copy(w, v);

      mzd_shift_left(r, v, i);
      mm256_shift_left(wr, wr, i);

      if (!mzd_local_equal(r, w)) {
        printf("lshift fail\n");
      }
    }

    for (unsigned int i = 0; i < 32; ++i) {
      mzd_randomize_ssl(v);
      mzd_local_copy(w, v);

      mzd_shift_right(r, v, i);
      mm256_shift_right(wr, wr, i);

      if (!mzd_local_equal(r, w)) {
        printf("rshift fail\n");
      }
    }

    mzd_local_free(w);
    mzd_local_free(v);
    mzd_local_free(r);
  }
#endif
#endif
}

int main() {
  test_mzd_local_equal();
  test_mzd_mul();
#ifdef WITH_AVX2
  if (CPU_SUPPORTS_AVX2) {
    test_mzd_mul_avx();
    test_mzd_mul_avx_128();
  }
#endif
  test_mzd_shift();
#ifdef WITH_NEON
  test_mzd_mul_vl_neon_192();
  test_mzd_mul_vl_neon_256();

  test_mzd_addmul_vl_neon_192();
  test_mzd_addmul_vl_neon_256();
#endif
}
