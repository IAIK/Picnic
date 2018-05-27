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
    mzd_local_t* b = mzd_local_copy(NULL, a);

    if (mzd_local_equal(a, b)) {
      printf("equal: ok [%u]\n", (i + 1) * 64);
    } else {
      printf("equal: fail [%u]\n", (i + 1) * 64);
      ret = -1;
    }

    mzd_xor(b, b, a);
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

static int test_mzd_mul_f(const char* n, unsigned int rows, unsigned int cols, mul_fn f, bool is_addmul) {
  int ret = 0;

  mzd_t* A    = mzd_init(rows, cols);
  mzd_t* v    = mzd_init(1, rows);
  mzd_t* c    = mzd_init(1, cols);

  mzd_randomize(A);
  mzd_randomize(v);
  mzd_randomize(c);

  mzd_local_t* Al = mzd_convert(A);
  mzd_local_t* vl = mzd_convert(v);
  mzd_local_t* c2 = mzd_convert(c);

  for (unsigned int k = 0; k < 3; ++k) {
    mzd_t* r  = is_addmul ? mzd_addmul_naive(c, v, A) : mzd_mul_naive(c, v, A);
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

static int test_mzd_mul_l_f(const char* n, unsigned int rows, unsigned int cols, mul_fn f, bool is_addmul) {
  int ret = 0;

  mzd_t* A    = mzd_init(rows, cols);
  mzd_t* v    = mzd_init(1, rows);
  mzd_t* c    = mzd_init(1, cols);

  mzd_randomize(A);
  mzd_randomize(v);
  mzd_randomize(c);

  mzd_local_t* Al = mzd_convert(A);
  mzd_local_t* All = mzd_precompute_matrix_lookup(Al);
  mzd_local_t* vl = mzd_convert(v);
  mzd_local_t* c2 = mzd_convert(c);

  for (unsigned int k = 0; k < 3; ++k) {
    mzd_t* r  = is_addmul ? mzd_addmul_naive(c, v, A) : mzd_mul_naive(c, v, A);
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

#ifdef WITH_AVX2
static int test_mzd_mul_avx(void) {
  return test_mzd_mul_f("mul avx", 192, 192, mzd_mul_v_avx, false);
}

static int test_mzd_mul_avx_128(void) {
  return test_mzd_mul_f("mul avx 128", 128, 128, mzd_mul_v_avx_128, false);
}

static int test_mzd_mul_avx_192(void) {
  return test_mzd_mul_f("mul avx 192", 192, 192, mzd_mul_v_avx_192, false);
}

static int test_mzd_mul_avx_256(void) {
  return test_mzd_mul_f("mul avx 256", 256, 256, mzd_mul_v_avx_256, false);
}

static int test_mzd_addmul_avx(void) {
  return test_mzd_mul_f("addmul avx", 192, 192, mzd_addmul_v_avx, true);
}

static int test_mzd_addmul_avx_128(void) {
  return test_mzd_mul_f("addmul avx 128", 128, 128, mzd_addmul_v_avx_128, true);
}

static int test_mzd_addmul_avx_192(void) {
  return test_mzd_mul_f("addmul avx 192", 192, 192, mzd_addmul_v_avx_192, true);
}

static int test_mzd_addmul_avx_256(void) {
  return test_mzd_mul_f("addmul avx 256", 256, 256, mzd_addmul_v_avx_256, true);
}

static int test_mzd_mull_avx(void) {
  return test_mzd_mul_l_f("mull avx", 192, 192, mzd_mul_vl_avx, false);
}

static int test_mzd_mull_avx_128(void) {
  return test_mzd_mul_l_f("mull avx 128", 128, 128, mzd_mul_vl_avx_128, false);
}

static int test_mzd_mull_avx_192(void) {
  return test_mzd_mul_l_f("mull avx 192", 192, 192, mzd_mul_vl_avx_192, false);
}

static int test_mzd_mull_avx_256(void) {
  return test_mzd_mul_l_f("mull avx 256", 256, 256, mzd_mul_vl_avx_256, false);
}

static int test_mzd_addmull_avx(void) {
  return test_mzd_mul_l_f("addmull avx", 192, 192, mzd_addmul_vl_avx, true);
}

static int test_mzd_addmull_avx_128(void) {
  return test_mzd_mul_l_f("addmull avx 128", 128, 128, mzd_addmul_vl_avx_128, true);
}

static int test_mzd_addmull_avx_192(void) {
  return test_mzd_mul_l_f("addmull avx 192", 192, 192, mzd_addmul_vl_avx_192, true);
}

static int test_mzd_addmull_avx_256(void) {
  return test_mzd_mul_l_f("addmull avx 256", 256, 256, mzd_addmul_vl_avx_256, true);
}
#endif

#ifdef WITH_SSE2
static int test_mzd_mul_sse(void) {
  return test_mzd_mul_f("mul sse", 192, 192, mzd_mul_v_sse, false);
}

static int test_mzd_mul_sse_128(void) {
  return test_mzd_mul_f("mul sse 128", 128, 128, mzd_mul_v_sse_128, false);
}

static int test_mzd_mul_sse_192(void) {
  return test_mzd_mul_f("mul sse 192", 192, 192, mzd_mul_v_sse_192, false);
}

static int test_mzd_mul_sse_256(void) {
  return test_mzd_mul_f("mul sse 256", 256, 256, mzd_mul_v_sse_256, false);
}

static int test_mzd_addmul_sse(void) {
  return test_mzd_mul_f("addmul sse", 192, 192, mzd_addmul_v_sse, true);
}

static int test_mzd_addmul_sse_128(void) {
  return test_mzd_mul_f("addmul sse 128", 128, 128, mzd_addmul_v_sse_128, true);
}

static int test_mzd_addmul_sse_192(void) {
  return test_mzd_mul_f("addmul sse 192", 192, 192, mzd_addmul_v_sse_192, true);
}

static int test_mzd_addmul_sse_256(void) {
  return test_mzd_mul_f("addmul sse 256", 256, 256, mzd_addmul_v_sse_256, true);
}

static int test_mzd_mull_sse(void) {
  return test_mzd_mul_l_f("mull sse", 192, 192, mzd_mul_vl_sse, false);
}

static int test_mzd_mull_sse_128(void) {
  return test_mzd_mul_l_f("mull sse 128", 128, 128, mzd_mul_vl_sse_128, false);
}

static int test_mzd_mull_sse_192(void) {
  return test_mzd_mul_l_f("mull sse 192", 192, 192, mzd_mul_vl_sse_192, false);
}

static int test_mzd_mull_sse_256(void) {
  return test_mzd_mul_l_f("mull sse 256", 256, 256, mzd_mul_vl_sse_256, false);
}

static int test_mzd_addmull_sse(void) {
  return test_mzd_mul_l_f("addmull sse", 192, 192, mzd_addmul_vl_sse, true);
}

static int test_mzd_addmull_sse_128(void) {
  return test_mzd_mul_l_f("addmull sse 128", 128, 128, mzd_addmul_vl_sse_128, true);
}

static int test_mzd_addmull_sse_192(void) {
  return test_mzd_mul_l_f("addmull sse 192", 192, 192, mzd_addmul_vl_sse_192, true);
}

static int test_mzd_addmull_sse_256(void) {
  return test_mzd_mul_l_f("addmull sse 256", 256, 256, mzd_addmul_vl_sse_256, true);
}
#endif

#ifdef WITH_NEON
static int test_mzd_mul_neon(void) {
  return test_mzd_mul_f("mul neon", 192, 192, mzd_mul_v_neon, false);
}

static int test_mzd_mul_neon_128(void) {
  return test_mzd_mul_f("mul neon 128", 128, 128, mzd_mul_v_neon_128, false);
}

static int test_mzd_mul_neon_192(void) {
  return test_mzd_mul_f("mul neon 192", 192, 192, mzd_mul_v_neon_192, false);
}

static int test_mzd_mul_neon_256(void) {
  return test_mzd_mul_f("mul neon 256", 256, 256, mzd_mul_v_neon_256, false);
}

static int test_mzd_addmul_neon(void) {
  return test_mzd_mul_f("addmul neon", 192, 192, mzd_addmul_v_neon, true);
}

static int test_mzd_addmul_neon_128(void) {
  return test_mzd_mul_f("addmul neon 128", 128, 128, mzd_addmul_v_neon_128, true);
}

static int test_mzd_addmul_neon_192(void) {
  return test_mzd_mul_f("addmul neon 192", 192, 192, mzd_addmul_v_neon_192, true);
}

static int test_mzd_addmul_neon_256(void) {
  return test_mzd_mul_f("addmul neon 256", 256, 256, mzd_addmul_v_neon_256, true);
}

static int test_mzd_mull_neon(void) {
  return test_mzd_mul_l_f("mull neon", 192, 192, mzd_mul_vl_neon, false);
}

static int test_mzd_mull_neon_128(void) {
  return test_mzd_mul_l_f("mull neon 128", 128, 128, mzd_mul_vl_neon /* _128 */, false);
}

static int test_mzd_mull_neon_192(void) {
  return test_mzd_mul_l_f("mull neon 192", 192, 192, mzd_mul_vl_neon /* _192 */, false);
}

static int test_mzd_mull_neon_256(void) {
  return test_mzd_mul_l_f("mull neon 256", 256, 256, mzd_mul_vl_neon /* _256 */, false);
}

static int test_mzd_addmull_neon(void) {
  return test_mzd_mul_l_f("addmull neon", 192, 192, mzd_addmul_vl_neon, true);
}

static int test_mzd_addmull_neon_128(void) {
  return test_mzd_mul_l_f("addmull neon 128", 128, 128, mzd_addmul_vl_neon /* _128 */, true);
}

static int test_mzd_addmull_neon_192(void) {
  return test_mzd_mul_l_f("addmull neon 192", 192, 192, mzd_addmul_vl_neon /* _192 */, true);
}

static int test_mzd_addmull_neon_256(void) {
  return test_mzd_mul_l_f("addmull neon 256", 256, 256, mzd_addmul_vl_neon /*_256 */, true);
}
#endif

static int test_mzd_mul(void) {
  int ret = 0;

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
        mzd_mul_v(cl, vl, Al);
        mzd_mul_vl(cll, vl, All);
        mzd_t* r2 = mzd_mul(c2, v, A, __M4RI_STRASSEN_MUL_CUTOFF);
        mzd_t* r3 = mzd_mul(c3, At, vt, __M4RI_STRASSEN_MUL_CUTOFF);

        if (!mzd_local_equal(cl, cll)) {
          printf("mul: fail [%u x %u]\n", i * 64, j * 64);
          ret = -1;
        }

        mzd_local_t* rc = mzd_convert(r2);
        if (!mzd_local_equal(cl, rc)) {
          printf("mul: fail [%u x %u]\n", i * 64, j * 64);
          ret = -1;
        }
        mzd_local_free(rc);

        mzd_t* r4 = mzd_transpose(NULL, r3);
        if (mzd_cmp(r4, r2) != 0) {
          printf("mul: fail [%u x %u]\n", i * 64, j * 64);
          ret = -1;
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

  return ret;
}

static int test_mzd_shift(void) {
  int ret = 0;
#ifdef WITH_CUSTOM_INSTANCES
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
        ret = -1;
      }
    }

    for (unsigned int i = 0; i < 32; ++i) {
      mzd_randomize_ssl(v);
      mzd_local_copy(w, v);

      mzd_shift_right(r, v, i);
      *wr = mm128_shift_right(*wr, i);

      if (!mzd_local_equal(r, w)) {
        printf("rshift fail\n");
        ret = -1;
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
        ret = -1;
      }
    }

    for (unsigned int i = 0; i < 32; ++i) {
      mzd_randomize_ssl(v);
      mzd_local_copy(w, v);

      mzd_shift_right(r, v, i);
      *wr = mm256_shift_right(*wr, i);

      if (!mzd_local_equal(r, w)) {
        printf("rshift fail\n");
        ret = -1;
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
        ret = -1;
      }
    }

    for (unsigned int i = 0; i < 32; ++i) {
      mzd_randomize_ssl(v);
      mzd_local_copy(w, v);

      mzd_shift_right(r, v, i);
      mm256_shift_right(wr, wr, i);

      if (!mzd_local_equal(r, w)) {
        printf("rshift fail\n");
        ret = -1;
      }
    }

    mzd_local_free(w);
    mzd_local_free(v);
    mzd_local_free(r);
  }
#endif
#endif
#endif
  return ret;
}

int main() {
  int ret = 0;

  ret |= test_mzd_local_equal();
  ret |= test_mzd_mul();
  ret |= test_mzd_mul_uint64_128();
  ret |= test_mzd_mul_uint64_192();
  ret |= test_mzd_mul_uint64_256();
  ret |= test_mzd_addmul_uint64_128();
  ret |= test_mzd_addmul_uint64_192();
  ret |= test_mzd_addmul_uint64_256();
  ret |= test_mzd_mull_uint64_128();
  ret |= test_mzd_mull_uint64_192();
  ret |= test_mzd_mull_uint64_256();
  ret |= test_mzd_addmull_uint64_128();
  ret |= test_mzd_addmull_uint64_192();
  ret |= test_mzd_addmull_uint64_256();
#ifdef WITH_AVX2
  if (CPU_SUPPORTS_AVX2) {
    ret |= test_mzd_mul_avx();
    ret |= test_mzd_mul_avx_128();
    ret |= test_mzd_mul_avx_192();
    ret |= test_mzd_mul_avx_256();
    ret |= test_mzd_addmul_avx();
    ret |= test_mzd_addmul_avx_128();
    ret |= test_mzd_addmul_avx_192();
    ret |= test_mzd_addmul_avx_256();
    ret |= test_mzd_mull_avx();
    ret |= test_mzd_mull_avx_128();
    ret |= test_mzd_mull_avx_192();
    ret |= test_mzd_mull_avx_256();
    ret |= test_mzd_addmull_avx();
    ret |= test_mzd_addmull_avx_128();
    ret |= test_mzd_addmull_avx_192();
    ret |= test_mzd_addmull_avx_256();
  }
#endif
#ifdef WITH_SSE2
  if (CPU_SUPPORTS_SSE2) {
    ret |= test_mzd_mul_sse();
    ret |= test_mzd_mul_sse_128();
    ret |= test_mzd_mul_sse_192();
    ret |= test_mzd_mul_sse_256();
    ret |= test_mzd_addmul_sse();
    ret |= test_mzd_addmul_sse_128();
    ret |= test_mzd_addmul_sse_192();
    ret |= test_mzd_addmul_sse_256();
    ret |= test_mzd_mull_sse();
    ret |= test_mzd_mull_sse_128();
    ret |= test_mzd_mull_sse_192();
    ret |= test_mzd_mull_sse_256();
    ret |= test_mzd_addmull_sse();
    ret |= test_mzd_addmull_sse_128();
    ret |= test_mzd_addmull_sse_192();
    ret |= test_mzd_addmull_sse_256();
  }
#endif
#ifdef WITH_NEON
  if (CPU_SUPPORTS_NEON) {
    ret |= test_mzd_mul_neon();
    ret |= test_mzd_mul_neon_128();
    ret |= test_mzd_mul_neon_192();
    ret |= test_mzd_mul_neon_256();
    ret |= test_mzd_addmul_neon();
    ret |= test_mzd_addmul_neon_128();
    ret |= test_mzd_addmul_neon_192();
    ret |= test_mzd_addmul_neon_256();
    ret |= test_mzd_mull_neon();
    ret |= test_mzd_mull_neon_128();
    ret |= test_mzd_mull_neon_192();
    ret |= test_mzd_mull_neon_256();
    ret |= test_mzd_addmull_neon();
    ret |= test_mzd_addmull_neon_128();
    ret |= test_mzd_addmull_neon_192();
    ret |= test_mzd_addmull_neon_256();
  }
#endif
  ret |= test_mzd_shift();
  return ret;
}
