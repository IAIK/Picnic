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
#endif
#include "utils.h"

#include <stdio.h>

static int test_uint64_128(void) {
  int ret = 0;

  mzd_local_t val = {{ 0x01, 0x02 }};
  mzd_local_t tmp;

  for (unsigned int i = 1; i <= 62; ++i) {
    mzd_shift_left_uint64_128(&tmp, &val, i);
    mzd_shift_right_uint64_128(&tmp, &tmp, i);

    if (!mzd_local_equal(&val, &tmp, 1, 128)) {
      printf("shift left right fail: %u\n", i);
      ret = -1;
    }
  }

  val.w64[0] = val.w64[1] = UINT64_C(0xfedcba9876543210);
  {
    const mzd_local_t cval = {{UINT64_C(0x7654321000000000), UINT64_C(0x76543210fedcba98)}};
    mzd_shift_left_uint64_128(&tmp, &val, 32);
    if (!mzd_local_equal(&cval, &tmp, 1, 128)) {
      printf("shift left 32\n");
      ret = -1;
    }
  }

  {
    const mzd_local_t cval = {{UINT64_C(0x76543210fedcba98), UINT64_C(0xfedcba98)}};
    mzd_shift_right_uint64_128(&tmp, &val, 32);
    if (!mzd_local_equal(&cval, &tmp, 1, 128)) {
      printf("shift right 32\n");
      ret = -1;
    }
  }

  return ret;
}

static int test_uint64_256(void) {
  int ret = 0;

  mzd_local_t val = {{ 0x01, 0x02, 0x03, 0x04 }};
  mzd_local_t tmp;

  for (unsigned int i = 1; i <= 61; ++i) {
    mzd_shift_left_uint64_256(&tmp, &val, i);
    mzd_shift_right_uint64_256(&tmp, &tmp, i);

    if (!mzd_local_equal(&val, &tmp, 1, 256)) {
      printf("shift left right fail: %u\n", i);
      ret = -1;
    }

    mzd_rotate_left_uint64_256(&tmp, &val, i);
    mzd_rotate_right_uint64_256(&tmp, &tmp, i);

    if (!mzd_local_equal(&val, &tmp, 1, 256)) {
      printf("rotate left right fail: %u\n", i);
      ret = -1;
    }
  }

  val.w64[0] = val.w64[1] = val.w64[2] = val.w64[3] = UINT64_C(0xfedcba9876543210);
  {
    const mzd_local_t cval = {{UINT64_C(0x7654321000000000), UINT64_C(0x76543210fedcba98),
                               UINT64_C(0x76543210fedcba98), UINT64_C(0x76543210fedcba98)}};
    mzd_shift_left_uint64_256(&tmp, &val, 32);
    if (!mzd_local_equal(&cval, &tmp, 1, 256)) {
      printf("shift left 32\n");
      ret = -1;
    }
  }

  {
    const mzd_local_t cval = {{UINT64_C(0x76543210fedcba98), UINT64_C(0x76543210fedcba98),
                               UINT64_C(0x76543210fedcba98), UINT64_C(0xfedcba98)}};
    mzd_shift_right_uint64_256(&tmp, &val, 32);
    if (!mzd_local_equal(&cval, &tmp, 1, 256)) {
      printf("shift right 32\n");
      ret = -1;
    }
  }

  return ret;
}

#if defined(WITH_OPT)
#if defined(WITH_SSE2) || defined(WITH_NEON)
ATTR_TARGET_S128
static int test_s128_128(void) {
  int ret = 0;

  mzd_local_t val = {{UINT64_C(0xfedcba9876543210), UINT64_C(0xfedcba9876543210) }};
  mzd_local_t cval, tmp;

  mzd_shift_left_uint64_128(&cval, &val, 1);
  tmp.w128[0] = mm128_shift_left(val.w128[0], 1);
  if (!mzd_local_equal(&cval, &tmp, 1, 128)) {
    printf("mm128 shift left fail: 1\n");
    ret = -1;
  }

  mzd_shift_left_uint64_128(&cval, &val, 2);
  tmp.w128[0] = mm128_shift_left(val.w128[0], 2);
  if (!mzd_local_equal(&cval, &tmp, 1, 128)) {
    printf("mm128 shift left fail: 2\n");
    ret = -1;
  }

  mzd_shift_right_uint64_128(&cval, &val, 1);
  tmp.w128[0] = mm128_shift_right(val.w128[0], 1);
  if (!mzd_local_equal(&cval, &tmp, 1, 128)) {
    printf("mm128 shift right fail: 1\n");
    ret = -1;
  }

  mzd_shift_right_uint64_128(&cval, &val, 2);
  tmp.w128[0] = mm128_shift_right(val.w128[0], 2);
  if (!mzd_local_equal(&cval, &tmp, 1, 128)) {
    printf("mm128 shift right fail: 2\n");
    ret = -1;
  }

  mzd_rotate_left_uint64_128(&cval, &val, 1);
  tmp.w128[0] = mm128_rotate_left(val.w128[0], 1);
  if (!mzd_local_equal(&cval, &tmp, 1, 128)) {
    printf("mm128 rotate left fail: 1\n");
    ret = -1;
  }

  mzd_rotate_left_uint64_128(&cval, &val, 2);
  tmp.w128[0] = mm128_rotate_left(val.w128[0], 2);
  if (!mzd_local_equal(&cval, &tmp, 1, 128)) {
    printf("mm128 rotate left fail: 2\n");
    ret = -1;
  }

  mzd_rotate_right_uint64_128(&cval, &val, 1);
  tmp.w128[0] = mm128_rotate_right(val.w128[0], 1);
  if (!mzd_local_equal(&cval, &tmp, 1, 128)) {
    printf("mm128 rotate right fail: 1\n");
    ret = -1;
  }

  mzd_rotate_right_uint64_128(&cval, &val, 2);
  tmp.w128[0] = mm128_rotate_right(val.w128[0], 2);
  if (!mzd_local_equal(&cval, &tmp, 1, 128)) {
    printf("mm128 rotate right fail: 2\n");
    ret = -1;
  }

  return ret;
}

ATTR_TARGET_S128
static int test_s128_256(void) {
  int ret = 0;

  mzd_local_t val = {{UINT64_C(0xfedcba9876543210), UINT64_C(0xfedcba9876543210),
                      UINT64_C(0xfedcba9876543210), UINT64_C(0xfedcba9876543210)}};
  mzd_local_t cval, tmp;

  mzd_shift_left_uint64_256(&cval, &val, 1);
  mm128_shift_left_256(tmp.w128, val.w128, 1);
  if (!mzd_local_equal(&cval, &tmp, 1, 256)) {
    printf("mm128_256 shift left fail: 1\n");
    ret = -1;
  }

  mzd_shift_left_uint64_256(&cval, &val, 2);
  mm128_shift_left_256(tmp.w128, val.w128, 2);
  if (!mzd_local_equal(&cval, &tmp, 1, 256)) {
    printf("mm128_256 shift left fail: 2\n");
    ret = -1;
  }

  mzd_shift_right_uint64_256(&cval, &val, 1);
  mm128_shift_right_256(tmp.w128, val.w128, 1);
  if (!mzd_local_equal(&cval, &tmp, 1, 256)) {
    printf("mm128_256 shift right fail: 1\n");
    ret = -1;
  }

  mzd_shift_right_uint64_256(&cval, &val, 2);
  mm128_shift_right_256(tmp.w128, val.w128, 2);
  if (!mzd_local_equal(&cval, &tmp, 1, 256)) {
    printf("mm128_256 shift right fail: 2\n");
    ret = -1;
  }

  mzd_rotate_left_uint64_256(&cval, &val, 1);
  mm128_rotate_left_256(tmp.w128, val.w128, 1);
  if (!mzd_local_equal(&cval, &tmp, 1, 256)) {
    printf("mm128_256 rotate left fail: 1\n");
    ret = -1;
  }

  mzd_rotate_left_uint64_256(&cval, &val, 2);
  mm128_rotate_left_256(tmp.w128, val.w128, 2);
  if (!mzd_local_equal(&cval, &tmp, 1, 256)) {
    printf("mm128_256 rotate left fail: 2\n");
    ret = -1;
  }

  mzd_rotate_right_uint64_256(&cval, &val, 1);
  mm128_rotate_right_256(tmp.w128, val.w128, 1);
  if (!mzd_local_equal(&cval, &tmp, 1, 256)) {
    printf("mm128_256 rotate right fail: 1\n");
    ret = -1;
  }

  mzd_rotate_right_uint64_256(&cval, &val, 2);
  mm128_rotate_right_256(tmp.w128, val.w128, 2);
  if (!mzd_local_equal(&cval, &tmp, 1, 256)) {
    printf("mm128_256 rotate right fail: 2\n");
    ret = -1;
  }

  return ret;
}
#endif /* WITH_SSE2 || WITH_NEON */

#if defined(WITH_AVX2)
ATTR_TARGET_AVX2
static int test_s256_256(void) {
  int ret = 0;

  mzd_local_t val = {{UINT64_C(0xfedcba9876543210), UINT64_C(0xfedcba9876543210),
                      UINT64_C(0xfedcba9876543210), UINT64_C(0xfedcba9876543210)}};
  mzd_local_t cval, tmp;

  mzd_shift_left_uint64_256(&cval, &val, 1);
  tmp.w256 = mm256_shift_left(val.w256, 1);
  if (!mzd_local_equal(&cval, &tmp, 1, 256)) {
    printf("mm256 shift left fail: 1\n");
    ret = -1;
  }

  mzd_shift_left_uint64_256(&cval, &val, 2);
  tmp.w256 = mm256_shift_left(val.w256, 2);
  if (!mzd_local_equal(&cval, &tmp, 1, 256)) {
    printf("mm256 shift left fail: 2\n");
    ret = -1;
  }

  mzd_shift_right_uint64_256(&cval, &val, 1);
  tmp.w256 = mm256_shift_right(val.w256, 1);
  if (!mzd_local_equal(&cval, &tmp, 1, 256)) {
    printf("mm256 shift right fail: 1\n");
    ret = -1;
  }

  mzd_shift_right_uint64_256(&cval, &val, 2);
  tmp.w256 = mm256_shift_right(val.w256, 2);
  if (!mzd_local_equal(&cval, &tmp, 1, 256)) {
    printf("mm256 shift right fail: 2\n");
    ret = -1;
  }

  mzd_rotate_left_uint64_256(&cval, &val, 1);
  tmp.w256 = mm256_rotate_left(val.w256, 1);
  if (!mzd_local_equal(&cval, &tmp, 1, 256)) {
    printf("mm256 rotate left fail: 1\n");
    ret = -1;
  }

  mzd_rotate_left_uint64_256(&cval, &val, 2);
  tmp.w256 = mm256_rotate_left(val.w256, 2);
  if (!mzd_local_equal(&cval, &tmp, 1, 256)) {
    printf("mm256 rotate left fail: 2\n");
    ret = -1;
  }

  mzd_rotate_right_uint64_256(&cval, &val, 1);
  tmp.w256 = mm256_rotate_right(val.w256, 1);
  if (!mzd_local_equal(&cval, &tmp, 1, 256)) {
    printf("mm256 rotate right fail: 1\n");
    ret = -1;
  }

  mzd_rotate_right_uint64_256(&cval, &val, 2);
  tmp.w256 = mm256_rotate_right(val.w256, 2);
  if (!mzd_local_equal(&cval, &tmp, 1, 256)) {
    printf("mm256 rotate right fail: 2\n");
    ret = -1;
  }

  return ret;
}
#endif /* WITH_AVX2 */
#endif /* WITH_OPT */

int main() {
  int ret = 0;

  ret |= test_uint64_128();
  ret |= test_uint64_256();
#if defined(WITH_OPT)
#if defined(WITH_SSE2) || defined(WITH_NEON)
  if (CPU_SUPPORTS_SSE2 || CPU_SUPPORTS_NEON) {
    ret |= test_s128_128();
    ret |= test_s128_256();
  }
#endif /* WITH_SSE2 || WITH_NEON */
#if defined(WITH_AVX2)
  if (CPU_SUPPORTS_AVX2) {
    ret |= test_s256_256();
  }
#endif /* WITH_AVX2 */
#endif

  return ret;
}
