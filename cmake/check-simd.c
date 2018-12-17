#if defined(__GNUC__) || __has_attribute(target)
#define ATTRIBUTE_TARGET(x) __attribute__((target((x))))
#else
#define ATTRIBUTE_TARGET(x)
#endif

#if defined(SSE2) || defined(AVX2) || defined(BMI2)
#include <immintrin.h>

#if defined(SSE2)
ATTRIBUTE_TARGET("sse2") void test(void) {
  (void)_mm_setzero_si128();
}
#endif

#if defined(AVX2)
ATTRIBUTE_TARGET("avx2") void test(void) {
  (void)_mm256_setzero_si256();
}
#endif

#if defined(BMI2)
ATTRIBUTE_TARGET("bmi2") void test(void) {
  (void)_pext_u32(0, 0);
}
#endif
#endif

#if defined(NEON)
#include <arm_neon.h>

void test(void) {
  (void)vmovq_n_u64(0);
}
#endif

#if defined(POPCNT)
#include <nmmintrin.h>

ATTRIBUTE_TARGET("popcnt") void test(void) {
  (void)_mm_popcnt_u64(0);
}
#endif

int main() {
  test();
}
