#if defined(__GNUC__) || __has_attribute(target)
#define ATTRIBUTE_TARGET(x) __attribute__((target((x))))
#else
#define ATTRIBUTE_TARGET(x)
#endif

#if defined(SSE2) || defined(AVX2)
#include <immintrin.h>

#if defined(SSE2)
ATTRIBUTE_TARGET("sse2") void test(void) {
  __m128i v = _mm_setzero_si128();
  (void)v;
}
#endif

#if defined(AVX2)
ATTRIBUTE_TARGET("avx2") void test(void) {
  __m256i v = _mm256_setzero_si256();
  (void)v;
}
#endif

#endif

#if defined(NEON)
#include <arm_neon.h>

void test(void) {
  uint32x4_t v = vmovq_n_u32(0);
  (void)v;
}
#endif

#if defined(SSE42)
#include <nmmintrin.h>

ATTRIBUTE_TARGET("sse4.2") void test(void) {
  unsigned long long v = _mm_popcnt_u64(0);
  (void)v;
}
#endif

int main() {
  test();
}
