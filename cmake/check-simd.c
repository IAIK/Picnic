#if defined(SSE2) || defined(SSE4_1) || defined(AVX2)
#include <immintrin.h>

#if defined(SSE2)
__attribute__((target("sse2"))) void test(void) {
  __m128i v = _mm_setzero_si128();
  (void)v;
}
#endif

#if defined(SSE4_1)
__attribute__((target("sse4.1"))) void test(void) {
  __m128i v = _mm_setzero_si128();
  (void)_mm_testz_si128(v, v);
}
#endif

#if defined(AVX2)
__attribute__((target("avx2"))) void test(void) {
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

int main() {
  test();
}
