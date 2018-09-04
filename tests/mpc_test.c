#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <m4ri/m4ri.h>

#include "../mpc.h"
#include "../mzd_additional.h"

#include "utils.h"

static mzd_local_t** mpc_init_empty_share_vector(uint32_t n, unsigned sc) {
  mzd_local_t** s = malloc(sc * sizeof(mzd_local_t*));
  mzd_local_init_multiple(s, sc, 1, n);
  return s;
}

static mzd_local_t* mpc_reconstruct_from_share(mzd_local_t* dst, mzd_local_t** shared_vec) {
  if (!dst) {
    dst = mzd_local_init(shared_vec[0]->nrows, shared_vec[0]->ncols);
  }

  mzd_xor(dst, shared_vec[0], shared_vec[1]);
  mzd_xor(dst, dst, shared_vec[2]);
  return dst;
}

static mzd_local_t* mzd_init_random_vector(rci_t n) {
  mzd_local_t* a = mzd_local_init(1, n);
  mzd_randomize_ssl(a);
  return a;
}

static mzd_local_t** mpc_init_share_vector(mzd_local_t const* v) {
  mzd_local_t** s = malloc(3 * sizeof(mzd_local_t*));
  mzd_local_init_multiple(s, 3, 1, v->ncols);

  mzd_randomize_ssl(s[0]);
  mzd_randomize_ssl(s[1]);

  mzd_xor(s[2], s[0], s[1]);
  mzd_xor(s[2], s[2], v);

  return s;
}

static void test_mpc_share(void) {
  mzd_local_t* t1    = mzd_init_random_vector(10);
  mzd_local_t** s1   = mpc_init_share_vector(t1);
  mzd_local_t* t1cmb = mpc_reconstruct_from_share(NULL, s1);

  if (mzd_local_equal(t1, t1cmb))
    printf("Share test successful.\n");

  mzd_local_free(t1);
  mzd_local_free_multiple(s1);
  mzd_local_free(t1cmb);
  free(s1);
}

static void test_mpc_add(void) {
  mzd_local_t* t1  = mzd_init_random_vector(10);
  mzd_local_t* t2  = mzd_init_random_vector(10);
  mzd_local_t* res = mzd_local_init(1, 10);
  mzd_xor(res, t1, t2);

  mzd_local_t** s1   = mpc_init_share_vector(t1);
  mzd_local_t** s2   = mpc_init_share_vector(t2);
  mzd_local_t** ress = mpc_init_empty_share_vector(10, 3);
  mpc_xor(ress, s1, s2, 3);

  mzd_local_t* cmp = mpc_reconstruct_from_share(NULL, ress);

  if (mzd_local_equal(res, cmp))
    printf("Shared add test successful.\n");

  mzd_local_free(t1);
  mzd_local_free(t2);
  mzd_local_free(res);
  mzd_local_free_multiple(s1);
  mzd_local_free_multiple(s2);
  mzd_local_free_multiple(ress);
  mzd_local_free(cmp);
  free(s1);
  free(s2);
  free(ress);
}

void run_tests(void) {
  test_mpc_share();
  test_mpc_add();
}

int main() {
  run_tests();

  return 0;
}
