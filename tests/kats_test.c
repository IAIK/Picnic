#include "picnic.h"

#include <memory.h>
#include <stdio.h>
#include <stdlib.h>

static int run_picnic_test(const uint8_t* msg, size_t msg_len, const uint8_t* pk, size_t pk_len,
                           const uint8_t* sk, size_t sk_len, const uint8_t* sig, size_t sig_len) {
  picnic_privatekey_t private_key;
  picnic_publickey_t public_key;
  size_t signature_len = sig_len + 5000;

  uint8_t* signature = malloc(signature_len);

  int ret = picnic_read_private_key(&private_key, sk, sk_len);
  if (ret != 0) {
    return 0;
  }

  ret = picnic_read_public_key(&public_key, pk, pk_len);
  if (ret != 0) {
    return 0;
  }

  ret = picnic_validate_keypair(&private_key, &public_key);
  if (ret != 0) {
    return 0;
  }

  /* Recreate the signature, check it matches */
  ret = picnic_sign(&private_key, msg, msg_len, signature, &signature_len);
  if (ret != 0) {
    return 0;
  }

  if (signature_len != sig_len) {
    return 0;
  }
  if (memcmp(sig, signature, signature_len) != 0) {
    return 0;
  }

  /* Verify the provided signature */
  ret = picnic_verify(&public_key, msg, msg_len, sig, sig_len);
  if (ret != 0) {
    return 0;
  }

  free(signature);
  return 1;
}

/* These are the KATs from the submission, the first in each rsp file */

static int picnic_test_vector_L1FS(void) {
#include "kat_l1_fs.c.i"
  return run_picnic_test(msg, sizeof(msg), pk, sizeof(pk), sk, sizeof(sk), sig, sizeof(sig));
}

static int picnic_test_vector_L1UR(void) {
#include "kat_l1_ur.c.i"
  return run_picnic_test(msg, sizeof(msg), pk, sizeof(pk), sk, sizeof(sk), sig, sizeof(sig));
}

static int picnic_test_vector_L3FS(void) {
#include "kat_l3_fs.c.i"
  return run_picnic_test(msg, sizeof(msg), pk, sizeof(pk), sk, sizeof(sk), sig, sizeof(sig));
}

static int picnic_test_vector_L3UR(void) {
#include "kat_l3_ur.c.i"
  return run_picnic_test(msg, sizeof(msg), pk, sizeof(pk), sk, sizeof(sk), sig, sizeof(sig));
}

static int picnic_test_vector_L5FS(void) {
#include "kat_l5_fs.c.i"
  return run_picnic_test(msg, sizeof(msg), pk, sizeof(pk), sk, sizeof(sk), sig, sizeof(sig));
}

static int picnic_test_vector_L5UR(void) {
#include "kat_l5_ur.c.i"
  return run_picnic_test(msg, sizeof(msg), pk, sizeof(pk), sk, sizeof(sk), sig, sizeof(sig));
}

typedef int (*test_fn_t)(void);

static const test_fn_t tests[] = {picnic_test_vector_L1FS, picnic_test_vector_L1UR,
                                  picnic_test_vector_L3FS, picnic_test_vector_L3UR,
                                  picnic_test_vector_L5FS, picnic_test_vector_L5UR};

static const size_t num_tests = sizeof(tests) / sizeof(tests[0]);

int main() {
  int ret = 0;
  for (size_t s = 0; s < num_tests; ++s) {
    const int t = tests[s]();
    if (!t) {
      printf("ERR: Picnic KAT test %zu FAILED (%d)\n", s, t);
      ret = -1;
    }
  }

  return ret;
}
