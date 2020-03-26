/*
 *  This file is part of the optimized implementation of the Picnic signature scheme.
 *  See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */

#if defined(HAVE_CONFIG_H)
#include <config.h>
#else
#error "The KATs require config.h to be available."
#endif

#include "picnic.h"
#include "macros.h"
#include "utils.h"

#include <memory.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>

typedef struct {
  size_t mlen;
  uint8_t* msg;
  uint8_t pk[PICNIC_MAX_PUBLICKEY_SIZE];
  uint8_t sk[PICNIC_MAX_PRIVATEKEY_SIZE];
  size_t smlen;
  uint8_t* sm;
} test_vector_t;

static void clear_test_vector(test_vector_t* tv) {
  free(tv->msg);
  free(tv->sm);
  memset(tv, 0, sizeof(*tv));
}

static uint8_t parse_hex_c(const char c) {
  if (c >= '0' && c <= '9') {
    return c - '0';
  } else if (c >= 'a' && c <= 'f') {
    return 10 + c - 'a';
  } else if (c >= 'A' && c <= 'F') {
    return 10 + c - 'A';
  } else {
    return UINT8_MAX;
  }
}

static int parse_hex(uint8_t* dst, const char* src, size_t len) {
  for (size_t s = 0; s < len; ++s, src += 2, ++dst) {
    uint8_t high = parse_hex_c(src[0]);
    uint8_t low  = parse_hex_c(src[1]);
    if (high == UINT8_MAX || low == UINT8_MAX) {
      printf("parse_hex failed\n");
      return -1;
    }
    *dst = high << 4 | low;
  }
  return 0;
}

#if !defined(HAVE_GETLINE)
#if !defined(HAVE_SSIZE_T)
// this should be good enough
typedef int64_t ssize_t;
#endif

static ssize_t getline(char** line, size_t* len, FILE* fp) {
  if (line == NULL || len == NULL || fp == NULL) {
    errno = EINVAL;
    return -1;
  }

  char chunk[4096];
  if (*line == NULL || *len < sizeof(chunk)) {
    *len  = sizeof(chunk);
    *line = malloc(*len);
    if (*line == NULL) {
      errno = ENOMEM;
      return -1;
    }
  }

  (*line)[0] = '\0';
  while (fgets(chunk, sizeof(chunk), fp)) {
    size_t len_used         = strlen(*line);
    const size_t chunk_used = strlen(chunk);

    if (*len - len_used < chunk_used) {
      if (*len > SIZE_MAX / 2) {
        errno = EOVERFLOW;
        return -1;
      } else {
        *len *= 2;
      }

      char* tmp_line = realloc(*line, *len);
      if (tmp_line == NULL) {
        errno = ENOMEM;
        return -1;
      }
      *line = tmp_line;
    }

    memcpy(*line + len_used, chunk, chunk_used);
    len_used += chunk_used;
    (*line)[len_used] = '\0';

    if ((*line)[len_used - 1] == '\n') {
      return len_used;
    }
  }

  return -1;
}
#endif

static int read_test_vector(FILE* file, test_vector_t* tv, size_t pks, size_t sks) {
  char* line = NULL;
  size_t len = 0;
  ssize_t nread;
  bool expect_data = false;

  while ((nread = getline(&line, &len, file)) != -1) {
    if (nread <= 1 || line[0] == '#' || (nread == 2 && line[0] == '\r' && line[1] == '\n')) {
      if (expect_data) {
        printf("Expected data.\n");
        goto err;
      }
      // skip empty lines and comments
      continue;
    }

    const size_t uread = nread;
    if (strncmp(line, "count = ", 8) == 0) {
      // skip count
      expect_data = true;
      continue;
    } else if (strncmp(line, "seed = ", 7) == 0) {
      // skip seed
      continue;
    } else if (strncmp(line, "mlen = ", 7) == 0) {
      // read message length
      if (sscanf(line + 7, SIZET_FMT, &tv->mlen) != 1) {
        printf("Unable to parse message length.\n");
        goto err;
      }
    } else if (strncmp(line, "msg = ", 6) == 0 && tv->mlen && uread >= 2 * tv->mlen + 6) {
      // read message
      tv->msg = calloc(1, tv->mlen);
      if (parse_hex(tv->msg, line + 6, tv->mlen) == -1) {
        printf("Unable to parse message.\n");
        goto err;
      }
    } else if (strncmp(line, "pk = ", 5) == 0 && uread >= 2 * pks + 5) {
      // read pk
      if (parse_hex(tv->pk, line + 5, pks) == -1) {
        printf("Unable to parse public key.\n");
        goto err;
      }
    } else if (strncmp(line, "sk = ", 5) == 0 && uread >= 2 * sks + 5) {
      // read sk
      if (parse_hex(tv->sk, line + 5, sks) == -1) {
        printf("Unable to parse secret key.\n");
        goto err;
      }
    } else if (strncmp(line, "smlen = ", 8) == 0) {
      // read signature length
      if (sscanf(line + 8, SIZET_FMT, &tv->smlen) != 1) {
        printf("Unable to parse signature length.\n");
        goto err;
      }
    } else if (strncmp(line, "sm = ", 5) == 0 && tv->smlen && uread >= 2 * tv->smlen + 5) {
      // read signature
      tv->sm = calloc(1, tv->smlen);
      if (parse_hex(tv->sm, line + 5, tv->smlen) == -1) {
        printf("Unable to parse signature.\n");
        goto err;
      }
      break;
    } else {
      printf("Do not know how handle line (len = " SIZET_FMT "): %s", uread, line);
      goto err;
    }
  }

  free(line);
  line = NULL;

  if (!tv->mlen || !tv->smlen || !tv->msg || !tv->sm) {
    if (expect_data) {
      printf("Missing test vector data (mlen: " SIZET_FMT " smlen: " SIZET_FMT ").\n", tv->mlen,
             tv->smlen);
    }
    goto err;
  }

  return 0;

err:
  free(line);
  clear_test_vector(tv);
  return -1;
}

static int run_picnic_test(const uint8_t* msg, size_t msg_len, const uint8_t* pk, size_t pk_len,
                           const uint8_t* sk, size_t sk_len, const uint8_t* sig, size_t sig_len) {
  picnic_privatekey_t private_key;
  picnic_publickey_t public_key;
  size_t signature_len = sig_len + 5000;

  uint8_t* signature = malloc(signature_len);

  int ret = picnic_read_private_key(&private_key, sk, sk_len);
  if (ret != 0) {
    printf("Unable to read private key.\n");
    goto err;
  }

  ret = picnic_read_public_key(&public_key, pk, pk_len);
  if (ret != 0) {
    printf("Unable to read public key.\n");
    goto err;
  }

  ret = picnic_validate_keypair(&private_key, &public_key);
  if (ret != 0) {
    printf("Key pair does not validate.\n");
    goto err;
  }

  /* Recreate the signature, check it matches */
  ret = picnic_sign(&private_key, msg, msg_len, signature, &signature_len);
  if (ret != 0) {
    printf("Unable to sign.\n");
    goto err;
  }

  if (signature_len != sig_len) {
    printf("Signature length does not match.\n");
    goto err;
  }
  if (memcmp(sig, signature, signature_len) != 0) {
    printf("Signature does not match.\n");
    goto err;
  }

  /* Verify the provided signature */
  ret = picnic_verify(&public_key, msg, msg_len, sig, sig_len);
  if (ret != 0) {
    printf("Signature does not verify.\n");
    goto err;
  }

  free(signature);
  return 1;

err:
  free(signature);
  return 0;
}

static int run_test_vectors_from_file(const char* path, size_t pks, size_t sks) {
  FILE* file = fopen(path, "r");
  if (!file) {
    return 0;
  }

  size_t vectors_run       = 0;
  size_t vectors_succeeded = 0;
  test_vector_t tv         = {0, NULL, {0}, {0}, 0, NULL};
  while (read_test_vector(file, &tv, pks, sks) != -1) {
    // Test vectors generated for NIST have message length and the message at the beginning.
    const size_t offset = tv.mlen + sizeof(uint32_t);

    ++vectors_run;
    vectors_succeeded +=
        run_picnic_test(tv.msg, tv.mlen, tv.pk, pks, tv.sk, sks, tv.sm + offset, tv.smlen - offset);
    clear_test_vector(&tv);
  };
  fclose(file);

  return (vectors_run && vectors_succeeded == vectors_run) ? 1 : 0;
}

typedef struct {
  const char* path;
  size_t pks;
  size_t sks;
} test_t;

static int run_test(picnic_params_t param) {
  static const test_t tests[] = {
      {NULL, 0, 0},
      {KATDIR "/kat_l1_fs.txt", PICNIC_PUBLIC_KEY_SIZE(Picnic_L1_FS),
       PICNIC_PRIVATE_KEY_SIZE(Picnic_L1_FS)},
      {KATDIR "/kat_l1_ur.txt", PICNIC_PUBLIC_KEY_SIZE(Picnic_L1_UR),
       PICNIC_PRIVATE_KEY_SIZE(Picnic_L1_FS)},
      {KATDIR "/kat_l3_fs.txt", PICNIC_PUBLIC_KEY_SIZE(Picnic_L3_FS),
       PICNIC_PRIVATE_KEY_SIZE(Picnic_L3_FS)},
      {KATDIR "/kat_l3_ur.txt", PICNIC_PUBLIC_KEY_SIZE(Picnic_L3_UR),
       PICNIC_PRIVATE_KEY_SIZE(Picnic_L3_FS)},
      {KATDIR "/kat_l5_fs.txt", PICNIC_PUBLIC_KEY_SIZE(Picnic_L5_FS),
       PICNIC_PRIVATE_KEY_SIZE(Picnic_L5_FS)},
      {KATDIR "/kat_l5_ur.txt", PICNIC_PUBLIC_KEY_SIZE(Picnic_L5_UR),
       PICNIC_PRIVATE_KEY_SIZE(Picnic_L5_FS)},
      {KATDIR "/kat_picnic2_l1_fs.txt", PICNIC_PUBLIC_KEY_SIZE(Picnic3_L1_FS),
       PICNIC_PRIVATE_KEY_SIZE(Picnic3_L1_FS)},
      {KATDIR "/kat_picnic2_l3_fs.txt", PICNIC_PUBLIC_KEY_SIZE(Picnic3_L3_FS),
       PICNIC_PRIVATE_KEY_SIZE(Picnic3_L3_FS)},
      {KATDIR "/kat_picnic2_l5_fs.txt", PICNIC_PUBLIC_KEY_SIZE(Picnic3_L5_FS),
       PICNIC_PRIVATE_KEY_SIZE(Picnic3_L5_FS)}};

  return run_test_vectors_from_file(tests[param].path, tests[param].pks, tests[param].sks);
}

int main(int argc, char** argv) {
  if (argc == 2) {
    const picnic_params_t param = argument_to_params(argv[1], false);
    if (param == PARAMETER_SET_INVALID) {
      printf("ERR: invalid test idx\n");
      return 1;
    }

    const int t = run_test(param);
    if (!t) {
      printf("ERR: Picnic KAT test %u FAILED (%d)\n", param, t);
      return -1;
    }
    return 0;
  }

  int ret = 0;
  for (size_t s = Picnic_L1_FS; s <= Picnic3_L5_FS; ++s) {
    const int t = run_test(s);
    if (!t) {
      printf("ERR: Picnic KAT test " SIZET_FMT " FAILED (%d)\n", s, t);
      ret = -1;
    }
  }

  return ret;
}
