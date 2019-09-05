/*
 *  This file is part of the optimized implementation of the Picnic signature scheme.
 *  See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "picnic.h"

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const size_t rep         = 32;
static const size_t message_inc = 512;

static int picnic_test_multiple_message_sizes(picnic_params_t parameters) {
  const size_t max_signature_size = picnic_signature_size(parameters);
  if (!max_signature_size) {
    /* not supported */
    return -2;
  }

  picnic_publickey_t pk;
  picnic_privatekey_t sk;

  int ret = picnic_keygen(parameters, &pk, &sk);
  if (ret) {
    return ret;
  }

  uint8_t* signature = malloc(max_signature_size);
  uint8_t* message   = malloc(rep * message_inc);
  if (!signature || !message) {
    ret = -1;
    goto end;
  }

  for (size_t c = 0; c < rep; ++c) {
    const size_t message_size = (c + 1) * message_inc;
    /* fill message with some data */
    memset(message, rand() & 0xFF, message_size);
    memset(signature, rand() & 0xFF, max_signature_size);

    size_t signature_len = max_signature_size;
    ret                  = picnic_sign(&sk, message, message_size, signature, &signature_len);
    if (ret) {
      goto end;
    }

    /* must verify */
    ret = picnic_verify(&pk, message, message_size, signature, signature_len);
    if (ret) {
      goto end;
    }

    /* must fail */
    ret = picnic_verify(&pk, message, message_size - 1, signature, signature_len);
    if (!ret) {
      ret = -1;
      goto end;
    }

    /* must fail */
    ret = picnic_verify(&pk, message, message_size, signature, signature_len - 1);
    if (!ret) {
      ret = -1;
      goto end;
    }
  }

  ret = 0;

end:
  free(message);
  free(signature);
  return ret;
}

static int picnic_test_with_read_write(picnic_params_t parameters) {
  const size_t max_signature_size = picnic_signature_size(parameters);
  if (!max_signature_size) {
    /* not supported */
    return -2;
  }

  picnic_publickey_t pk;
  picnic_privatekey_t sk;

  int ret = picnic_keygen(parameters, &pk, &sk);
  if (ret) {
    return ret;
  }

  size_t signature_len = max_signature_size;
  uint8_t* signature   = malloc(max_signature_size);
  if (!signature) {
    return -1;
  }

  uint8_t message[256];
  memset(message, 0x12, sizeof(message));

  ret = picnic_sign(&sk, message, sizeof(message), signature, &signature_len);
  if (ret) {
    return ret;
  }

  ret = picnic_verify(&pk, message, sizeof(message), signature, signature_len);
  if (ret) {
    return ret;
  }

  uint8_t pk_buf[PICNIC_MAX_PUBLICKEY_SIZE + 1];
  ret = picnic_write_public_key(&pk, pk_buf, sizeof(pk_buf));
  if (ret <= 0) {
    return ret;
  }

  memset(&pk, 0x00, sizeof(picnic_publickey_t));

  ret = picnic_read_public_key(&pk, pk_buf, sizeof(pk_buf));
  if (ret) {
    return ret;
  }

  ret = picnic_verify(&pk, message, sizeof(message), signature, signature_len);
  if (ret) {
    return ret;
  }

  uint8_t sk_buf[PICNIC_MAX_PRIVATEKEY_SIZE + 1];
  ret = picnic_write_private_key(&sk, sk_buf, sizeof(sk_buf));
  if (ret <= 0) {
    return ret;
  }

  memset(&sk, 0x00, sizeof(picnic_privatekey_t));
  ret = picnic_read_private_key(&sk, sk_buf, sizeof(sk_buf));
  if (ret) {
    return ret;
  }

  ret = picnic_validate_keypair(&sk, &pk);
  if (ret) {
    return ret;
  }

  free(signature);
  return 0;
}

int main(void) {
  int ret = 0;
  for (picnic_params_t params = 1; params < PARAMETER_SET_MAX_INDEX; params++) {
    printf("testing with read/write: %s ... ", picnic_get_param_name(params));
    int r = picnic_test_with_read_write(params);
    if (r == -2) {
      printf("SKIPPED\n");
    } else if (r) {
      printf("FAILED\n");
      ret = -1;
    } else {
      printf("OK\n");
    }

    printf("testing multiple message sizes: %s ... ", picnic_get_param_name(params));
    r = picnic_test_multiple_message_sizes(params);
    if (r == -2) {
      printf("SKIPPED\n");
    } else if (r) {
      printf("FAILED\n");
      ret = -1;
    } else {
      printf("OK\n");
    }
  }
  return ret;
}
