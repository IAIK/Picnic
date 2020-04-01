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
#include "picnic_instances.h"
#include "io.h"
#include "utils.h"

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const size_t rep         = 32;
static const size_t message_inc = 512;

static int picnic_test_keys(picnic_params_t parameters) {
  const picnic_instance_t* instance = picnic_instance_get(parameters);
  if (!instance) {
    /* not supported */
    return -2;
  }

  const size_t diff_key = instance->input_size * 8 - instance->lowmc.k;
  const size_t diff_block = instance->output_size * 8 - instance->lowmc.n;
  /* instances with key size properly aligned */
  if (!diff_key && !diff_block) {
    return -2;
  }

  const size_t lowmc_blocksize = picnic_get_lowmc_block_size(parameters);

  for (size_t c = 0; c < rep; ++c) {
    picnic_publickey_t pk  = {0};
    picnic_privatekey_t sk = {0};

    int ret = picnic_keygen(parameters, &pk, &sk);
    if (ret) {
      return ret;
    }

    // Public and private keys are serialized as follows:
    // - public key: instance || C || p
    // - secret key: instance || sk || C || p

    for (size_t i = 0; i < diff_key; ++i) {
      if (getBit(sk.data, 8 + lowmc_blocksize * 8 - i - 1)) {
        return -1;
      }
    }
    for (size_t i = 0; i < diff_block; ++i) {
      if (getBit(sk.data, 8 + 2 * lowmc_blocksize * 8 - i - 1)) {
        return -1;
      }
      if (getBit(sk.data, 8 + 3 * lowmc_blocksize * 8 - i - 1)) {
        return -1;
      }
    }
  }

  return 0;
}

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

static int picnic_test_modified_signatures(picnic_params_t parameters) {
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

  const size_t message_size = message_inc;
  /* fill message with some data */
  memset(message, rand() & 0xFF, message_size);

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

  for (size_t l = 0; l < signature_len; ++l) {
    signature[l] += 1;

    /* must fail */
    ret = picnic_verify(&pk, message, message_size - 1, signature, signature_len);
    if (!ret) {
      ret = -1;
      goto end;
    }

    signature[l] -= 1;
  }

  ret = 0;

end:
  free(message);
  free(signature);
  return ret;
}

static int perform_test(picnic_params_t param) {
  int ret = 0;

  printf("testing keys: %s ...", picnic_get_param_name(param));
  int r   = picnic_test_keys(param);
  if (r == -2) {
    printf("SKIPPED\n");
  } else if (r) {
    printf("FAILED\n");
    ret = -1;
  } else {
    printf("OK\n");
  }

  printf("testing with read/write: %s ... ", picnic_get_param_name(param));
  r = picnic_test_with_read_write(param);
  if (r == -2) {
    printf("SKIPPED\n");
  } else if (r) {
    printf("FAILED\n");
    ret = -1;
  } else {
    printf("OK\n");
  }

  printf("testing multiple message sizes: %s ... ", picnic_get_param_name(param));
  r = picnic_test_multiple_message_sizes(param);
  if (r == -2) {
    printf("SKIPPED\n");
  } else if (r) {
    printf("FAILED\n");
    ret = -1;
  } else {
    printf("OK\n");
  }
  return ret;

  printf("testing with modified signatures: %s ... ", picnic_get_param_name(param));
  r = picnic_test_modified_signatures(param);
  if (r == -2) {
    printf("SKIPPED\n");
  } else if (r) {
    printf("FAILED\n");
    ret = -1;
  } else {
    printf("OK\n");
  }
  return ret;
}

int main(int argc, char** argv) {
  if (argc == 2) {
    const picnic_params_t param = argument_to_params(argv[1], true);
    if (param == PARAMETER_SET_INVALID) {
      printf("ERR: invalid test idx\n");
      return 1;
    }

    return perform_test(param);
  }

  int ret = 0;
  for (unsigned int param = Picnic_L1_FS; param < PARAMETER_SET_MAX_INDEX; ++param) {
    const int r = perform_test(param);
    if (r) {
      ret = -1;
    }
  }

  return ret;
}
