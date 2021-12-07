#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "../picnic_impl.h"
#include "../io.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int test_vector(const picnic_params_t param) {
  if (param <= PARAMETER_SET_INVALID || param >= PARAMETER_SET_MAX_INDEX) {
    printf("invalid parameter set\n");
    return -1;
  }

  printf("Picnic test vectors for parameter set: %s\n", picnic_get_param_name(param));

  const picnic_instance_t* instance = picnic_instance_get(param);
  if (!instance) {
    printf("unsupported instance\n");
    return -1;
  }

  picnic_privatekey_t sk = {{0}};
  sk.data[0]             = param;
  for (unsigned int i = 0; i < instance->lowmc.k; i += 2) {
    setBit(&sk.data[1], i, 1);
  }
  for (unsigned int i = 0; i < instance->lowmc.n; i += 3) {
    setBit(&sk.data[1 + instance->input_output_size + instance->input_output_size], i, 1);
  }

  picnic_publickey_t pk = {{0}};
  picnic_sk_to_pk(&sk, &pk);
  memcpy(&sk.data[1 + instance->input_output_size], &pk.data[1], instance->input_output_size);

  const uint8_t msg[] = {0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09};

  if (picnic_validate_keypair(&sk, &pk)) {
    return -1;
  }
  picnic_visualize_keys(stdout, &sk, &pk);

  const size_t max_signature_size = picnic_signature_size(param);
  if (!max_signature_size) {
    return -1;
  }

  uint8_t* sig  = malloc(max_signature_size);
  size_t siglen = max_signature_size;
  int ret       = 0;
  if (!picnic_sign(&sk, msg, sizeof(msg), sig, &siglen)) {
    picnic_visualize(stdout, &pk, msg, sizeof(msg), sig, siglen);

    if (picnic_verify(&pk, msg, sizeof(msg), sig, siglen)) {
      ret = -1;
      printf("verify:  failed\n");
    } else {
      printf("verify:  success\n");
    }
  } else {
    ret = -1;
  }

  free(sig);
  return ret;
}

int main(int argc, char** argv) {
  if (argc != 2) {
    printf("provide an integer specifying the parameter set\n");
    return -1;
  }

  return test_vector(atoi(argv[1]));
}
