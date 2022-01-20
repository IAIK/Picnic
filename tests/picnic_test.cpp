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

#include <iostream>
#include <vector>

#include <picnic.h>
#include "utils.h"

static int picnic_sign_verify(const picnic_params_t param) {
  static constexpr uint8_t m[] = "test message";

  const size_t max_signature_size = picnic_signature_size(param);
  if (!max_signature_size) {
    // not supported
    return -2;
  }

  // Create a key pair
  picnic_privatekey_t private_key;
  picnic_publickey_t public_key;
  std::cout << "Creating key pair ... ";
  if (picnic_keygen(param, &public_key, &private_key)) {
    std::cout << "FAILED!" << std::endl;
    return -1;
  }
  std::cout << "OK" << std::endl;

  // Valid key pair
  std::cout << "Validating key pair ... ";
  if (picnic_get_public_key_param(&public_key) != param ||
      picnic_get_private_key_param(&private_key) != param ||
      picnic_validate_keypair(&private_key, &public_key)) {
    std::cout << "FAILED!" << std::endl;
    return -1;
  }
  std::cout << "OK" << std::endl;

  std::vector<uint8_t> sig;
  sig.resize(max_signature_size);
  size_t siglen = max_signature_size;

  // Sign a message
  std::cout << "Signing message ... ";
  if (picnic_sign(&private_key, m, sizeof(m), sig.data(), &siglen)) {
    std::cout << "FAILED!" << std::endl;
    return -1;
  }
  // Verify signature
  std::cout << "OK\nVerifying signature ... ";
  if (picnic_verify(&public_key, m, sizeof(m), sig.data(), siglen)) {
    std::cout << "FAILED!" << std::endl;
    return -1;
  }
  std::cout << "OK" << std::endl;

  return 0;
}

static int perform_test(const picnic_params_t param) {
  std::cout << "testing: " << picnic_get_param_name(param) << " ...";
  const int r = picnic_sign_verify(param);
  if (r == -2) {
    std::cout << "SKIPPED\n";
  } else if (r) {
    std::cout << "FAILED\n";
    return -1;
  } else {
    std::cout << "OK\n";
  }
  return 0;
}

int main(int argc, char** argv) {
  if (argc == 2) {
    const picnic_params_t param = argument_to_params(argv[1]);
    if (param == PARAMETER_SET_INVALID) {
      std::cout << "ERR: invalid test idx" << std::endl;
      return 1;
    }

    return perform_test(param);
  }

  int ret = 0;
  for (const auto param : all_parameters) {
    const int r = perform_test(param);
    if (r) {
      ret = -1;
    }
  }

  return ret;
}
