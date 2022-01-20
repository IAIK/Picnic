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

#include <picnic.h>

#include "utils.h"

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <array>
#include <iostream>
#include <vector>

namespace {
  constexpr size_t rep         = 32;
  constexpr size_t message_inc = 512;

  uint8_t getBit(const uint8_t* array, size_t bitNumber) {
    return (array[bitNumber / 8] >> (7 - (bitNumber % 8))) & 0x01;
  }

  int picnic_test_keys(picnic_params_t parameters) {
    const size_t max_signature_size = picnic_signature_size(parameters);
    if (!max_signature_size) {
      // not supported
      return -2;
    }

    size_t diff            = 0;
    size_t lowmc_blocksize = 0;
    switch (parameters) {
    case Picnic_L1_full:
    case Picnic3_L1:
      diff            = 7;
      lowmc_blocksize = LOWMC_BLOCK_SIZE_Picnic_L1_full;
      break;
    case Picnic_L5_full:
    case Picnic3_L5:
      diff            = 1;
      lowmc_blocksize = LOWMC_BLOCK_SIZE_Picnic_L5_full;
      break;
    default:
      // instances with key size properly aligned
      return 0;
    }

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
      for (size_t i = 0; i < diff; ++i) {
        if (getBit(sk.data, 8 + lowmc_blocksize * 8 - i - 1)) {
          return -1;
        }
      }
      for (size_t i = 0; i < diff; ++i) {
        if (getBit(sk.data, 8 + 2 * lowmc_blocksize * 8 - i - 1) ||
            getBit(sk.data, 8 + 3 * lowmc_blocksize * 8 - i - 1)) {
          return -1;
        }
      }
    }

    return 0;
  }

  int picnic_test_multiple_message_sizes(picnic_params_t parameters) {
    const size_t max_signature_size = picnic_signature_size(parameters);
    if (!max_signature_size) {
      // not supported
      return -2;
    }

    picnic_publickey_t pk;
    picnic_privatekey_t sk;
    int ret = picnic_keygen(parameters, &pk, &sk);
    if (ret) {
      return ret;
    }

    std::vector<uint8_t> signature, message;
    for (size_t c = 0; c < rep; ++c) {
      const size_t message_size = (c + 1) * message_inc;

      signature.resize(max_signature_size);
      message.resize(message_size);
      // fill with random data
      randomize_container(message);
      randomize_container(signature);

      size_t signature_len = max_signature_size;
      ret = picnic_sign(&sk, message.data(), message_size, signature.data(), &signature_len);
      if (ret) {
        return ret;
      }

      // must verify
      ret = picnic_verify(&pk, message.data(), message_size, signature.data(), signature_len);
      if (ret) {
        return ret;
      }

      // must fail
      if (!picnic_verify(&pk, message.data(), message_size - 1, signature.data(), signature_len)) {
        return -1;
      }

      // must fail
      if (!picnic_verify(&pk, message.data(), message_size, signature.data(), signature_len - 1)) {
        return -1;
      }
    }

    return 0;
  }

  int picnic_test_with_read_write(picnic_params_t parameters) {
    const size_t max_signature_size = picnic_signature_size(parameters);
    if (!max_signature_size) {
      // not supported
      return -2;
    }

    picnic_publickey_t pk;
    picnic_privatekey_t sk;

    int ret = picnic_keygen(parameters, &pk, &sk);
    if (ret) {
      return ret;
    }

    std::vector<uint8_t> signature;
    signature.resize(max_signature_size);
    size_t signature_len = max_signature_size;

    std::array<uint8_t, 256> message;
    randomize_container(message);

    ret = picnic_sign(&sk, message.data(), message.size(), signature.data(), &signature_len);
    if (ret) {
      return ret;
    }

    ret = picnic_verify(&pk, message.data(), message.size(), signature.data(), signature_len);
    if (ret) {
      return ret;
    }

    std::vector<uint8_t> buf;
    buf.resize(picnic_get_public_key_size(parameters));
    randomize_container(buf);
    ret = picnic_write_public_key(&pk, buf.data(), buf.size());
    if (ret <= 0) {
      return ret;
    }

    picnic_publickey_t pk2;
    ret = picnic_read_public_key(&pk2, buf.data(), buf.size());
    if (ret) {
      return ret;
    }

    ret = picnic_verify(&pk2, message.data(), message.size(), signature.data(), signature_len);
    if (ret) {
      return ret;
    }

    buf.resize(picnic_get_private_key_size(parameters));
    randomize_container(buf);
    ret = picnic_write_private_key(&sk, buf.data(), buf.size());
    if (ret <= 0) {
      return ret;
    }

    picnic_privatekey_t sk2;
    ret = picnic_read_private_key(&sk2, buf.data(), buf.size());
    if (ret) {
      return ret;
    }

    ret = picnic_validate_keypair(&sk2, &pk2);
    if (ret) {
      return ret;
    }

    return 0;
  }

  int picnic_test_modified_signatures(picnic_params_t parameters) {
    const size_t max_signature_size = picnic_signature_size(parameters);
    if (!max_signature_size) {
      // not supported
      return -2;
    }

    picnic_publickey_t pk;
    picnic_privatekey_t sk;

    int ret = picnic_keygen(parameters, &pk, &sk);
    if (ret) {
      return ret;
    }

    std::vector<uint8_t> signature, message;
    signature.resize(max_signature_size);

    const size_t message_size = message_inc;
    message.resize(message_size);
    // fill message with some data
    randomize_container(message);

    size_t signature_len = max_signature_size;
    ret = picnic_sign(&sk, message.data(), message_size, signature.data(), &signature_len);
    if (ret) {
      return ret;
    }
    signature.resize(signature_len);

    // must verify
    ret = picnic_verify(&pk, message.data(), message_size, signature.data(), signature.size());
    if (ret) {
      return ret;
    }

    std::uniform_int_distribution<uint8_t> dist_diff{1};
    std::uniform_int_distribution<size_t> dist_pos{1, signature.size() / 8};
    std::random_device rnd;
    std::default_random_engine eng{rnd()};

    for (size_t l = dist_pos(eng), rounds = rep; rounds; l += dist_pos(eng), --rounds) {
      std::vector<uint8_t> signature2{signature};
      signature2[l % signature2.size()] += dist_diff(eng);

      // must fail
      if (!picnic_verify(&pk, message.data(), message_size, signature2.data(), signature2.size()))
        return -1;
    }

    return 0;
  }

  int perform_test(picnic_params_t param) {
    int ret = 0;

    std::cout << "testing keys: " << picnic_get_param_name(param) << " ...";
    int r = picnic_test_keys(param);
    if (r == -2) {
      std::cout << "SKIPPED";
    } else if (r) {
      std::cout << "FAILED";
      ret = -1;
    } else {
      std::cout << "OK";
    }

    std::cout << "\ntesting with read/write: " << picnic_get_param_name(param) << " ...";
    r = picnic_test_with_read_write(param);
    if (r == -2) {
      std::cout << "SKIPPED";
    } else if (r) {
      std::cout << "FAILED";
      ret = -1;
    } else {
      std::cout << "OK";
    }

    std::cout << "testing multiple message sizes: " << picnic_get_param_name(param) << " ...";
    r = picnic_test_multiple_message_sizes(param);
    if (r == -2) {
      std::cout << "SKIPPED";
    } else if (r) {
      std::cout << "FAILED";
      ret = -1;
    } else {
      std::cout << "OK";
    }

    std::cout << "testing with modified signatures " << picnic_get_param_name(param) << " ...";
    r = picnic_test_modified_signatures(param);
    if (r == -2) {
      std::cout << "SKIPPED" << std::endl;
    } else if (r) {
      std::cout << "FAILED" << std::endl;
      ret = -1;
    } else {
      std::cout << "OK" << std::endl;
    }
    return ret;
  }
} // namespace

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
