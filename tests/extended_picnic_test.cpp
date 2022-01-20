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

#include <array>
#include <vector>

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

namespace {
  constexpr size_t rep         = 32;
  constexpr size_t message_inc = 512;

  uint8_t getBit(const uint8_t* array, size_t bitNumber) {
    return (array[bitNumber / 8] >> (7 - (bitNumber % 8))) & 0x01;
  }
} // namespace

BOOST_DATA_TEST_CASE(test_keys, all_supported_parameters(), parameters) {
  BOOST_TEST_CONTEXT("Parameter set: " << picnic_get_param_name(parameters)) {
    const size_t max_signature_size = picnic_signature_size(parameters);
    BOOST_TEST(max_signature_size);

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
      return;
    }

    for (size_t c = 0; c < rep; ++c) {
      picnic_publickey_t pk  = {0};
      picnic_privatekey_t sk = {0};

      BOOST_TEST(!picnic_keygen(parameters, &pk, &sk));

      // Public and private keys are serialized as follows:
      // - public key: instance || C || p
      // - secret key: instance || sk || C || p
      for (size_t i = 0; i < diff; ++i) {
        BOOST_TEST_INFO("sk bit " << (8 + lowmc_blocksize * 8 - i - 1));
        BOOST_TEST(!getBit(sk.data, 8 + lowmc_blocksize * 8 - i - 1));
      }
      for (size_t i = 0; i < diff; ++i) {
        BOOST_TEST_INFO("sk bit " << (8 + 2 * lowmc_blocksize * 8 - i - 1));
        BOOST_TEST(!getBit(sk.data, 8 + 2 * lowmc_blocksize * 8 - i - 1));
        BOOST_TEST_INFO("sk bit " << (8 + 3 * lowmc_blocksize * 8 - i - 1));
        BOOST_TEST(!getBit(sk.data, 8 + 3 * lowmc_blocksize * 8 - i - 1));
      }
    }
  }
}

BOOST_DATA_TEST_CASE(multiple_messages, all_supported_parameters(), parameters) {
  BOOST_TEST_CONTEXT("Parameter set: " << picnic_get_param_name(parameters)) {
    const size_t max_signature_size = picnic_signature_size(parameters);
    BOOST_TEST(max_signature_size);

    picnic_publickey_t pk;
    picnic_privatekey_t sk;
    BOOST_TEST(!picnic_keygen(parameters, &pk, &sk));

    std::vector<uint8_t> signature, message;
    for (size_t c = 0; c < rep; ++c) {
      const size_t message_size = (c + 1) * message_inc;

      signature.resize(max_signature_size);
      message.resize(message_size);
      // fill with random data
      randomize_container(message);
      randomize_container(signature);

      size_t signature_len = max_signature_size;
      BOOST_TEST(!picnic_sign(&sk, message.data(), message_size, signature.data(), &signature_len));
      BOOST_TEST(
          !picnic_verify(&pk, message.data(), message_size, signature.data(), signature_len));
      // must fail
      BOOST_TEST(
          picnic_verify(&pk, message.data(), message_size - 1, signature.data(), signature_len));
      // must fail
      BOOST_TEST(
          picnic_verify(&pk, message.data(), message_size, signature.data(), signature_len - 1));
    }
  }
}

BOOST_DATA_TEST_CASE(read_write, all_supported_parameters(), parameters) {
  BOOST_TEST_CONTEXT("Parameter set: " << picnic_get_param_name(parameters)) {
    const size_t max_signature_size = picnic_signature_size(parameters);
    BOOST_TEST(max_signature_size);

    picnic_publickey_t pk;
    picnic_privatekey_t sk;
    BOOST_TEST(!picnic_keygen(parameters, &pk, &sk));

    std::vector<uint8_t> signature;
    signature.resize(max_signature_size);
    size_t signature_len = max_signature_size;

    std::array<uint8_t, 256> message;
    randomize_container(message);

    BOOST_TEST(!picnic_sign(&sk, message.data(), message.size(), signature.data(), &signature_len));
    BOOST_TEST(
        !picnic_verify(&pk, message.data(), message.size(), signature.data(), signature_len));

    std::vector<uint8_t> buf;
    buf.resize(picnic_get_public_key_size(parameters));
    randomize_container(buf);

    const auto pk_written = picnic_write_public_key(&pk, buf.data(), buf.size());
    BOOST_TEST(pk_written > 0);
    BOOST_TEST(pk_written <= buf.size());
    buf.resize(pk_written);

    picnic_publickey_t pk2;
    BOOST_TEST(!picnic_read_public_key(&pk2, buf.data(), buf.size()));
    BOOST_TEST(
        !picnic_verify(&pk2, message.data(), message.size(), signature.data(), signature_len));

    buf.resize(picnic_get_private_key_size(parameters));
    randomize_container(buf);
    const auto sk_written = picnic_write_private_key(&sk, buf.data(), buf.size());
    BOOST_TEST(sk_written > 0);
    BOOST_TEST(sk_written <= buf.size());
    buf.resize(sk_written);

    picnic_privatekey_t sk2;
    BOOST_TEST(!picnic_read_private_key(&sk2, buf.data(), buf.size()));
    BOOST_TEST(!picnic_validate_keypair(&sk2, &pk2));
  }
}

BOOST_DATA_TEST_CASE(modified_signature, all_supported_parameters(), parameters) {
  BOOST_TEST_CONTEXT("Parameter set: " << picnic_get_param_name(parameters)) {
    const size_t max_signature_size = picnic_signature_size(parameters);
    BOOST_TEST(max_signature_size);

    picnic_publickey_t pk;
    picnic_privatekey_t sk;
    BOOST_TEST(!picnic_keygen(parameters, &pk, &sk));

    std::vector<uint8_t> signature, message;
    signature.resize(max_signature_size);

    const size_t message_size = message_inc;
    message.resize(message_size);
    // fill message with some data
    randomize_container(message);

    size_t signature_len = max_signature_size;
    BOOST_TEST(!picnic_sign(&sk, message.data(), message_size, signature.data(), &signature_len));
    signature.resize(signature_len);
    BOOST_TEST(
        !picnic_verify(&pk, message.data(), message_size, signature.data(), signature.size()));

    std::uniform_int_distribution<uint8_t> dist_diff{1};
    std::uniform_int_distribution<size_t> dist_pos{1, signature.size() / 8};
    std::random_device rnd;
    std::default_random_engine eng{rnd()};

    for (size_t l = dist_pos(eng), rounds = rep; rounds; l += dist_pos(eng), --rounds) {
      std::vector<uint8_t> signature2{signature};
      signature2[l % signature2.size()] += dist_diff(eng);

      // must fail
      BOOST_TEST(
          picnic_verify(&pk, message.data(), message_size, signature2.data(), signature2.size()));
    }
  }
}
