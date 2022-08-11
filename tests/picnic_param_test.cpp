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

#include <array>
#include <vector>

#if !defined(PARAM) || !defined(HEADER)
#error "Missing PARAM and HEADER definitions"
#endif

#include HEADER
#include "utils.h"

#define BOOST_TEST_MODULE PICNIC_CONCAT(PARAM, test)
#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

BOOST_AUTO_TEST_CASE(sign_verify) {
  const size_t max_signature_size = PICNIC_CONCAT(PARAM, signature_size)();
  BOOST_TEST(max_signature_size);

  std::array<uint8_t, 32> m;
  randomize_container(m);

  // Create a key pair
  PICNIC_CONCAT(PARAM, privatekey_t) private_key;
  PICNIC_CONCAT(PARAM, publickey_t) public_key;
  BOOST_TEST(!PICNIC_CONCAT(PARAM, keygen)(&public_key, &private_key));

  // Valid key pair
  BOOST_TEST(!PICNIC_CONCAT(PARAM, validate_keypair)(&private_key, &public_key));

  std::vector<uint8_t> sig;
  sig.resize(max_signature_size);
  size_t siglen = max_signature_size;

  // Sign a message
  BOOST_TEST(!PICNIC_CONCAT(PARAM, sign)(&private_key, m.data(), m.size(), sig.data(), &siglen));
  BOOST_TEST(siglen > 0);
  BOOST_TEST(siglen <= max_signature_size);
  // Verify signature
  BOOST_TEST(!PICNIC_CONCAT(PARAM, verify)(&public_key, m.data(), m.size(), sig.data(), siglen));
}

BOOST_AUTO_TEST_CASE(sign_verify_generic) {
  const size_t max_signature_size = PICNIC_CONCAT(PARAM, signature_size)();
  BOOST_TEST(max_signature_size);

  std::array<uint8_t, 32> m;
  randomize_container(m);

  // Create a key pair
  PICNIC_CONCAT(PARAM, privatekey_t) private_key;
  PICNIC_CONCAT(PARAM, publickey_t) public_key;
  BOOST_TEST(!PICNIC_CONCAT(PARAM, keygen)(&public_key, &private_key));

  // Valid key pair
  BOOST_TEST(!PICNIC_CONCAT(PARAM, validate_keypair)(&private_key, &public_key));

  // Serialize key pair
  std::vector<uint8_t> serialized_private_key;
  serialized_private_key.resize(PICNIC_CONCAT(PARAM, get_private_key_size()));
  BOOST_TEST(PICNIC_CONCAT(PARAM, write_private_key)(&private_key, serialized_private_key.data(),
                                                     serialized_private_key.size()) ==
             serialized_private_key.size());

  std::vector<uint8_t> serialized_public_key;
  serialized_public_key.resize(PICNIC_CONCAT(PARAM, get_public_key_size()));
  BOOST_TEST(PICNIC_CONCAT(PARAM, write_public_key)(&public_key, serialized_public_key.data(),
                                                    serialized_public_key.size()) ==
             serialized_public_key.size());

  // Deserialize key pair
  picnic_privatekey_t gen_private_key;
  picnic_publickey_t gen_public_key;
  BOOST_TEST(!picnic_read_private_key(&gen_private_key, serialized_private_key.data(),
                                      serialized_private_key.size()));
  BOOST_TEST(!picnic_read_public_key(&gen_public_key, serialized_public_key.data(),
                                     serialized_public_key.size()));

  std::vector<uint8_t> sig;
  sig.resize(max_signature_size);
  size_t siglen = max_signature_size;

  // Sign a message
  BOOST_TEST(!PICNIC_CONCAT(PARAM, sign)(&private_key, m.data(), m.size(), sig.data(), &siglen));
  BOOST_TEST(siglen > 0);
  BOOST_TEST(siglen <= max_signature_size);
  // Verify signature
  BOOST_TEST(!PICNIC_CONCAT(PARAM, verify)(&public_key, m.data(), m.size(), sig.data(), siglen));
  BOOST_TEST(!picnic_verify(&gen_public_key, m.data(), m.size(), sig.data(), siglen));

  siglen = max_signature_size;
  // Sign another message
  BOOST_TEST(!picnic_sign(&gen_private_key, m.data(), m.size(), sig.data(), &siglen));
  BOOST_TEST(siglen > 0);
  BOOST_TEST(siglen <= max_signature_size);
  // Verify signature
  BOOST_TEST(!PICNIC_CONCAT(PARAM, verify)(&public_key, m.data(), m.size(), sig.data(), siglen));
  BOOST_TEST(!picnic_verify(&gen_public_key, m.data(), m.size(), sig.data(), siglen));
}
