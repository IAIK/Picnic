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

#include <vector>

#include <picnic.h>
#include "utils.h"

#define BOOST_TEST_MODULE picnic_test
#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

BOOST_DATA_TEST_CASE(sign_verify, all_supported_parameters(), param) {
  BOOST_TEST_CONTEXT("Parameter set: " << picnic_get_param_name(param)) {
    static constexpr uint8_t m[] = "test message";

    const size_t max_signature_size = picnic_signature_size(param);
    BOOST_TEST(max_signature_size);

    // Create a key pair
    picnic_privatekey_t private_key;
    picnic_publickey_t public_key;
    BOOST_TEST(!picnic_keygen(param, &public_key, &private_key));

    // Valid key pair
    BOOST_TEST(picnic_get_public_key_param(&public_key) == param);
    BOOST_TEST(picnic_get_private_key_param(&private_key) == param);
    BOOST_TEST(!picnic_validate_keypair(&private_key, &public_key));

    std::vector<uint8_t> sig;
    sig.resize(max_signature_size);
    size_t siglen = max_signature_size;

    // Sign a message
    BOOST_TEST(!picnic_sign(&private_key, m, sizeof(m), sig.data(), &siglen));
    BOOST_TEST(siglen > 0);
    BOOST_TEST(siglen <= max_signature_size);
    // Verify signature
    BOOST_TEST(!picnic_verify(&public_key, m, sizeof(m), sig.data(), siglen));
  }
}
