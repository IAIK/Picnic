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
#include "bench_utils.h"

#include <random>
#include <algorithm>
#include <array>
#include <vector>
#include <iostream>
#include <chrono>

using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::microseconds;

namespace {
  struct timing_and_size_t {
    microseconds keygen, sign, verify;
    std::size_t size, max_size;
  };

  void print_timings(const std::vector<timing_and_size_t>& timings) {
    for (const auto& timing : timings) {
      std::cout << timing.keygen.count() << ',' << timing.sign.count() << ','
                << timing.verify.count() << ',' << timing.size << ',' << timing.max_size
                << std::endl;
    }
  }

  void bench_sign_and_verify(const bench_options_t& options) {
    const size_t max_signature_size = picnic_signature_size(options.params);
    if (!max_signature_size) {
      std::cout << "Failed to create Picnic instance." << std::endl;
      return;
    }

    std::vector<timing_and_size_t> timings;
    timings.reserve(options.iter);
    std::vector<uint8_t> sig;
    sig.resize(max_signature_size);

    std::array<uint8_t, 32> m;
    {
      std::uniform_int_distribution<uint8_t> dist;
      std::random_device rnd;
      std::default_random_engine eng(rnd());
      std::generate(m.begin(), m.end(), [&dist, &eng] { return dist(eng); });
    }

    for (unsigned int i = 0; i != options.iter; ++i) {
      timing_and_size_t timing;
      timing.max_size = max_signature_size;

      auto start_time = high_resolution_clock::now();
      picnic_privatekey_t private_key;
      picnic_publickey_t public_key;

      if (picnic_keygen(options.params, &public_key, &private_key)) {
        std::cout << "picnic_keygen: failed." << std::endl;
        break;
      }
      timing.keygen = duration_cast<microseconds>(high_resolution_clock::now() - start_time);

      start_time    = high_resolution_clock::now();
      size_t siglen = max_signature_size;
      if (picnic_sign(&private_key, m.data(), m.size(), sig.data(), &siglen)) {
        std::cout << "picnic_sign: failed" << std::endl;
        break;
      }
      timing.sign = duration_cast<microseconds>(high_resolution_clock::now() - start_time);
      timing.size = siglen;

      start_time = high_resolution_clock::now();
      if (picnic_verify(&public_key, m.data(), m.size(), sig.data(), siglen)) {
        std::cout << "picnic_verify: failed" << std::endl;
        break;
      }
      timing.verify = duration_cast<microseconds>(high_resolution_clock::now() - start_time);
      timings.emplace_back(timing);
    }

    print_timings(timings);
  }
} // namespace

int main(int argc, char** argv) {
  bench_options_t opts = {PARAMETER_SET_INVALID, 0};
  int ret              = parse_args(&opts, argc, argv) ? 0 : -1;

  if (!ret) {
    bench_sign_and_verify(opts);
  }

  return ret;
}
