/*
 *  This file is part of the optimized implementation of the Picnic signature scheme.
 *  See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */
#ifndef PICNIC_TEST_UTILS_H
#define PICNIC_TEST_UTILS_H

#include <picnic.h>

#if defined(__cplusplus)
#include <algorithm>
#include <random>
#include <vector>
#include <limits>

namespace {
  template <typename C>
  void randomize_container(C& container) {
    std::uniform_int_distribution<unsigned int> dist{
        0, std::numeric_limits<typename C::value_type>::max()};
    std::random_device rnd;
    std::default_random_engine eng(rnd());
    std::generate(container.begin(), container.end(), [&dist, &eng] { return dist(eng); });
  }
} // namespace

std::vector<picnic_params_t> all_supported_parameters();

extern "C" {
#endif
picnic_params_t argument_to_params(const char* arg);
#if defined(__cplusplus)
}
#endif

#endif
