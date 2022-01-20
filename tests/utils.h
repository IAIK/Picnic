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

namespace {
  constexpr picnic_params_t all_parameters[] = {
      Picnic_L1_FS,   Picnic_L1_UR, Picnic_L1_full, Picnic3_L1,   Picnic_L3_FS,   Picnic_L3_UR,
      Picnic_L3_full, Picnic3_L3,   Picnic_L5_FS,   Picnic_L5_UR, Picnic_L5_full, Picnic3_L5,
  };

  template <typename C>
  void randomize_container(C& container) {
    std::uniform_int_distribution<typename C::value_type> dist;
    std::random_device rnd;
    std::default_random_engine eng(rnd());
    std::generate(container.begin(), container.end(), [&dist, &eng] { return dist(eng); });
  }
} // namespace

extern "C" {
#endif
picnic_params_t argument_to_params(const char* arg);
#if defined(__cplusplus)
}
#endif

#endif
