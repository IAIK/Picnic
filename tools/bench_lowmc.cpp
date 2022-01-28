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

// FIXME: this should not be needed
#define restrict

#include "bench_utils.h"
#include "../io.h"
#include "../lowmc.h"
#include "../picnic_instances.h"

#include <algorithm>
#include <chrono>
#include <iostream>
#include <random>
#include <vector>

using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::microseconds;

namespace {
  void print_timings(const std::vector<microseconds>& timings) {
    for (const auto& ms : timings) {
      std::cout << ms.count() << std::endl;
    }
  }

  void bench_lowmc(const bench_options_t& options) {
    const picnic_instance_t* pp = picnic_instance_get(options.params);
    if (!pp) {
      std::cout << "Failed to create Picnic instance." << std::endl;
      return;
    }

    std::vector<microseconds> timings;
    timings.reserve(options.iter);

    const lowmc_parameters_t* lowmc         = &pp->lowmc;
    const lowmc_implementation_f lowmc_impl = lowmc_get_implementation(lowmc);

    mzd_local_t* sk = mzd_local_init(1, lowmc->n);
    mzd_local_t* pt = mzd_local_init(1, lowmc->n);
    mzd_local_t* ct = mzd_local_init(1, lowmc->n);

    const size_t input_size = (lowmc->n + 7) / 8;
    std::vector<uint8_t> rand_buffer;
    rand_buffer.resize(2 * input_size);
    {
      std::uniform_int_distribution<unsigned int> dist{0, 255};
      std::random_device rnd;
      std::default_random_engine eng(rnd());
      std::generate(rand_buffer.begin(), rand_buffer.end(), [&dist, &eng] { return dist(eng); });
    }

    mzd_from_char_array(sk, rand_buffer.data(), input_size);
    mzd_from_char_array(pt, rand_buffer.data() + input_size, input_size);

    for (unsigned int i = 0; i != options.iter; ++i) {
      auto start_time = high_resolution_clock::now();
      lowmc_impl(sk, pt, ct);
      timings.emplace_back(duration_cast<microseconds>(high_resolution_clock::now() - start_time));

      std::swap(pt, ct);
    }

    mzd_local_free(ct);
    mzd_local_free(pt);
    mzd_local_free(sk);

    print_timings(timings);
  }
} // namespace

int main(int argc, char** argv) {
  bench_options_t opts = {PARAMETER_SET_INVALID, 0};
  int ret              = parse_args(&opts, argc, argv) ? 0 : -1;

  if (!ret) {
    bench_lowmc(opts);
  }

  return ret;
}
