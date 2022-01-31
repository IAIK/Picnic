/*
 *  This file is part of the optimized implementation of the Picnic signature scheme.
 *  See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */

#ifndef BENCH_UTILS_H
#define BENCH_UTILS_H

#include <boost/program_options.hpp>
#include <iostream>
#include <tuple>
#include <string>

#include "utils.h"

namespace {
  std::istream& operator>>(std::istream& in, picnic_params_t& param) {
    std::string str;
    in >> str;
    param = argument_to_params(str.c_str());
    if (param == PARAMETER_SET_INVALID) {
      in.setstate(std::ios::failbit);
    }
    return in;
  }

  std::tuple<picnic_params_t, unsigned int> parse_args(int argc, char** argv) {
    using namespace boost::program_options;

    options_description options{"Options"};
    options.add_options()("help", "produce help message");
    options.add_options()("iter,i", value<unsigned int>()->default_value(100),
                          "set number of iterations");
    options.add_options()("params,p", value<picnic_params_t>()->required(), "set parameter set");

    variables_map vm;
    try {
      store(parse_command_line(argc, argv, options), vm);
      notify(vm);

      if (vm.count("help")) {
        std::cout << options << std::endl;
        return std::make_tuple(PARAMETER_SET_INVALID, 0);
      }

      return std::make_tuple(vm["params"].as<picnic_params_t>(), vm["iter"].as<unsigned int>());
    } catch (const boost::exception& e) {
      std::cout << options << std::endl;
      return std::make_tuple(PARAMETER_SET_INVALID, 0);
    }
  }
} // namespace

#endif
