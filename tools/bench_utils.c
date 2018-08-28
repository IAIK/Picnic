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

#include "bench_utils.h"

#include <errno.h>
#if !defined(_MSC_VER)
#include <getopt.h>
#endif
#include <limits.h>
#include <stdlib.h>

static bool parse_long(long* value, const char* arg) {
  errno        = 0;
  const long v = strtol(arg, NULL, 10);

  if ((errno == ERANGE && (v == LONG_MAX || v == LONG_MIN)) || (errno != 0 && v == 0)) {
    return false;
  }
  *value = v;

  return true;
}

static bool parse_uint32_t(uint32_t* value, const char* arg) {
  long tmp = 0;
  if (!parse_long(&tmp, arg)) {
    return false;
  }

  if (tmp < 0 || (unsigned long)tmp > UINT32_MAX) {
    return false;
  }

  *value = tmp;
  return true;
}

static void print_usage(const char* arg0) {
#if defined(_MSC_VER)
  printf("usage: %s iterations instance\n", arg0);
#else
  printf("usage: %s [-i iterations] instance\n", arg0);
#endif
}

bool parse_args(bench_options_t* options, int argc, char** argv) {
  if (argc <= 1) {
    print_usage(argv[0]);
    return false;
  }

  options->params = PARAMETER_SET_INVALID;
  options->iter   = 10;

#if !defined(_MSC_VER)
  static const struct option long_options[] = {
    {"iter", required_argument, NULL, 'i'},
#if defined(WITH_CUSTOM_INSTANCES)
    {"lowmc", required_argument, NULL, 'l'},
#endif
    {0, 0, 0, 0}
  };

  int c            = -1;
  int option_index = 0;

  while ((c = getopt_long(argc, argv, "i:l:", long_options, &option_index)) != -1) {
    switch (c) {
    case 'i':
      if (!parse_uint32_t(&options->iter, optarg)) {
        printf("Failed to parse argument as positive base-10 number!\n");
        return false;
      }
      break;

#if defined(WITH_CUSTOM_INSTANCES)
    case 'l':
      options->lowmc_file = optarg;
      break;
#endif

    case '?':
    default:
#if defined(WITH_CUSTOM_INSTANCES)
      printf("usage: %s [-i iter] [-l lowmc] param\n", argv[0]);
#else
      printf("usage: %s [-i iter] param\n", argv[0]);
#endif
      return false;
    }
  }

  if (optind == argc - 1) {
    uint32_t p = -1;
    if (!parse_uint32_t(&p, argv[optind])) {
      printf("Failed to parse argument as positive base-10 number!\n");
      return false;
    }

    if (p <= PARAMETER_SET_INVALID || p >= PARAMETER_SET_MAX_INDEX) {
      printf("Invalid parameter set selected!\n");
      return false;
    }
    options->params = p;
  } else {
    print_usage(argv[0]);
    return false;
  }
#else
  if (argc != 3) {
    print_usage(argv[0]);
    return false;
  }

  uint32_t p = -1;
  if (!parse_uint32_t(&options->iter, argv[1]) || !parse_uint32_t(&p, argv[2])) {
    printf("Failed to parse argument as positive base-10 number!\n");
    return false;
  }

  if (p <= PARAMETER_SET_INVALID || p >= PARAMETER_SET_MAX_INDEX) {
    printf("Invalid parameter set selected!\n");
    return false;
  }
  options->params = p;
#endif

  return true;
}
