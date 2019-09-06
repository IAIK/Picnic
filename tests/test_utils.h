#ifndef PICNIC_TEST_UTILS_H
#define PICNIC_TEST_UTILS_H

#include "picnic.h"

#include <errno.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

static inline picnic_params_t argument_to_params(const char* arg, bool support_m1) {
  for (unsigned int param = Picnic_L1_FS; param < PARAMETER_SET_MAX_INDEX; ++param) {
    if (!strcasecmp(arg, picnic_get_param_name(param))) {
      return param;
    }
  }

  const long idx = strtol(arg, NULL, 10);
  if ((errno == ERANGE && (idx == LONG_MAX || idx == LONG_MIN)) || (errno != 0 && idx == 0) ||
      idx < 1 || (size_t)idx >= PARAMETER_SET_MAX_INDEX) {
    return PARAMETER_SET_INVALID;
  }
  if (!support_m1 && (size_t)idx > Picnic2_L5_FS) {
    return PARAMETER_SET_INVALID;
  }

  return idx;
}

#endif
