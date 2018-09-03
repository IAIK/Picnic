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

#include "io.h"

#include <string.h>
#include "compat.h"

void mzd_to_char_array(uint8_t* dst, const mzd_local_t* data, unsigned len) {
  const size_t word_count = len / sizeof(uint64_t);
  const uint64_t* rows    = &CONST_FIRST_ROW(data)[word_count - 1];

  for (size_t i = word_count; i; --i, --rows, dst += sizeof(uint64_t)) {
    const uint64_t tmp = htobe64(*rows);
    memcpy(dst, &tmp, sizeof(tmp));
  }
}

void mzd_from_char_array(mzd_local_t* result, const uint8_t* data, unsigned len) {
  const size_t word_count = len / sizeof(uint64_t);
  uint64_t* rows          = &FIRST_ROW(result)[word_count - 1];

  for (size_t i = word_count; i; --i, --rows, data += sizeof(uint64_t)) {
    uint64_t tmp;
    memcpy(&tmp, data, sizeof(tmp));
    *rows = be64toh(tmp);
  }
}

void print_hex(FILE* out, const uint8_t* data, size_t len) {
  for (size_t i = len; i; --i, ++data) {
    fprintf(out, "%02X", *data);
  }
}
