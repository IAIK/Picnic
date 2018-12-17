/*
 *  This file is part of the optimized implementation of the Picnic signature scheme.
 *  See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */

#ifndef BENCH_TIMING_H
#define BENCH_TIMING_H

#include <stdint.h>
#include <stdbool.h>

typedef struct timing_context_s timing_context_t;

typedef uint64_t (*timing_read_f)(timing_context_t* ctx);
typedef void (*timing_close_f)(timing_context_t* ctx);

struct timing_context_s {
  timing_read_f read;
  timing_close_f close;

  union {
    int fd;
  } data;
};

bool timing_init(timing_context_t* ctx);

static inline uint64_t timing_read(timing_context_t* ctx) {
  return ctx->read(ctx);
}

static inline void timing_close(timing_context_t* ctx) {
  ctx->close(ctx);
}

#endif
