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

#include "bench_timing.h"

#include <time.h>
#include <limits.h>
#include <string.h>

#if defined(__linux__) && defined(__aarch64__)
#include <setjmp.h>
#include <signal.h>

/* Based on code from https://github.com/IAIK/armageddon/tree/master/libflush
 *
 * Copyright (c) 2015-2016 Moritz Lipp
 *
 * This software is provided 'as-is', without any express or implied
 * warranty. In no event will the authors be held liable for any damages
 * arising from the use of this software.
 *
 * Permission is granted to anyone to use this software for any purpose,
 * including commercial applications, and to alter it and redistribute it
 * freely, subject to the following restrictions:
 *
 *   1. The origin of this software must not be misrepresented; you must not
 *   claim that you wrote the original software. If you use this software
 *   in a product, an acknowledgment in the product documentation would be
 *   appreciated but is not required.
 *
 *   2. Altered source versions must be plainly marked as such, and must not be
 *   misrepresented as being the original software.
 *
 *   3. This notice may not be removed or altered from any source
 *   distribution. */

#define ARMV8_PMCR_E (1 << 0) /* Enable all counters */
#define ARMV8_PMCR_P (1 << 1) /* Reset all counters */
#define ARMV8_PMCR_C (1 << 2) /* Cycle counter reset */

#define ARMV8_PMUSERENR_EN (1 << 0) /* EL0 access enable */
#define ARMV8_PMUSERENR_CR (1 << 2) /* Cycle counter read enable */
#define ARMV8_PMUSERENR_ER (1 << 3) /* Event counter read enable */

#define ARMV8_PMCNTENSET_EL0_EN (1 << 31) /* Performance Monitors Count Enable Set register */

static void armv8_close(timing_context_t* ctx) {
  (void)ctx;
  uint32_t value = 0;
  uint32_t mask  = 0;

  /* Disable Performance Counter */
  asm volatile("MRS %0, PMCR_EL0" : "=r"(value));
  mask = 0;
  mask |= ARMV8_PMCR_E; /* Enable */
  mask |= ARMV8_PMCR_C; /* Cycle counter reset */
  mask |= ARMV8_PMCR_P; /* Reset all counters */
  asm volatile("MSR PMCR_EL0, %0" : : "r"(value & ~mask));

  /* Disable cycle counter register */
  asm volatile("MRS %0, PMCNTENSET_EL0" : "=r"(value));
  mask = 0;
  mask |= ARMV8_PMCNTENSET_EL0_EN;
  asm volatile("MSR PMCNTENSET_EL0, %0" : : "r"(value & ~mask));
}

static uint64_t armv8_read(timing_context_t* ctx) {
  (void)ctx;
  uint64_t result = 0;
  asm volatile("MRS %0, PMCCNTR_EL0" : "=r"(result));
  return result;
}

static sigjmp_buf jmpbuf;
static volatile sig_atomic_t armv8_sigill = 0;

static void armv8_sigill_handler(int sig) {
  (void)sig;
  armv8_sigill = 1;
  // Return to sigsetjump
  siglongjmp(jmpbuf, 1);
}

static bool armv8_init(timing_context_t* ctx) {
  if (armv8_sigill) {
    return false;
  }

  struct sigaction act, oldact;
  memset(&act, 0, sizeof(act));
  act.sa_handler = &armv8_sigill_handler;
  if (sigaction(SIGILL, &act, &oldact) < 0) {
    return false;
  }

  if (sigsetjmp(jmpbuf, 1)) {
    // Returned from armv8_sigill_handler
    sigaction(SIGILL, &oldact, NULL);
    return false;
  }

  uint32_t value = 0;

  /* Enable Performance Counter */
  asm volatile("MRS %0, PMCR_EL0" : "=r"(value));
  value |= ARMV8_PMCR_E; /* Enable */
  value |= ARMV8_PMCR_C; /* Cycle counter reset */
  value |= ARMV8_PMCR_P; /* Reset all counters */
  asm volatile("MSR PMCR_EL0, %0" : : "r"(value));

  /* Enable cycle counter register */
  asm volatile("MRS %0, PMCNTENSET_EL0" : "=r"(value));
  value |= ARMV8_PMCNTENSET_EL0_EN;
  asm volatile("MSR PMCNTENSET_EL0, %0" : : "r"(value));

  // Restore old signal handler
  sigaction(SIGILL, &oldact, NULL);

  ctx->read  = armv8_read;
  ctx->close = armv8_close;

  return true;
}
#endif

#if defined(__linux__)
#include <linux/perf_event.h>
#include <linux/version.h>
#include <sys/syscall.h>
#include <unistd.h>

static void perf_close(timing_context_t* ctx) {
  if (ctx->data.fd != -1) {
    close(ctx->data.fd);
    ctx->data.fd = -1;
  }
}

static uint64_t perf_read(timing_context_t* ctx) {
  uint64_t tmp_time;
  if (read(ctx->data.fd, &tmp_time, sizeof(tmp_time)) != sizeof(tmp_time)) {
    return UINT64_MAX;
  }

  return tmp_time;
}

static int perf_event_open(struct perf_event_attr* event, pid_t pid, int cpu, int gfd,
                           unsigned long flags) {
  const long fd = syscall(__NR_perf_event_open, event, pid, cpu, gfd, flags);
  if (fd > INT_MAX) {
    /* too large to handle, but should never happen */
    return -1;
  }

  return fd;
}

static bool perf_init(timing_context_t* ctx) {
  struct perf_event_attr pea;
  memset(&pea, 0, sizeof(pea));

  pea.size           = sizeof(pea);
  pea.type           = PERF_TYPE_HARDWARE;
  pea.config         = PERF_COUNT_HW_CPU_CYCLES;
  pea.disabled       = 0;
  pea.exclude_kernel = 1;
  pea.exclude_hv     = 1;
#if LINUX_VERSION_CODE >= KERNEL_VERSION(3, 7, 0)
  pea.exclude_callchain_kernel = 1;
  pea.exclude_callchain_user   = 1;
#endif

  const int fd = perf_event_open(&pea, 0, -1, -1, 0);
  if (fd == -1) {
    return false;
  }

  ctx->read    = perf_read;
  ctx->close   = perf_close;
  ctx->data.fd = fd;
  return true;
}
#endif

static void clock_close(timing_context_t* ctx) {
  (void)ctx;
}

static uint64_t clock_read(timing_context_t* ctx) {
  (void)ctx;
  return clock() * (1000000 / CLOCKS_PER_SEC);
}

static bool clock_init(timing_context_t* ctx) {
  ctx->read  = clock_read;
  ctx->close = clock_close;
  return true;
}

bool timing_init(timing_context_t* ctx) {
#if defined(__linux__) && defined(__aarch64__)
  if (armv8_init(ctx)) {
    return true;
  }
#endif
#if defined(__linux__)
  if (perf_init(ctx)) {
    return true;
  }
#endif
  return clock_init(ctx);
}
