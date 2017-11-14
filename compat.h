/*
 *  This file is part of the optimized implementation of the Picnic signature scheme.
 *  See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */

#ifndef PICNIC_COMPAT_H
#define PICNIC_COMPAT_H

/* in case cmake checks failed or were not run, define HAVE_* for known good
 * configurations */

#if !defined(HAVE_ALIGNED_ALLOC) &&                                                                \
    (defined(_ISOC11_SOURCE) ||                                                                    \
     (defined(__STDC_VERSION__) && __STDC_VERSION__ >= 201112L && !defined(__MINGW32__) &&         \
      !defined(__MINGW64__) && !defined(__APPLE__)))
#define HAVE_ALIGNED_ALLOC
#endif

#if !defined(HAVE_POSIX_MEMALIGN) && defined(_POSIX_C_SOURCE) && _POSIX_C_SOURCE >= 200112L
#define HAVE_POSIX_MEMALIGN
#endif

#if !defined(HAVE_MEMALIGN) && defined(__linux__)
#define HAVE_MEMALIGN
#endif

#if defined(HAVE_ALIGNED_ALLOC)
#include <stdlib.h>

#define aligned_free free
#else
#include <stddef.h>

/* compat implementation of aligned_alloc from ISO C 2011 */
void* aligned_alloc(size_t alignment, size_t size);
/* some aligned alloc implementations require custom free functions, so we
 * provide one too */
void aligned_free(void* ptr);
#endif

/* backwards compatibility macros for GCC 4.8 and 4.9
 *
 * bs{l,r}i was introduced in GCC 5 and in clang as macros sometime in 2015.
 * */
#ifdef WITH_OPT
#if (!defined(__clang__) && defined(__GNUC__) && __GNUC__ < 5) ||                                  \
    (defined(__clang__) && !defined(_mm_bslli_si128))
#define _mm_bslli_si128(a, imm) _mm_slli_si128((a), (imm))
#define _mm_bsrli_si128(a, imm) _mm_srli_si128((a), (imm))
#endif
#endif

#include "endian_compat.h"

#endif
