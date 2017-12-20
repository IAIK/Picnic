/*
 *  This file is part of the optimized implementation of the Picnic signature scheme.
 *  See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */

#ifndef PICNIC_COMPAT_ENDIAN_H
#define PICNIC_COMPAT_ENDIAN_H

#include <stdint.h>

#if defined(__GCC__) || defined(__clang__)
#define bswap16(x) __builtin_bswap16(x)
#define bswap32(x) __builtin_bswap32(x)
#define bswap64(x) __builtin_bswap64(x)
#elif defined(_MSC_VER)
#include <stdlib.h>

#define bswap16(x) _byteswap_ushort(x)
#define bswap32(x) _byteswap_ulong(x)
#define bswap64(x) _byteswap_uint64(x)
#else
static inline uint16_t bswap16(uint16_t x) {
  return ((x & 0xff00) >> 8) | ((x & 0x00ff) << 8);
}

static inline uint32_t bswap32(uint32_t x) {
  return ((x & 0xff000000) >> 24) | ((x & 0x00ff0000) >> 8) | ((x & 0x0000ff00) << 8) |
         ((x & 0x000000ff) << 24);
}

static inline uint64_t bswap64(uint64_t x) {
  return ((x & UINT64_C(0xff00000000000000)) >> 56) | ((x & UINT64_C(0x00ff000000000000)) >> 40) |
         ((x & UINT64_C(0x0000ff0000000000)) >> 24) | ((x & UINT64_C(0x000000ff00000000)) >> 8) |
         ((x & UINT64_C(0x00000000ff000000)) << 8) | ((x & UINT64_C(0x0000000000ff0000)) << 24) |
         ((x & UINT64_C(0x000000000000ff00)) << 40) | ((x & UINT64_C(0x00000000000000ff)) << 56);
}
#endif

/* Linux / GLIBC */
#if defined(__linux__) || defined(__GLIBC__)
#include <endian.h>
#define HAVE_HOSTSWAP
#define _BYTE_ORDER __BYTE_ORDER
#define _LITTLE_ENDIAN __LITTLE_ENDIAN
#define _BIG_ENDIAN __BIG_ENDIAN
#endif

/* Windows */
#if defined(_WIN16) || defined(_WIN32) || defined(_WIN64)
#if defined(__MINGW32__) || defined(__MINGW64__)
#include <sys/param.h>

#define _LITTLE_ENDIAN LITTLE_ENDIAN
#define _BIG_ENDIAN BIG_ENDIAN
#define _BYTE_ORDER BYTE_ORDER
#else

#define _LITTLE_ENDIAN 1234
#define _BIG_ENDIAN 4321
/* X-Box 360 is big-endian, but we simply ignore that. */
#define _BYTE_ORDER _LITTLE_ENDIAN
#endif
#endif

/* OS X */
#if defined(__APPLE__)
#include <machine/endian.h>

#define _BYTE_ORDER BYTE_ORDER
#define _LITTLE_ENDIAN LITTLE_ENDIAN
#define _BIG_ENDIAN BIG_ENDIAN
#endif

#if !defined(HAVE_HOSTSWAP)
#if _BYTE_ORDER == _LITTLE_ENDIAN
#define htobe16(x) bswap16((x))
#define htole16(x) ((uint16_t)(x))
#define be16toh(x) bswap16((x))
#define le16toh(x) ((uint16_t)(x))

#define htobe32(x) bswap32((x))
#define htole32(x) ((uint32_t)(x))
#define be32toh(x) bswap32((x))
#define le32toh(x) ((uint32_t)(x))

#define htobe64(x) bswap64((x))
#define htole64(x) ((uint64_t)(x))
#define be64toh(x) bswap64((x))
#define le64toh(x) ((uint64_t)(x))
#elif _BYTE_ORDER == _BIG_ENDIAN
#define htobe16(x) ((uint16_t)(x))
#define htole16(x) bswap16((x))
#define be16toh(x) ((uint16_t)(x))
#define le16toh(x) bswap16((x))

#define htobe32(x) ((uint32_t)(x))
#define htole32(x) bswap32((x))
#define be32toh(x) ((uint32_t)(x))
#define le32toh(x) bswap32((x))

#define htobe64(x) ((uint64_t)(x))
#define htole64(x) bswap64((x))
#define be64toh(x) ((uint64_t)(x))
#define le64toh(x) bswap64((x))
#endif
#endif

#endif
