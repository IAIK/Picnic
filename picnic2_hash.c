/*! @file hash.c
 *  @brief Wraps the SHA-3 implementation.
 *
 *  This file is part of the reference implementation of the Picnic signature scheme.
 *  See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */

#include "picnic2_hash.h"
#include <stdio.h>
#include <assert.h>
#if defined(__WINDOWS__)
    #include <Windows.h>
    #include <bcrypt.h>
#else
    #include <endian.h>
#endif


void HashUpdate(HashInstance* ctx, const uint8_t* data, size_t byteLen)
{
    HashReturn ret = Keccak_HashUpdate(ctx, data, byteLen * 8);

    if (ret != SUCCESS) {
        fprintf(stderr, "%s: Keccak_HashUpdate failed (returned %d)\n", __func__, ret);
        assert(!"Keccak_HashUpdate failed");
    }
}

void HashInit(HashInstance* ctx, const picnic_instance_t* params, uint8_t hashPrefix)
{
    if (params->lowmc->n == 128) {         /* L1 */
        Keccak_HashInitialize_SHAKE128(ctx);
    }
    else {                                      /* L3, L5 */
        Keccak_HashInitialize_SHAKE256(ctx);
    }

    if (hashPrefix != HASH_PREFIX_NONE) {
        HashUpdate(ctx, &hashPrefix, 1);
    }
}

void HashFinal(HashInstance* ctx)
{
    HashReturn ret = Keccak_HashFinal(ctx, NULL);

    if (ret != SUCCESS) {
        fprintf(stderr, "%s: Keccak_HashFinal failed (returned %d)\n", __func__, ret);
    }
}


void HashSqueeze(HashInstance* ctx, uint8_t* digest, size_t byteLen)
{
    HashReturn ret = Keccak_HashSqueeze(ctx, digest, byteLen * 8);

    if (ret != SUCCESS) {
        fprintf(stderr, "%s: Keccak_HashSqueeze failed (returned %d)\n", __func__, ret);
    }
}

uint16_t toLittleEndian(uint16_t x)
{
#if defined(__WINDOWS__)
    #if BYTE_ORDER == LITTLE_ENDIAN
        return x;
    #else
        return __builtin_bswap16(x);
    #endif
#else
    return htole16(x);
#endif
}

uint16_t fromLittleEndian(uint16_t x)
{
#if defined(__WINDOWS__)
    #if BYTE_ORDER == LITTLE_ENDIAN
        return x;
    #else
        return __builtin_bswap16(x);
    #endif
#else
    return le16toh(x);
#endif
}

void HashUpdateIntLE(HashInstance* ctx, uint16_t x)
{
    uint16_t outputBytesLE = toLittleEndian(x);

    HashUpdate(ctx, (uint8_t*)&outputBytesLE, sizeof(uint16_t));
}

