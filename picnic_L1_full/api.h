/*
 *  This file is part of the optimized implementation of the Picnic signature scheme.
 *  See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */

#ifndef PICNIC_L1_FULL_API_H
#define PICNIC_L1_FULL_API_H

#define CRYPTO_SECRETKEYBYTES (1 + 2 * 17 + 17)
#define CRYPTO_PUBLICKEYBYTES (1 + 2 * 17)
#define CRYPTO_BYTES (4 + 32061)
#define CRYPTO_ALGNAME "picnicl1full"
#define CRYPTO_DETERMINISTIC 1

int crypto_sign_keypair(unsigned char* pk, unsigned char* sk);
int crypto_sign(unsigned char* sm, unsigned long long* smlen, const unsigned char* m,
                unsigned long long mlen, const unsigned char* sk);
int crypto_sign_open(unsigned char* m, unsigned long long* mlen, const unsigned char* sm,
                     unsigned long long smlen, const unsigned char* pk);

#endif
