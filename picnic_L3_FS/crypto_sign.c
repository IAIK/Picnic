/*
 *  This file is part of the optimized implementation of the Picnic signature scheme.
 *  See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */

#ifdef SUPERCOP
#include "crypto_sign.h"
#else
#include "api.h"
#endif

#include "endian_compat.h"
#include "picnic.h"

#include <string.h>

// This file provides the implementation of the SUPERCOP / NIST competition API for a specific
// Picnic instance. For this implementation, we read/write the keys directly from the user provider
// memory areas.
//
// As those interfaces require the "signature" to include message, the signature is serialized as
// * 4-byte signature length (little endian)
// * the message message
// * the actual signature

int crypto_sign_keypair(unsigned char* pk, unsigned char* sk) {
  picnic_publickey_t ppk;
  picnic_privatekey_t psk;

  int ret = picnic_keygen(Picnic_L3_FS, &ppk, &psk);
  if (ret) {
    return ret;
  }

  ret = picnic_write_public_key(&ppk, pk, PICNIC_PUBLIC_KEY_SIZE(Picnic_L3_FS));
  if (ret < 0) {
    picnic_clear_private_key(&psk);
    return ret;
  }

  ret = picnic_write_private_key(&psk, sk, PICNIC_PRIVATE_KEY_SIZE(Picnic_L3_FS));
  picnic_clear_private_key(&psk);
  if (ret < 0) {
    return ret;
  }

  return 0;
}

int crypto_sign(unsigned char* sm, unsigned long long* smlen, const unsigned char* m,
                unsigned long long mlen, const unsigned char* sk) {
  // The first byte encodes the parameter set and is public.
  picnic_declassify(&sk[0], sizeof(unsigned char));
  if (sk[0] != Picnic_L3_FS) {
    return -3;
  }

  size_t signature_len = PICNIC_SIGNATURE_SIZE(Picnic_L3_FS);
  uint32_t len         = 0;

  picnic_privatekey_t psk;
  int ret = picnic_read_private_key(&psk, sk, PICNIC_PRIVATE_KEY_SIZE(Picnic_L3_FS));
  if (ret < 0) {
    picnic_clear_private_key(&psk);
    return ret;
  }

  ret = picnic_sign(&psk, m, mlen, sm + sizeof(len) + mlen, &signature_len);
  picnic_clear_private_key(&psk);
  if (ret) {
    return ret;
  }

  len    = htole32(signature_len);
  *smlen = sizeof(len) + mlen + signature_len;
  // Move the message first in case m and sm overlap.
  memmove(sm + sizeof(len), m, mlen);
  memcpy(sm, &len, sizeof(len));

  return 0;
}

int crypto_sign_open(unsigned char* m, unsigned long long* mlen, const unsigned char* sm,
                     unsigned long long smlen, const unsigned char* pk) {
  if (pk[0] != Picnic_L3_FS) {
    return -3;
  }

  uint32_t signature_len = 0;
  // The signature is too short to hold the signature length.
  if (smlen < sizeof(signature_len)) {
    return -2;
  }

  memcpy(&signature_len, sm, sizeof(signature_len));
  signature_len = le32toh(signature_len);
  // The signature is too short to hold the signature
  if (signature_len + sizeof(signature_len) > smlen) {
    return -2;
  }

  picnic_publickey_t ppk;
  int ret = picnic_read_public_key(&ppk, pk, PICNIC_PUBLIC_KEY_SIZE(Picnic_L3_FS));
  if (ret < 0) {
    return ret;
  }

  const size_t message_len = smlen - signature_len - sizeof(signature_len);
  const uint8_t* message   = sm + sizeof(signature_len);
  const uint8_t* sig       = sm + sizeof(signature_len) + message_len;

  ret = picnic_verify(&ppk, message, message_len, sig, signature_len);
  if (ret) {
    return ret;
  }

  // The signature is valid, so copy the message to its destination.
  memmove(m, message, message_len);
  *mlen = message_len;

  return 0;
}

// vim: ft=c
