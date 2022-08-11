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

#include "picnic_l1_full.h"

#include <stdlib.h>
#include <string.h>

#include "compat.h"
#include "io.h"
#include "lowmc.h"
#include "lowmc_129_129_4.h"
#include "picnic_impl.h"
#include "picnic_instances.h"
#include "randomness.h"

#define LOWMC_BLOCK_BITS LOWMC_129_129_4_N
#define LOWMC_BLOCK_SZ LOWMC_BLOCK_SIZE_Picnic_L1_full
#define PRIVATE_KEY_SIZE PICNIC_PRIVATE_KEY_SIZE_Picnic_L1_full
#define PUBLIC_KEY_SIZE PICNIC_PUBLIC_KEY_SIZE_Picnic_L1_full
#define SIGNATURE_SIZE PICNIC_SIGNATURE_SIZE_Picnic_L1_full
#define PARAM Picnic_L1_full

// Public and private keys are serialized as follows:
// - public key: instance || C || p
// - secret key: instance || sk || C || p

#define SK_SK(sk) &(sk)->data[0]
#define SK_C(sk) &(sk)->data[LOWMC_BLOCK_SZ]
#define SK_PT(sk) &(sk)->data[2 * LOWMC_BLOCK_SZ]

#define PK_C(pk) &(pk)->data[0]
#define PK_PT(pk) &(pk)->data[LOWMC_BLOCK_SZ]

size_t PICNIC_CALLING_CONVENTION picnic_l1_full_signature_size(void) {
  return SIGNATURE_SIZE;
}

size_t PICNIC_CALLING_CONVENTION picnic_l1_full_get_private_key_size(void) {
  return PRIVATE_KEY_SIZE;
}

size_t PICNIC_CALLING_CONVENTION picnic_l1_full_get_public_key_size(void) {
  return PUBLIC_KEY_SIZE;
}

const char* PICNIC_CALLING_CONVENTION picnic_l1_full_get_param_name(void) {
  return "picnic_l1_full";
}

int PICNIC_CALLING_CONVENTION picnic_l1_full_keygen(picnic_l1_full_publickey_t* pk,
                                                  picnic_l1_full_privatekey_t* sk) {
  if (!pk || !sk) {
    return -1;
  }

  uint8_t* sk_sk = SK_SK(sk);
  uint8_t* sk_pt = SK_PT(sk);
  uint8_t* sk_c  = SK_C(sk);

  // generate private key
  // random secret key
  if (rand_bits(sk_sk, LOWMC_BLOCK_BITS)) {
    return -1;
  }
  // random plain text
  if (rand_bits(sk_pt, LOWMC_BLOCK_BITS)) {
    return -1;
  }
  // encrypt plaintext under secret key
  if (picnic_l1_full_sk_to_pk(sk, pk)) {
    return -1;
  }
  // copy ciphertext to secret key
  memcpy(sk_c, PK_C(pk), LOWMC_BLOCK_SZ);
  return 0;
}

int PICNIC_CALLING_CONVENTION picnic_l1_full_sk_to_pk(const picnic_l1_full_privatekey_t* sk,
                                                    picnic_l1_full_publickey_t* pk) {
  if (!sk || !pk) {
    return -1;
  }

  const picnic_instance_t* instance = picnic_instance_get(PARAM);

  const uint8_t* sk_sk = SK_SK(sk);
  uint8_t* pk_c        = PK_C(pk);
  uint8_t* pk_pt       = PK_PT(pk);
  const uint8_t* sk_pt = SK_PT(sk);

  mzd_local_t plaintext[(MAX_LOWMC_BLOCK_SIZE_BITS + 255) / 256];
  mzd_local_t privkey[(MAX_LOWMC_BLOCK_SIZE_BITS + 255) / 256];
  mzd_local_t ciphertext[(MAX_LOWMC_BLOCK_SIZE_BITS + 255) / 256];

  mzd_from_char_array(plaintext, sk_pt, LOWMC_BLOCK_SZ);
  mzd_from_char_array(privkey, sk_sk, LOWMC_BLOCK_SZ);

  // compute public key
  lowmc_compute(&instance->lowmc, privkey, plaintext, ciphertext);

  memcpy(pk_pt, sk_pt, LOWMC_BLOCK_SZ);
  mzd_to_char_array(pk_c, ciphertext, LOWMC_BLOCK_SZ);

  return 0;
}

int PICNIC_CALLING_CONVENTION picnic_l1_full_validate_keypair(const picnic_l1_full_privatekey_t* sk,
                                                            const picnic_l1_full_publickey_t* pk) {
  if (!sk || !pk) {
    return -1;
  }

  const picnic_instance_t* instance = picnic_instance_get(PARAM);

  const uint8_t* sk_sk = SK_SK(sk);
  const uint8_t* sk_pt = SK_PT(sk);
  const uint8_t* sk_c  = SK_C(sk);
  const uint8_t* pk_pt = PK_PT(pk);
  const uint8_t* pk_c  = PK_C(pk);

  // check param and plaintext
  if (memcmp(sk_pt, pk_pt, LOWMC_BLOCK_SZ) != 0 || memcmp(sk_c, pk_c, LOWMC_BLOCK_SZ) != 0) {
    return -1;
  }

  mzd_local_t plaintext[(MAX_LOWMC_BLOCK_SIZE_BITS + 255) / 256];
  mzd_local_t privkey[(MAX_LOWMC_BLOCK_SIZE_BITS + 255) / 256];
  mzd_local_t ciphertext[(MAX_LOWMC_BLOCK_SIZE_BITS + 255) / 256];

  mzd_from_char_array(plaintext, sk_pt, LOWMC_BLOCK_SZ);
  mzd_from_char_array(privkey, sk_sk, LOWMC_BLOCK_SZ);

  // compute public key
  lowmc_compute(&instance->lowmc, privkey, plaintext, ciphertext);

  uint8_t buffer[MAX_LOWMC_BLOCK_SIZE];
  mzd_to_char_array(buffer, ciphertext, LOWMC_BLOCK_SZ);

  return memcmp(buffer, pk_c, LOWMC_BLOCK_SZ);
}

int PICNIC_CALLING_CONVENTION picnic_l1_full_sign(const picnic_l1_full_privatekey_t* sk,
                                                const uint8_t* message, size_t message_len,
                                                uint8_t* signature, size_t* signature_len) {
  if (!sk || !signature || !signature_len) {
    return -1;
  }

  const picnic_instance_t* instance = picnic_instance_get(PARAM);
  if (!instance) {
    return -1;
  }

  const uint8_t* sk_sk = SK_SK(sk);
  const uint8_t* sk_c  = SK_C(sk);
  const uint8_t* sk_pt = SK_PT(sk);

  picnic_context_t context;
  mzd_from_char_array(context.m_plaintext, sk_pt, LOWMC_BLOCK_SZ);
  mzd_from_char_array(context.m_key, sk_sk, LOWMC_BLOCK_SZ);
  context.plaintext   = sk_pt;
  context.private_key = sk_sk;
  context.public_key  = sk_c;
  context.msg         = message;
  context.msglen      = message_len;
#if defined(WITH_UNRUH)
  context.unruh = (PARAM == Picnic_L1_UR || PARAM == Picnic_L3_UR || PARAM == Picnic_L5_UR);
#endif
  int ret = picnic_impl_sign(instance, &context, signature, signature_len);
  picnic_explicit_bzero(context.m_key, sizeof(context.m_key));
  return ret;
}

int PICNIC_CALLING_CONVENTION picnic_l1_full_verify(const picnic_l1_full_publickey_t* pk,
                                                  const uint8_t* message, size_t message_len,
                                                  const uint8_t* signature, size_t signature_len) {
  if (!pk || !signature || !signature_len) {
    return -1;
  }

  const picnic_instance_t* instance = picnic_instance_get(PARAM);
  if (!instance) {
    return -1;
  }

  const uint8_t* pk_c  = PK_C(pk);
  const uint8_t* pk_pt = PK_PT(pk);

  picnic_context_t context;
  mzd_from_char_array(context.m_plaintext, pk_pt, LOWMC_BLOCK_SZ);
  mzd_from_char_array(context.m_key, pk_c, LOWMC_BLOCK_SZ);
  context.plaintext   = PK_PT(pk);
  context.private_key = NULL;
  context.public_key  = PK_C(pk);
  context.msg         = message;
  context.msglen      = message_len;
#if defined(WITH_UNRUH)
  context.unruh = (PARAM == Picnic_L1_UR || PARAM == Picnic_L3_UR || PARAM == Picnic_L5_UR);
#endif
  return picnic_impl_verify(instance, &context, signature, signature_len);
}

int PICNIC_CALLING_CONVENTION picnic_l1_full_write_public_key(const picnic_l1_full_publickey_t* key,
                                                            uint8_t* buf, size_t buflen) {
  if (!key || !buf) {
    return -1;
  }

  if (buflen < PUBLIC_KEY_SIZE) {
    return -1;
  }

  buf[0] = PARAM;
  memcpy(buf + 1, key->data, PUBLIC_KEY_SIZE - 1);
  return PUBLIC_KEY_SIZE;
}

int PICNIC_CALLING_CONVENTION picnic_l1_full_read_public_key(picnic_l1_full_publickey_t* key,
                                                           const uint8_t* buf, size_t buflen) {
  if (!key || !buf || buflen < 1) {
    return -1;
  }

  const picnic_params_t param = buf[0];
  if (param != PARAM) {
    return -1;
  }
  const picnic_instance_t* instance = picnic_instance_get(PARAM);
  if (!instance) {
    return -1;
  }

  if (buflen < PUBLIC_KEY_SIZE) {
    return -1;
  }

#if (LOWMC_BLOCK_BITS & 0x7) != 0
  static const unsigned int diff = LOWMC_BLOCK_SZ * 8 - LOWMC_BLOCK_BITS;
  if (check_padding_bits(buf[1 + LOWMC_BLOCK_SZ - 1], diff) ||
      check_padding_bits(buf[1 + 2 * LOWMC_BLOCK_SZ - 1], diff)) {
    return -1;
  }
#endif

  memcpy(key->data, buf + 1, PUBLIC_KEY_SIZE - 1);
  return 0;
}

int PICNIC_CALLING_CONVENTION picnic_l1_full_write_private_key(const picnic_l1_full_privatekey_t* key,
                                                             uint8_t* buf, size_t buflen) {
  if (!key || !buf) {
    return -1;
  }

  if (buflen < PRIVATE_KEY_SIZE) {
    return -1;
  }

  buf[0] = PARAM;
  memcpy(buf + 1, &key->data, PRIVATE_KEY_SIZE - 1);
  return PRIVATE_KEY_SIZE;
}

int PICNIC_CALLING_CONVENTION picnic_l1_full_read_private_key(picnic_l1_full_privatekey_t* key,
                                                            const uint8_t* buf, size_t buflen) {
  if (!key || !buf || buflen < 1) {
    return -1;
  }

  const picnic_params_t param = buf[0];
  if (param != PARAM) {
    return -1;
  }

  const picnic_instance_t* instance = picnic_instance_get(PARAM);
  if (!instance) {
    return -1;
  }

  if (buflen < PRIVATE_KEY_SIZE) {
    return -1;
  }

#if (LOWMC_BLOCK_BITS & 0x7) != 0
  static const unsigned int diff = LOWMC_BLOCK_SZ * 8 - LOWMC_BLOCK_BITS;
  /* sanity check of public data: padding bits need to be 0 */
  const int check = check_padding_bits(buf[1 + LOWMC_BLOCK_SZ - 1], diff) |
                    check_padding_bits(buf[1 + 2 * LOWMC_BLOCK_SZ - 1], diff) |
                    check_padding_bits(buf[1 + 3 * LOWMC_BLOCK_SZ - 1], diff);
  picnic_declassify(&check, sizeof(check));
  if (check) {
    return -1;
  }
#endif

  memcpy(key->data, buf + 1, PRIVATE_KEY_SIZE - 1);
  return 0;
}

void PICNIC_CALLING_CONVENTION picnic_l1_full_clear_private_key(picnic_l1_full_privatekey_t* key) {
  picnic_explicit_bzero(key, sizeof(*key));
}
