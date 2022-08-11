/*
 *  This file is part of the optimized implementation of the Picnic signature scheme.
 *  See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */

#include "picnic.h"

#ifndef Picnic3_L1_H
#define Picnic3_L1_H

#ifdef __cplusplus
extern "C" {
#endif

/** Public key */
typedef struct {
  uint8_t data[PICNIC_PUBLIC_KEY_SIZE_Picnic3_L1 - 1];
} picnic3_l1_publickey_t;

/** Private key */
typedef struct {
  uint8_t data[PICNIC_PRIVATE_KEY_SIZE_Picnic3_L1 - 1];
} picnic3_l1_privatekey_t;

/**
 * Get a string representation of the parameter set.
 *
 * @return A null-terminated string describing the parameter set.
 */
PICNIC_EXPORT const char* PICNIC_CALLING_CONVENTION picnic3_l1_get_param_name(void);

/**
 * Get the size of a private key for serialization
 *
 * @return The size of serialized private key, or 0 on error.
 */
PICNIC_EXPORT size_t PICNIC_CALLING_CONVENTION picnic3_l1_get_private_key_size(void);

/**
 * Get the size of a public key for serialization
 *
 * @return The size of serialized public key, or 0 on error.
 */
PICNIC_EXPORT size_t PICNIC_CALLING_CONVENTION picnic3_l1_get_public_key_size(void);

/* Signature API */

/**
 * Key generation function.
 * Generates a public and private key pair, for the specified parameter set.
 *
 * @param[out] pk         The new public key.
 * @param[out] sk         The new private key.
 *
 * @return Returns 0 for success, or a nonzero value indicating an error.
 */
PICNIC_EXPORT int PICNIC_CALLING_CONVENTION picnic3_l1_keygen(picnic3_l1_publickey_t* pk,
                                                                picnic3_l1_privatekey_t* sk);

/**
 * Signature function.
 * Signs a message with the given keypair.
 *
 * @param[in] sk      The signer's private key.
 * @param[in] message The message to be signed.
 * @param[in] message_len The length of the message, in bytes.
 * @param[out] signature A buffer to hold the signature. The specific max number of
 * bytes required for a parameter set is given by picnic3_l1_signature_size(). Note
 * that the length of each signature varies slightly, for the parameter sets
 * using the FS transform.  The parameter sets using the Unruh transform have a
 * fixed length.
 * @param[in,out] signature_len The length of the provided signature buffer.
 * On success, this is set to the number of bytes written to the signature buffer.
 *
 * @return Returns 0 for success, or a nonzero value indicating an error.
 *
 * @see picnic_verify(), picnic_keygen(), picnic_signature_size()
 */
PICNIC_EXPORT int PICNIC_CALLING_CONVENTION picnic3_l1_sign(const picnic3_l1_privatekey_t* sk,
                                                              const uint8_t* message,
                                                              size_t message_len,
                                                              uint8_t* signature,
                                                              size_t* signature_len);

/**
 * Get the number of bytes required to hold a signature.
 *
 * @param[in] parameters The parameter set of the signature.
 *
 * @return The number of bytes required to hold the signature created by
 * picnic_sign
 *
 * @note The size of signatures with parameter sets using the FS transform vary
 *       slightly based on the random choices made during signing.  This function
 *       will return a suffcient number of bytes to hold a signature, and the
 *       picnic_sign() function returns the exact number used for a given signature.
 *
 * @see picnic_sign()
 */
PICNIC_EXPORT size_t PICNIC_CALLING_CONVENTION picnic3_l1_signature_size();

/**
 * Verification function.
 * Verifies a signature is valid with respect to a public key and message.
 *
 * @param[in] pk      The signer's public key.
 * @param[in] message The message the signature purpotedly signs.
 * @param[in] message_len The length of the message, in bytes.
 * @param[in] signature The signature to verify.
 * @param[in] signature_len The length of the signature.
 *
 * @return Returns 0 for success, indicating a valid signature, or a nonzero
 * value indicating an error or an invalid signature.
 *
 * @see picnic_sign(), picnic_keygen()
 */
PICNIC_EXPORT int PICNIC_CALLING_CONVENTION picnic3_l1_verify(const picnic3_l1_publickey_t* pk,
                                                                const uint8_t* message,
                                                                size_t message_len,
                                                                const uint8_t* signature,
                                                                size_t signature_len);

/**
 * Serialize a public key.
 *
 * @param[in]  key The public key to serialize
 * @param[out] buf The buffer to write the key to.
 *                 Must have size at least PICNIC_MAX_PUBLIC_KEY_SIZE_picnic3_l1 bytes.
 * @param[in]  buflen The length of buf, in bytes
 *
 * @return Returns the number of bytes written.
 */
PICNIC_EXPORT int PICNIC_CALLING_CONVENTION
picnic3_l1_write_public_key(const picnic3_l1_publickey_t* key, uint8_t* buf, size_t buflen);

/**
 * De-serialize a public key.
 *
 * @param[out] key The public key object to be populated.
 * @param[in]  buf The buffer to read the public key from.
 *                 Must be at least PICNIC_MAX_PUBLIC_KEY_SIZE_picnic3_l1 bytes.
 * @param[in]  buflen The length of buf, in bytes
 *
 * @return Returns 0 on success, or a nonzero value indicating an error.
 */
PICNIC_EXPORT int PICNIC_CALLING_CONVENTION
picnic3_l1_read_public_key(picnic3_l1_publickey_t* key, const uint8_t* buf, size_t buflen);

/**
 * Serialize a private key.
 *
 * @param[in]  key The private key to serialize
 * @param[out] buf The buffer to write the key to.
 *                 Must have size at least PICNIC_MAX_PRIVATE_KEY_SIZE_picnic3_l1 bytes.
 * @param[in]  buflen The length of buf, in bytes
 *
 * @return Returns the number of bytes written.
 */
PICNIC_EXPORT int PICNIC_CALLING_CONVENTION
picnic3_l1_write_private_key(const picnic3_l1_privatekey_t* key, uint8_t* buf, size_t buflen);

/**
 * De-serialize a private key.
 *
 * @param[out] key The private key object to be populated
 * @param[in]  buf The buffer to read the key from.
 *                 Must have size at least PICNIC_MAX_PRIVATE_KEY_SIZE_picnic3_l1 bytes.
 * @param[in]  buflen The length of buf, in bytes
 *
 * @return Returns 0 on success, or a nonzero value indicating an error.
 */
PICNIC_EXPORT int PICNIC_CALLING_CONVENTION
picnic3_l1_read_private_key(picnic3_l1_privatekey_t* key, const uint8_t* buf, size_t buflen);

/**
 * Check that a key pair is valid.
 *
 * @param[in] privatekey The private key to check
 * @param[in] publickey The public key to check
 *
 * @return Returns 0 if the key pair is valid, or a nonzero value indicating an error
 */
PICNIC_EXPORT int PICNIC_CALLING_CONVENTION picnic3_l1_validate_keypair(
    const picnic3_l1_privatekey_t* privatekey, const picnic3_l1_publickey_t* publickey);

/**
 * Clear data of a private key.
 *
 * @param[out] key The private key to clear
 */
PICNIC_EXPORT void PICNIC_CALLING_CONVENTION
picnic3_l1_clear_private_key(picnic3_l1_privatekey_t* key);

/**
 * Compute public key from private key.
 *
 * @param[in] privatekey The private key
 * @param[out] publickey The public key to be populated
 * @return Returns 0 on success, or a nonzero value indicating an error.
 **/
PICNIC_EXPORT int PICNIC_CALLING_CONVENTION picnic3_l1_sk_to_pk(
    const picnic3_l1_privatekey_t* privatekey, picnic3_l1_publickey_t* publickey);

#ifdef __cplusplus
}
#endif

#endif