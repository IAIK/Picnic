/*
 *  This file is part of the optimized implementation of the Picnic signature scheme.
 *  See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */

#ifndef PICNIC3_SIMULATE_MUL_H
#define PICNIC3_SIMULATE_MUL_H

#include "picnic3_types.h"
#include "io.h"

void mpc_matrix_mul_uint64_128(mzd_local_t** output, mzd_local_t** vec, const mzd_local_t* matrix,
                               shares_t* mask_shares);
void mpc_matrix_mul_uint64_129(mzd_local_t** output, mzd_local_t** vec, const mzd_local_t* matrix,
                               shares_t* mask_shares);
void mpc_matrix_mul_uint64_192(mzd_local_t** output, mzd_local_t** vec, const mzd_local_t* matrix,
                               shares_t* mask_shares);
void mpc_matrix_mul_uint64_256(mzd_local_t** output, mzd_local_t** vec, const mzd_local_t* matrix,
                               shares_t* mask_shares);
void mpc_matrix_mul_s128_128(mzd_local_t** output, mzd_local_t** vec, const mzd_local_t* matrix,
                             shares_t* mask_shares);
void mpc_matrix_mul_s128_129(mzd_local_t** output, mzd_local_t** vec, const mzd_local_t* matrix,
                             shares_t* mask_shares);
void mpc_matrix_mul_s128_192(mzd_local_t** output, mzd_local_t** vec, const mzd_local_t* matrix,
                             shares_t* mask_shares);
void mpc_matrix_mul_s128_256(mzd_local_t** output, mzd_local_t** vec, const mzd_local_t* matrix,
                             shares_t* mask_shares);
void mpc_matrix_mul_s256_128(mzd_local_t** output, mzd_local_t** vec, const mzd_local_t* matrix,
                             shares_t* mask_shares);
void mpc_matrix_mul_s256_129(mzd_local_t** output, mzd_local_t** vec, const mzd_local_t* matrix,
                             shares_t* mask_shares);
void mpc_matrix_mul_s256_192(mzd_local_t** output, mzd_local_t** vec, const mzd_local_t* matrix,
                             shares_t* mask_shares);
void mpc_matrix_mul_s256_256(mzd_local_t** output, mzd_local_t** vec, const mzd_local_t* matrix,
                             shares_t* mask_shares);


/* helper functions */
void copyShares(shares_t* dst, shares_t* src);
void xor_byte_array(uint8_t* out, const uint8_t* in1, const uint8_t* in2, uint32_t length);
uint64_t tapesToWord(randomTape_t* tapes);
uint64_t tapesToParityOfWord(randomTape_t* tapes, uint8_t without_last);
void reconstructShares(uint32_t* output, shares_t* shares);
void reconstructSharesN(uint32_t** output, shares_t* shares);
#if defined(PICNIC_STATIC)
void transpose_64_64_lsb(const uint64_t* in, uint64_t* out);
#endif
void transpose_64_64(const uint64_t* in, uint64_t* out);

#endif