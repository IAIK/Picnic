/*! @file picnic2_impl.c
 *  @brief This is the main file of the signature scheme for the Picnic2
 *  parameter sets.
 *
 *  This file is part of the reference implementation of the Picnic signature scheme.
 *  See the accompanying documentation for complete details.
 *
 *  The code is provided under the MIT license, see LICENSE for
 *  more details.
 *  SPDX-License-Identifier: MIT
 */

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "kdf_shake.h"
#include "macros.h"
#include "picnic_impl.h"
#include "picnic2_impl.h"
#include "picnic.h"
#include "picnic2_types.h"
#include "picnic2_tree.h"
#include "io.h"

#define LOWMC_MAX_STATE_SIZE 64
#define LOWMC_MAX_KEY_BITS 256
#define LOWMC_MAX_AND_GATES (3*38*10 + 4)   /* Rounded to nearest byte */
#define MAX_AUX_BYTES ((LOWMC_MAX_AND_GATES + LOWMC_MAX_KEY_BITS) / 8 + 1)

/* Number of leading zeroes of x.
 * From the book
 * H.S. Warren, *Hacker's Delight*, Pearson Education, 2003.
 * http://www.hackersdelight.org/hdcodetxt/nlz.c.txt
 */

static int32_t nlz(uint32_t x)
{
    uint32_t n;

    if (x == 0) return (32);
    n = 1;
    if((x >> 16) == 0) {n = n + 16; x = x << 16;}
    if((x >> 24) == 0) {n = n + 8;  x = x << 8;}
    if((x >> 28) == 0) {n = n + 4;  x = x << 4;}
    if((x >> 30) == 0) {n = n + 2;  x = x << 2;}
    n = n - (x >> 31);

    return n;
}

/* Helper functions */

void printHex(const char* s, uint8_t* data, size_t len)
{
    printf("%s: ", s);
    for (size_t i = 0; i < len; i++) {
        printf("%02X", data[i]);
    }
    printf("\n");
}

/* Get one bit from a byte array */
uint8_t getBit(const uint8_t* array, uint32_t bitNumber)
{
    return (array[bitNumber / 8] >> (7 - (bitNumber % 8))) & 0x01;
}

/* Get one bit from a 32-bit int array */
uint8_t getBitFromWordArray(const uint32_t* array, uint32_t bitNumber)
{
    return getBit((uint8_t*)array, bitNumber);
}

/* Set a specific bit in a byte array to a given value */
void setBit(uint8_t* bytes, uint32_t bitNumber, uint8_t val)
{
    bytes[bitNumber / 8] = (bytes[bitNumber >> 3]
                            & ~(1 << (7 - (bitNumber % 8)))) | (val << (7 - (bitNumber % 8)));
}

/* Set a specific bit in a byte array to a given value */
void setBitInWordArray(uint32_t* array, uint32_t bitNumber, uint8_t val)
{
    setBit((uint8_t*)array, bitNumber, val);
}

uint8_t parity(uint32_t* data, size_t len)
{
    uint32_t x = data[0];

    for (size_t i = 1; i < len; i++) {
        x ^= data[i];
    }

    /* Compute parity of x using code from Section 5-2 of
     * H.S. Warren, *Hacker's Delight*, Pearson Education, 2003.
     * http://www.hackersdelight.org/hdcodetxt/parity.c.txt
     */
    uint32_t y = x ^ (x >> 1);
    y ^= (y >> 2);
    y ^= (y >> 4);
    y ^= (y >> 8);
    y ^= (y >> 16);
    return y & 1;
}

uint32_t numBytes(uint32_t numBits)
{
    return (numBits == 0) ? 0 : ((numBits - 1) / 8 + 1);
}
uint32_t ceil_log2(uint32_t x)
{
    if (x == 0) {
        return 0;
    }
    return 32 - nlz(x - 1);
}

static void createRandomTapes(randomTape_t* tapes, uint8_t** seeds, uint8_t* salt, size_t t, const picnic_instance_t* params)
{
    Keccak_HashInstancetimes4 ctx;


    size_t tapeSizeBytes = 2 * params->view_size + params->input_size;

    allocateRandomTape(tapes, params);
    assert(params->num_MPC_parties % 4 == 0);
    for (size_t i = 0; i < params->num_MPC_parties; i+=4) {
        hash_init_x4(&ctx, params);

        const uint8_t* seeds_ptr[4] = {seeds[i], seeds[i+1], seeds[i+2], seeds[i+3]};
        hash_update_x4(&ctx, seeds_ptr, params->seed_size);
        const uint8_t* salt_ptr[4] = {salt, salt, salt, salt};
        hash_update_x4(&ctx, salt_ptr, params->seed_size);
        uint16_t tLE = htole16((uint16_t)t);
        const uint8_t* tLE_ptr[4] = {(uint8_t*)&tLE, (uint8_t*)&tLE, (uint8_t*)&tLE, (uint8_t*)&tLE};
        hash_update_x4(&ctx, tLE_ptr, sizeof(uint16_t));
        uint16_t iLE0 = htole16((uint16_t)(i+0));
        uint16_t iLE1 = htole16((uint16_t)(i+1));
        uint16_t iLE2 = htole16((uint16_t)(i+2));
        uint16_t iLE3 = htole16((uint16_t)(i+3));
        const uint8_t* iLE_ptr[4] = {(uint8_t*)&iLE0, (uint8_t*)&iLE1, (uint8_t*)&iLE2, (uint8_t*)&iLE3};
        hash_update_x4(&ctx, iLE_ptr, sizeof(uint16_t));
        hash_final_x4(&ctx);

        uint8_t* out_ptr[4] = {tapes->tape[i], tapes->tape[i+1], tapes->tape[i+2], tapes->tape[i+3]};
        hash_squeeze_x4(&ctx, out_ptr, tapeSizeBytes);
    }
}

#if 0
static uint64_t tapesToWordSimple(randomTape_t* tapes)
{
    uint64_t shares;

    for (size_t i = 0; i < 64; i++) {
        uint8_t bit = getBit(tapes->tape[i], tapes->pos);
        setBit((uint8_t*)&shares, i, bit);
    }
    tapes->pos++;
    return shares;
}
#endif

// TODO: hide behind architecture defines
#include <emmintrin.h>
void sse_trans_64_64(uint64_t const *inp, uint64_t *out) {
    const int nrows = 64;
    const int ncols = 64;
    const uint8_t* I = (uint8_t*) inp;
    uint8_t* O = (uint8_t*) out;
#   define INP(x, y) I[(x)*ncols/8 + (y)/8]
#   define OUT(x, y) O[(y)*nrows/8 + (x)/8]
    int rr, cc, i;
    union {
        __m128i x;
        uint8_t b[16];
    } tmp;

    // Do the main body in 16x8 blocks:
    for (rr = 0; rr <= nrows - 16; rr += 16) {
        for (cc = 0; cc < ncols; cc += 8) {
            for (i = 0; i < 16; ++i)
                tmp.b[i] = INP(rr + i, cc);
            for (i = 8; --i >= 0; tmp.x = _mm_slli_epi64(tmp.x, 1))
                *(uint16_t *) &OUT(rr, cc + i) = _mm_movemask_epi8(tmp.x);
        }
    }
}

static uint64_t tapesToWord(randomTape_t* tapes)
{
    uint64_t shares;

    if(tapes->pos % 64 == 0) {
        uint64_t buffer[64];
        for(size_t i = 0; i < 64; i++) {
            buffer[i/8*8 + 7 - i%8] = ((uint64_t*)tapes->tape[i])[tapes->pos / 64];
        }
        sse_trans_64_64(buffer, tapes->buffer);
    }

    shares = tapes->buffer[(tapes->pos % 64)/8*8 + 7 - (tapes->pos%64)%8];
    tapes->pos++;
    return shares;
}

/* Read one bit from each tape and assemble them into a word.  The tapes form a
 * z by N matrix, we'll transpose it, then the first "count" N-bit rows forms
 * an output word.  In the current implementation N is 64 so the words are
 * uint64_t. The return value must be freed with freeShares().
 */
static void tapesToWords(shares_t* shares, randomTape_t* tapes)
{
    for (size_t w = 0; w < shares->numWords; w++) {
        shares->shares[w] = tapesToWord(tapes);
    }
}

static void copyShares(shares_t* dst, shares_t* src)
{
    assert(dst->numWords == src->numWords);
    memcpy(dst->shares, src->shares, dst->numWords * sizeof(dst->shares[0]));
}

void xor_array(uint32_t* out, const uint32_t * in1, const uint32_t * in2, uint32_t length)
{
    for (uint32_t i = 0; i < length; i++) {
        out[i] = in1[i] ^ in2[i];
    }
}

void xor_array_RC(uint8_t* out, const uint8_t * in1, const uint8_t * in2, uint32_t length)
{
    for (uint32_t i = 0; i < length; i++) {
        out[i] = in1[i] ^ in2[length-1-i];
    }
}

/* For an input bit b = 0 or 1, return the word of all b bits, i.e.,
 * extend(1) = 0xFFFFFFFFFFFFFFFF
 * extend(0) = 0x0000000000000000
 * Assumes inputs are always 0 or 1.  If this doesn't hold, add "& 1" to the
 * input.
 */
ATTR_ALWAYS_INLINE inline uint64_t extend(uint8_t bit)
{
    return ~(bit - 1);
}

static uint64_t aux_mpc_AND(uint64_t a, uint64_t b, randomTape_t* tapes)
{
    uint64_t mask_a = parity64_uint64(a);
    uint64_t mask_b = parity64_uint64(b);
    uint64_t fresh_output_mask = tapesToWord(tapes);

    uint64_t and_helper = tapesToWord(tapes);

    /* Zero the last party's share of the helper value, compute it based on the
     * input masks; then update the tape. */
    setBit((uint8_t*)&and_helper, 63, 0);
    uint64_t aux_bit = (mask_a & mask_b) ^ parity64_uint64(and_helper);
    size_t lastParty = tapes->nTapes - 1;
    setBit(tapes->tape[lastParty], tapes->pos - 1, (uint8_t)aux_bit);

    return fresh_output_mask;
}


/**
 * S-box for m = 10, for Picnic2 aux computation
 */
void sbox_layer_10_uint64_aux(uint64_t* d, randomTape_t* tapes) {
    *d = __bswap_64(*d);
    uint8_t* state = (uint8_t*)d;
    for (uint32_t i = 0; i < 30; i += 3) {
        uint8_t a = getBit(state, i+2);
        uint8_t b = getBit(state, i+1);
        uint8_t c = getBit(state, i+0);

        uint8_t ab = parity64_uint64(aux_mpc_AND(a, b, tapes));
        uint8_t bc = parity64_uint64(aux_mpc_AND(b, c, tapes));
        uint8_t ca = parity64_uint64(aux_mpc_AND(c, a, tapes));

        setBit(state, i+2, a ^ bc);
        setBit(state, i+1, a ^ b ^ ca);
        setBit(state, i+0, a ^ b ^ c ^ ab);
    }
    *d = __bswap_64(*d);
}

static void mpc_xor_masks(shares_t* out, const shares_t* a, const shares_t* b)
{
    assert(out->numWords == a->numWords && a->numWords == b->numWords);

    for (size_t i = 0; i < out->numWords; i++) {
        out->shares[i] = a->shares[i] ^ b->shares[i];
    }
}

static void mpc_xor_masks_nl(shares_t* out, const shares_t* a, const shares_t* b, size_t index, size_t num)
{
    for (size_t i = 0; i < num; i++) {
        out->shares[i] = a->shares[i] ^ b->shares[index + num - 1 - i];
    }
}

#if 0

static void aux_mpc_sbox(shares_t* state, randomTape_t* tapes, const picnic_instance_t* params)
{
    for (size_t i = 0; i < params->lowmc->m * 3; i += 3) {
        uint64_t a = state->shares[i + 2];
        uint64_t b = state->shares[i + 1];
        uint64_t c = state->shares[i];

        uint64_t ab = aux_mpc_AND(a, b, tapes);
        uint64_t bc = aux_mpc_AND(b, c, tapes);
        uint64_t ca = aux_mpc_AND(c, a, tapes);

        state->shares[i + 2] = a ^ bc;
        state->shares[i + 1] = a ^ b ^ ca;
        state->shares[i] = a ^ b ^ c ^ ab;
    }
}

static void aux_matrix_mul(shares_t* output, const shares_t* vec, const uint64_t* matrix, shares_t* tmp_output, const picnic_instance_t* params)
{
    const uint32_t rowstride = ((params->lowmc->n+127)/128*128)/8;
    memset(tmp_output->shares, 0, sizeof(uint64_t)*tmp_output->numWords);
    for (size_t i = 0; i < params->lowmc->n; i++) {
        for (uint32_t j = 0; j < params->lowmc->n; j+=8) {
            uint8_t matrix_byte = ((uint8_t*)matrix)[i * rowstride + (params->lowmc->n -1 - j) / 8];
            tmp_output->shares[j+0] ^= vec->shares[params->lowmc->n - 1 - i] & extend((matrix_byte >> 7) & 1);
            tmp_output->shares[j+1] ^= vec->shares[params->lowmc->n - 1 - i] & extend((matrix_byte >> 6) & 1);
            tmp_output->shares[j+2] ^= vec->shares[params->lowmc->n - 1 - i] & extend((matrix_byte >> 5) & 1);
            tmp_output->shares[j+3] ^= vec->shares[params->lowmc->n - 1 - i] & extend((matrix_byte >> 4) & 1);
            tmp_output->shares[j+4] ^= vec->shares[params->lowmc->n - 1 - i] & extend((matrix_byte >> 3) & 1);
            tmp_output->shares[j+5] ^= vec->shares[params->lowmc->n - 1 - i] & extend((matrix_byte >> 2) & 1);
            tmp_output->shares[j+6] ^= vec->shares[params->lowmc->n - 1 - i] & extend((matrix_byte >> 1) & 1);
            tmp_output->shares[j+7] ^= vec->shares[params->lowmc->n - 1 - i] & extend((matrix_byte >> 0) & 1);
        }
    }

    copyShares(output, tmp_output);
}

// 1x128 x 128x30
static void aux_matrix_mul_z(shares_t* output, const shares_t* vec, const uint64_t* matrix, shares_t* tmp_output, const picnic_instance_t* params)
{
    const uint32_t rowstride = ((params->lowmc->n+127)/128*128)/8;
    memset(tmp_output->shares, 0, sizeof(uint64_t) * tmp_output->numWords);
    for (size_t i = 0; i < params->lowmc->m*3; i++) {
        uint64_t new_mask_i = 0;
        for (uint32_t j = 0; j < params->lowmc->n / 8; j++) {
            uint8_t matrix_byte = ((uint8_t*)matrix)[i * rowstride + (params->lowmc->n/8) - 1 - j];
            new_mask_i ^= vec->shares[j * 8 + 0] & extend((matrix_byte >> 7) & 1);
            new_mask_i ^= vec->shares[j * 8 + 1] & extend((matrix_byte >> 6) & 1);
            new_mask_i ^= vec->shares[j * 8 + 2] & extend((matrix_byte >> 5) & 1);
            new_mask_i ^= vec->shares[j * 8 + 3] & extend((matrix_byte >> 4) & 1);
            new_mask_i ^= vec->shares[j * 8 + 4] & extend((matrix_byte >> 3) & 1);
            new_mask_i ^= vec->shares[j * 8 + 5] & extend((matrix_byte >> 2) & 1);
            new_mask_i ^= vec->shares[j * 8 + 6] & extend((matrix_byte >> 1) & 1);
            new_mask_i ^= vec->shares[j * 8 + 7] & extend((matrix_byte >> 0) & 1);
        }
        tmp_output->shares[params->lowmc->m*3-1-i] = new_mask_i;
    }

    copyShares(output, tmp_output);
}

static void aux_shuffle(shares_t* state, uint64_t r_mask) {
    for(int i = 63; i >= 0 && r_mask != 0xFFFFFFFC00000000UL; i--) {
        if(!((r_mask >> i) & 1)) { // bit is not set
           //find next 1 and swap all entries until then
           for(int j = i-1; j >= 0; j--) {
               if((r_mask >> j) & 1) {
                  for(int k = j; k < i; k++) {
                      uint64_t t = state->shares[63-k];
                      state->shares[63-k] = state->shares[63-k-1];
                      state->shares[63-k-1] = t;
                  }
                  r_mask |= (1UL << i);  // set bit i
                  r_mask &= ~(1UL << j); // clear bit j
                  break;
               }
           }
        }
    }
}

//1x30 x 30x98
static void aux_matrix_addmul_r(shares_t* output, const shares_t* vec, const uint64_t* matrix, shares_t* tmp_output, const picnic_instance_t* params)
{
    const uint32_t rowstride = ((params->lowmc->n+127)/128*128)/8;
    copyShares(tmp_output, output);
    for (size_t i = 0; i < params->lowmc->m*3; i++) {
        for (uint32_t j = 0; j < params->lowmc->n; j+=8) {
            uint8_t matrix_byte = ((uint8_t*)matrix)[i * rowstride + (params->lowmc->n -1 - j) / 8];
            tmp_output->shares[j+0] ^= vec->shares[params->lowmc->m*3 - 1 - i] & extend((matrix_byte >> 7) & 1);
            tmp_output->shares[j+1] ^= vec->shares[params->lowmc->m*3 - 1 - i] & extend((matrix_byte >> 6) & 1);
            tmp_output->shares[j+2] ^= vec->shares[params->lowmc->m*3 - 1 - i] & extend((matrix_byte >> 5) & 1);
            tmp_output->shares[j+3] ^= vec->shares[params->lowmc->m*3 - 1 - i] & extend((matrix_byte >> 4) & 1);
            tmp_output->shares[j+4] ^= vec->shares[params->lowmc->m*3 - 1 - i] & extend((matrix_byte >> 3) & 1);
            tmp_output->shares[j+5] ^= vec->shares[params->lowmc->m*3 - 1 - i] & extend((matrix_byte >> 2) & 1);
            tmp_output->shares[j+6] ^= vec->shares[params->lowmc->m*3 - 1 - i] & extend((matrix_byte >> 1) & 1);
            tmp_output->shares[j+7] ^= vec->shares[params->lowmc->m*3 - 1 - i] & extend((matrix_byte >> 0) & 1);
        }
    }
    copyShares(output, tmp_output);
}

//assumes nl_part to be zero'd
static void aux_matrix_mul_nl_part(shares_t* nl_part, const shares_t* key, const uint64_t* matrix, const picnic_instance_t* params)
{
    const uint32_t rowstride = ((params->lowmc->r*32+255)/256*256)/8;
    for (size_t i = 0; i < params->lowmc->n; i++) {
        for (uint32_t j = 0; j < params->lowmc->r*32; j+=8) {
            uint8_t matrix_byte = ((uint8_t*)matrix)[i * rowstride + j/8];
            nl_part->shares[j+0] ^= key->shares[params->lowmc->n - 1 - i] & extend((matrix_byte >> 0) & 1);
            nl_part->shares[j+1] ^= key->shares[params->lowmc->n - 1 - i] & extend((matrix_byte >> 1) & 1);
            nl_part->shares[j+2] ^= key->shares[params->lowmc->n - 1 - i] & extend((matrix_byte >> 2) & 1);
            nl_part->shares[j+3] ^= key->shares[params->lowmc->n - 1 - i] & extend((matrix_byte >> 3) & 1);
            nl_part->shares[j+4] ^= key->shares[params->lowmc->n - 1 - i] & extend((matrix_byte >> 4) & 1);
            nl_part->shares[j+5] ^= key->shares[params->lowmc->n - 1 - i] & extend((matrix_byte >> 5) & 1);
            nl_part->shares[j+6] ^= key->shares[params->lowmc->n - 1 - i] & extend((matrix_byte >> 6) & 1);
            nl_part->shares[j+7] ^= key->shares[params->lowmc->n - 1 - i] & extend((matrix_byte >> 7) & 1);
        }
    }
}

#if 0
/* Simpler version of aux_matrix_mul, closer to the description in the spec */
static void aux_matrix_mul_simple(shares_t* output, const shares_t* vec, const uint32_t* matrix, shares_t* tmp_output, const picnic_instance_t* params)
{
    for (size_t i = 0; i < params->lowmc->n; i++) {

        uint64_t new_mask_i = 0;
        for (uint32_t j = 0; j < params->lowmc->n; j++) {
            uint8_t matrix_bit = getBit((uint8_t*)matrix, i * params->lowmc->n + j);
            new_mask_i ^= vec->shares[j] & extend(matrix_bit);
        }
        tmp_output->shares[i] = new_mask_i;
    }

    copyShares(output, tmp_output);
}
#endif


/* Input is the tapes for one parallel repitition; i.e., tapes[t]
 * Updates the random tapes of all players with the mask values for the output of
 * AND gates, and computes the N-th party's share such that the AND gate invariant
 * holds on the mask values.
 */
static void computeAuxTapeSimple(randomTape_t* tapes, const picnic_instance_t* params)
{
    shares_t* state = allocateShares(params->lowmc->n);
    shares_t* state2 = allocateShares(params->lowmc->n);
    shares_t* roundKey = allocateShares(params->lowmc->n);
    shares_t* key = allocateShares(params->lowmc->n);
    shares_t* tmp1 = allocateShares(params->lowmc->n);
    shares_t* nl_part = allocateShares(params->lowmc->r * (32));

    tapesToWords(key, tapes);
#if defined(REDUCED_ROUND_KEY_COMPUTATION)
    // The next line is the combination of two operations, simplified because XORs by
    // a constant is a NOP during preprocssing.
    // roundKey = key * KMatrix[0]
    // state = roundKey + plaintext
    aux_matrix_mul(state, key, params->lowmc->k0_matrix->w64, tmp1, params);
    aux_matrix_mul_nl_part(nl_part, key, params->lowmc->precomputed_non_linear_part_matrix, params);

#if defined(OPTIMIZED_LINEAR_LAYER_EVALUATION)
    for (uint32_t r = 0; r < params->lowmc->r-1; r++) {
        aux_mpc_sbox(state, tapes, params);
        // NOP: state =  state XOR RConstant(r - 1)
        mpc_xor_masks_nl(state, state, nl_part, r* 32+2, params->lowmc->m * 3);

        //MUL_Z
        aux_matrix_mul_z(state2, state, params->lowmc->rounds[r].z_matrix->w64, tmp1, params); // state = state * LMatrix(r-1)
        //SHUFFLE state
        aux_shuffle(state, params->lowmc->rounds[r].r_mask);
        //ADDMUL_R
        aux_matrix_addmul_r(state2, state, params->lowmc->rounds[r].r_matrix->w64, tmp1, params); // state = state * LMatrix(r-1)
        //XOR TOGETHER
        for(uint32_t i = 0; i < params->lowmc->m*3; i++) {
            state->shares[i] = 0;
        }
        mpc_xor_masks(state, state, state2);
    }
    aux_mpc_sbox(state, tapes, params);
    //TODO: can we skip this in aux tape computation? seems like
    mpc_xor_masks_nl(state, state, nl_part, (params->lowmc->r-1)* 32+2, params->lowmc->m * 3);
    aux_matrix_mul(state, state, params->lowmc->zr_matrix->w64, tmp1, params); // state = state * ZR_matrix
#else
    for (uint32_t r = 0; r < params->lowmc->r; r++) {
        aux_mpc_sbox(state, tapes, params);
        // NOP: state =  state XOR RConstant(r - 1)
        mpc_xor_masks_nl(state, state, nl_part, r* 32+2, params->lowmc->m * 3);
        aux_matrix_mul(state, state, params->lowmc->rounds[r].l_matrix->w64, tmp1, params); // state = state * LMatrix(r-1)
    }
#endif
#else
    // The next line is the combination of two operations, simplified because XORs by
    // a constant is a NOP during preprocssing.
    // roundKey = key * KMatrix[0]
    // state = roundKey + plaintext
    aux_matrix_mul(state, key, params->lowmc->k0_matrix->w64, tmp1, params);

    for (uint32_t r = 0; r < params->lowmc->r; r++) {
        aux_matrix_mul(roundKey, key, params->lowmc->rounds[r].k_matrix->w64, tmp1, params);    // roundKey = key * KMatrix(r)
        aux_mpc_sbox(state, tapes, params);
        aux_matrix_mul(state, state, params->lowmc->rounds[r].l_matrix->w64, tmp1, params); // state = state * LMatrix(r-1)
        // NOP: state =  state XOR RConstant(r - 1)
        mpc_xor_masks(state, state, roundKey);
    }
#endif
    // Reset the random tape counter so that the online execution uses the
    // same random bits as when computing the aux shares
    tapes->pos = 0;

    freeShares(tmp1);
    freeShares(key);
    freeShares(roundKey);
    freeShares(state);
    freeShares(state2);
    freeShares(nl_part);
}
#endif

/* Input is the tapes for one parallel repitition; i.e., tapes[t]
 * Updates the random tapes of all players with the mask values for the output of
 * AND gates, and computes the N-th party's share such that the AND gate invariant
 * holds on the mask values.
 */
static void computeAuxTape(randomTape_t* tapes, const picnic_instance_t* params) {
    shares_t *key = allocateShares(params->lowmc->n);
    mzd_local_t* lowmc_key = mzd_local_init_ex(params->lowmc->n, 1, true);

    tapesToWords(key, tapes);

    uint8_t temp[32] = {0,};
    // combine into key shares and calculate lowmc evaluation in plain
    for(uint32_t i = 0; i < params->lowmc->n; i++) {
        uint8_t key_bit = parity64_uint64(key->shares[i]);
        setBit(temp, i, key_bit);
    }
    mzd_from_char_array(lowmc_key, temp, params->lowmc->n/8);

    lowmc_compute_aux_implementation_f lowmc_aux_impl = params->impls.lowmc_aux;
    // Perform LowMC evaluation and record state before AND gates
    lowmc_aux_impl(lowmc_key, tapes);

    // Reset the random tape counter so that the online execution uses the
    // same random bits as when computing the aux shares
    tapes->pos = 0;

    freeShares(key);
    mzd_local_free(lowmc_key);
}

static void commit(uint8_t* digest, uint8_t* seed, uint8_t* aux, uint8_t* salt, size_t t, size_t j, const picnic_instance_t* params)
{
    /* Compute C[t][j];  as digest = H(seed||[aux]) aux is optional */
    Keccak_HashInstance ctx;

    hash_init(&ctx, params);
    hash_update(&ctx, seed, params->seed_size);
    if (aux != NULL) {
        size_t tapeLenBytes = params->view_size;
        hash_update(&ctx, aux, tapeLenBytes);
    }
    hash_update(&ctx, salt, params->seed_size);
    uint16_t tLE = htole16((uint16_t)t);
    hash_update(&ctx, (uint8_t*)&tLE, sizeof(uint16_t));
    uint16_t jLE = htole16((uint16_t)j);
    hash_update(&ctx, (uint8_t*)&jLE, sizeof(uint16_t));
    hash_final(&ctx);
    hash_squeeze(&ctx, digest, params->digest_size);
}

static void commit_h(uint8_t* digest, commitments_t* C, const picnic_instance_t* params)
{
    Keccak_HashInstance ctx;

    hash_init(&ctx, params);
    for (size_t i = 0; i < params->num_MPC_parties; i++) {
        hash_update(&ctx, C->hashes[i], params->seed_size);
    }
    hash_final(&ctx);
    hash_squeeze(&ctx, digest, params->digest_size);
}

// Commit to the views for one parallel rep
static void commit_v(uint8_t* digest, uint8_t* input, msgs_t* msgs, const picnic_instance_t* params)
{
    Keccak_HashInstance ctx;

    hash_init(&ctx, params);
    hash_update(&ctx, input, params->input_size);
    for (size_t i = 0; i < params->num_MPC_parties; i++) {
        size_t msgs_size = numBytes(msgs->pos);
        hash_update(&ctx, msgs->msgs[i], msgs_size);
    }
    hash_final(&ctx);
    hash_squeeze(&ctx, digest, params->digest_size);
}

static void reconstructShares(uint32_t* output, shares_t* shares)
{
    for (size_t i = 0; i < shares->numWords; i++) {
        setBitInWordArray(output, i, parity64_uint64(shares->shares[i]));
    }
}

static void wordToMsgsNoTranspose(uint64_t w, msgs_t* msgs)
{
    ((uint64_t*)msgs->msgs[msgs->pos % 64])[msgs->pos / 64] = w;
    msgs->pos++;
}

static void msgsTranspose(msgs_t* msgs)
{
    uint64_t buffer_in[64];
    uint64_t buffer_out[64];
    size_t pos;
    for(pos = 0; pos < msgs->pos/64; pos++) {
        for(size_t i = 0; i < 64; i++) {
            buffer_in[i/8*8 + 7 - i%8] = ((uint64_t*)msgs->msgs[i])[pos];
        }
        sse_trans_64_64(buffer_in, buffer_out);
        for(size_t i = 0; i < 64; i++) {
            ((uint64_t*)msgs->msgs[i])[pos] = buffer_out[(i)/8*8 + 7 - (i)%8];
        }
    }
    memset(&buffer_in, 0, 64*sizeof(uint64_t));
    for(size_t i = 0; i < msgs->pos%64; i++) {
        buffer_in[i/8*8 + 7 - i%8] = ((uint64_t*)msgs->msgs[i])[pos];
    }
    sse_trans_64_64(buffer_in, buffer_out);
    for(size_t i = 0; i < 64; i++) {
        ((uint64_t*)msgs->msgs[i])[pos] = buffer_out[(i)/8*8 + 7 - (i)%8];
    }
}

#if 0
static void wordToMsgs(uint64_t w, msgs_t* msgs)
{
    for (size_t i = 0; i < 64; i++) {
        uint8_t w_i = getBit((uint8_t*)&w, i);
        setBit(msgs->msgs[i], msgs->pos, w_i);
    }
    msgs->pos++;
}
#endif

static uint8_t mpc_AND(uint8_t a, uint8_t b, uint64_t mask_a, uint64_t mask_b, randomTape_t* tapes, msgs_t* msgs,
        uint64_t* out, uint8_t* unopened_msg)
{
    uint64_t output_mask = tapesToWord(tapes); // A fresh random mask to hide the result

    *out = output_mask;
    uint64_t and_helper = tapesToWord(tapes);   // The special mask value setup during preprocessing for each AND gate
    uint64_t s_shares = (extend(a) & mask_b) ^ (extend(b) & mask_a) ^ and_helper ^ output_mask;

    if (msgs->unopened >= 0) {
        uint8_t unopenedPartyBit = getBit(unopened_msg, msgs->pos);
        setBit((uint8_t*)&s_shares, msgs->unopened, unopenedPartyBit);
    }

    // Broadcast each share of s
    wordToMsgsNoTranspose(s_shares, msgs);

    return (uint8_t)(parity64_uint64(s_shares) ^ (a & b));
}

static void mpc_sbox(uint32_t* state, shares_t* state_masks, randomTape_t* tapes, msgs_t* msgs, uint8_t* unopenened_msg, const picnic_instance_t* params)
{
    for (size_t i = 0; i < params->lowmc->m * 3; i += 3) {
        uint8_t a = getBitFromWordArray(state, i + 2);
        uint64_t mask_a = state_masks->shares[i + 2];

        uint8_t b = getBitFromWordArray(state, i + 1);
        uint64_t mask_b = state_masks->shares[i + 1];

        uint8_t c = getBitFromWordArray(state, i);
        uint64_t mask_c = state_masks->shares[i];

        uint64_t bc_mask, ab_mask, ca_mask; // Fresh output masks used for the AND gate

        uint8_t ab = mpc_AND(a, b, mask_a, mask_b, tapes, msgs, &ab_mask, unopenened_msg);
        uint8_t bc = mpc_AND(b, c, mask_b, mask_c, tapes, msgs, &bc_mask, unopenened_msg);
        uint8_t ca = mpc_AND(c, a, mask_c, mask_a, tapes, msgs, &ca_mask, unopenened_msg);

        setBitInWordArray(state, i + 2, a ^ bc);
        state_masks->shares[i + 2] = mask_a ^ bc_mask;
        setBitInWordArray(state, i + 1, a ^ b ^ ca);
        state_masks->shares[i + 1] = mask_b ^ mask_a ^ ca_mask;
        setBitInWordArray(state, i, a ^ b ^ c ^ ab);
        state_masks->shares[i] = mask_a ^ mask_b ^ mask_c ^ ab_mask;
    }
}

/* For each word in shares; write player i's share to their stream of msgs */
static void broadcast(shares_t* shares, msgs_t* msgs)
{
    for (size_t w = 0; w < shares->numWords; w++) {
        wordToMsgsNoTranspose(shares->shares[w], msgs);
    }
}

static void mpc_matrix_mul(uint32_t* output, const uint32_t* vec, const uint64_t* matrix, shares_t* mask_shares, const picnic_instance_t* params)
{
    uint8_t temp[LOWMC_MAX_STATE_SIZE] = {0, };

    const uint32_t rowstride = ((params->lowmc->n+127)/128*128)/8;
    shares_t* tmp_mask = allocateShares(mask_shares->numWords);

    for (size_t i = 0; i < params->lowmc->n; i++) {
        uint8_t vec_bit = extend(getBit((uint8_t*)vec, params->lowmc->n - 1 - i))  & 0xFF;

        for (uint32_t j = 0; j < params->lowmc->n; j+=8) {
            uint8_t matrix_byte = ((uint8_t*)matrix)[(i * rowstride) + (params->lowmc->n - 1 - j) / 8];
            temp[j/8] ^= matrix_byte & vec_bit;

            tmp_mask->shares[j+0] ^= mask_shares->shares[params->lowmc->n - 1 - i] & extend((matrix_byte >> 7) & 1);
            tmp_mask->shares[j+1] ^= mask_shares->shares[params->lowmc->n - 1 - i] & extend((matrix_byte >> 6) & 1);
            tmp_mask->shares[j+2] ^= mask_shares->shares[params->lowmc->n - 1 - i] & extend((matrix_byte >> 5) & 1);
            tmp_mask->shares[j+3] ^= mask_shares->shares[params->lowmc->n - 1 - i] & extend((matrix_byte >> 4) & 1);
            tmp_mask->shares[j+4] ^= mask_shares->shares[params->lowmc->n - 1 - i] & extend((matrix_byte >> 3) & 1);
            tmp_mask->shares[j+5] ^= mask_shares->shares[params->lowmc->n - 1 - i] & extend((matrix_byte >> 2) & 1);
            tmp_mask->shares[j+6] ^= mask_shares->shares[params->lowmc->n - 1 - i] & extend((matrix_byte >> 1) & 1);
            tmp_mask->shares[j+7] ^= mask_shares->shares[params->lowmc->n - 1 - i] & extend((matrix_byte >> 0) & 1);

        }
    }
    memcpy(output, temp, params->lowmc->n/8);

    copyShares(mask_shares, tmp_mask);
    freeShares(tmp_mask);
}

static void mpc_matrix_mul_z(uint32_t* state2, const uint32_t* state, shares_t* mask2_shares, const shares_t* mask_shares, const uint64_t* matrix, const picnic_instance_t* params)
{
    const uint32_t rowstride = ((params->lowmc->n+127)/128*128)/8;
    memset(mask2_shares->shares, 0, sizeof(uint64_t) * mask2_shares->numWords);
    memset(state2, 0, params->lowmc->n/8);
    for (size_t i = 0; i < params->lowmc->m*3; i++) {
        uint8_t prod = 0;
        uint64_t new_mask_i = 0;
        for (uint32_t j = 0; j < params->lowmc->n / 8; j++) {
            uint8_t matrix_byte = ((uint8_t*)matrix)[i * rowstride + (params->lowmc->n/8) - 1 - j];
            uint8_t vec_byte = ((uint8_t*)state)[j];

            prod ^= matrix_byte & vec_byte;

            new_mask_i ^= mask_shares->shares[j * 8 + 0] & extend((matrix_byte >> 7) & 1);
            new_mask_i ^= mask_shares->shares[j * 8 + 1] & extend((matrix_byte >> 6) & 1);
            new_mask_i ^= mask_shares->shares[j * 8 + 2] & extend((matrix_byte >> 5) & 1);
            new_mask_i ^= mask_shares->shares[j * 8 + 3] & extend((matrix_byte >> 4) & 1);
            new_mask_i ^= mask_shares->shares[j * 8 + 4] & extend((matrix_byte >> 3) & 1);
            new_mask_i ^= mask_shares->shares[j * 8 + 5] & extend((matrix_byte >> 2) & 1);
            new_mask_i ^= mask_shares->shares[j * 8 + 6] & extend((matrix_byte >> 1) & 1);
            new_mask_i ^= mask_shares->shares[j * 8 + 7] & extend((matrix_byte >> 0) & 1);
        }
        //byte parity from: https://graphics.stanford.edu/~seander/bithacks.html#ParityWith64Bits
        uint8_t parity = (((prod * 0x0101010101010101ULL) & 0x8040201008040201ULL) % 0x1FF) & 1;
        setBit((uint8_t*)state2, params->lowmc->m*3-1-i, parity);
        mask2_shares->shares[params->lowmc->m*3-1-i] = new_mask_i;
    }
}

static void mpc_shuffle(uint8_t* state, shares_t* mask_shares, uint64_t r_mask) {
    for(int i = 63; i >= 0 && r_mask != 0xFFFFFFFC00000000UL; i--) {
        if(!((r_mask >> i) & 1)) { // bit is not set
            //find next 1 and swap all entries until then
            for(int j = i-1; j >= 0; j--) {
                if((r_mask >> j) & 1) {
                    for(int k = j; k < i; k++) {
                        uint64_t t = mask_shares->shares[63-k];
                        mask_shares->shares[63-k] = mask_shares->shares[63-k-1];
                        mask_shares->shares[63-k-1] = t;

                        uint8_t bit = getBit(state, 63-k);
                        setBit(state, 63-k, getBit(state, 63-k-1));
                        setBit(state, 63-k-1, bit);
                    }
                    r_mask |= (1UL << i);  // set bit i
                    r_mask &= ~(1UL << j); // clear bit j
                    break;
                }
            }
        }
    }
}

static void mpc_matrix_addmul_r(uint32_t* state2, const uint32_t* state, shares_t* mask2_shares,
        shares_t* mask_shares, const uint64_t* matrix, const picnic_instance_t* params)
{
    uint8_t temp[LOWMC_MAX_STATE_SIZE] = {0, };
    memcpy(temp, state2, params->lowmc->n/8);

    const uint32_t rowstride = ((params->lowmc->n+127)/128*128)/8;
    shares_t* tmp_mask = allocateShares(mask_shares->numWords);
    copyShares(tmp_mask, mask2_shares);

    for (size_t i = 0; i < params->lowmc->m*3; i++) {
        uint8_t vec_bit = extend(getBit((uint8_t*)state, params->lowmc->m*3 - 1 - i))  & 0xFF;

        for (uint32_t j = 0; j < params->lowmc->n; j+=8) {
            uint8_t matrix_byte = ((uint8_t*)matrix)[(i * rowstride) + (params->lowmc->n - 1 - j) / 8];

            tmp_mask->shares[j+0] ^= mask_shares->shares[params->lowmc->m*3 - 1 - i] & extend((matrix_byte >> 7) & 1);
            tmp_mask->shares[j+1] ^= mask_shares->shares[params->lowmc->m*3 - 1 - i] & extend((matrix_byte >> 6) & 1);
            tmp_mask->shares[j+2] ^= mask_shares->shares[params->lowmc->m*3 - 1 - i] & extend((matrix_byte >> 5) & 1);
            tmp_mask->shares[j+3] ^= mask_shares->shares[params->lowmc->m*3 - 1 - i] & extend((matrix_byte >> 4) & 1);
            tmp_mask->shares[j+4] ^= mask_shares->shares[params->lowmc->m*3 - 1 - i] & extend((matrix_byte >> 3) & 1);
            tmp_mask->shares[j+5] ^= mask_shares->shares[params->lowmc->m*3 - 1 - i] & extend((matrix_byte >> 2) & 1);
            tmp_mask->shares[j+6] ^= mask_shares->shares[params->lowmc->m*3 - 1 - i] & extend((matrix_byte >> 1) & 1);
            tmp_mask->shares[j+7] ^= mask_shares->shares[params->lowmc->m*3 - 1 - i] & extend((matrix_byte >> 0) & 1);

            //matrix_byte = ((matrix_byte * 0x80200802ULL) & 0x0884422110ULL) * 0x0101010101ULL >> 32;
            temp[j/8] ^= matrix_byte & vec_bit;

        }
    }
    memcpy(state2, temp, params->lowmc->n/8);

    copyShares(mask2_shares, tmp_mask);
    freeShares(tmp_mask);
}

static void mpc_matrix_mul_nl_part(uint32_t* nl_part, const uint32_t* key, const uint64_t* precomputed_nl_matrix,
        const uint64_t* precomputed_constant_nl, const shares_t* nl_part_masks, const shares_t* key_masks, const picnic_instance_t* params)
{

    const uint32_t rowstride = ((params->lowmc->r*32+255)/256*256)/8;
    uint8_t temp[rowstride];
    memset(temp, 0, rowstride);
    for (size_t i = 0; i < params->lowmc->n; i++) {
        uint8_t key_bit = extend(getBit((uint8_t*)key, params->lowmc->n - 1 - i))  & 0xFF;
        for (uint32_t j = 0; j < params->lowmc->r * 32; j+=8) {

            uint8_t matrix_byte = ((uint8_t*)precomputed_nl_matrix)[i * rowstride + j/ 8];
            temp[j/8] ^= matrix_byte & key_bit;

            nl_part_masks->shares[j+0] ^= key_masks->shares[params->lowmc->n - 1 - i] & extend((matrix_byte >> 0) & 1);
            nl_part_masks->shares[j+1] ^= key_masks->shares[params->lowmc->n - 1 - i] & extend((matrix_byte >> 1) & 1);
            nl_part_masks->shares[j+2] ^= key_masks->shares[params->lowmc->n - 1 - i] & extend((matrix_byte >> 2) & 1);
            nl_part_masks->shares[j+3] ^= key_masks->shares[params->lowmc->n - 1 - i] & extend((matrix_byte >> 3) & 1);
            nl_part_masks->shares[j+4] ^= key_masks->shares[params->lowmc->n - 1 - i] & extend((matrix_byte >> 4) & 1);
            nl_part_masks->shares[j+5] ^= key_masks->shares[params->lowmc->n - 1 - i] & extend((matrix_byte >> 5) & 1);
            nl_part_masks->shares[j+6] ^= key_masks->shares[params->lowmc->n - 1 - i] & extend((matrix_byte >> 6) & 1);
            nl_part_masks->shares[j+7] ^= key_masks->shares[params->lowmc->n - 1 - i] & extend((matrix_byte >> 7) & 1);
        }
    }
    xor_array(nl_part, (uint32_t*)temp, (uint32_t*)precomputed_constant_nl, params->lowmc->r);
}

#if 0
/* Alternative, simpler implementation of mpc_matrix_mul, closer to the description in the spec */
static void mpc_matrix_mul_simple(uint32_t* output, const uint32_t* vec, const uint32_t* matrix, shares_t* mask_shares, const picnic_instance_t* params)
{
    uint32_t prod[LOWMC_MAX_STATE_SIZE];
    uint32_t temp[LOWMC_MAX_STATE_SIZE];

    shares_t* tmp_mask = allocateShares(mask_shares->numWords);

    for (size_t i = 0; i < params->lowmc->n; i++) {
        tmp_mask->shares[i] = 0;
        for (uint32_t j = 0; j < params->lowmc->n; j++) {
            uint8_t matrix_bit = getBit((uint8_t*)matrix, i * params->lowmc->n + j);
            uint8_t vec_bit = getBit((uint8_t*)vec, j);
            setBit((uint8_t*)prod, j, matrix_bit & vec_bit);
            tmp_mask->shares[i] ^= mask_shares->shares[j] & extend(matrix_bit);
        }
        uint8_t output_bit_i = parity(&prod[0], (params->input_size / 4));
        setBit((uint8_t*)temp, i, output_bit_i);
    }

    memcpy(output, &temp, params->input_size);
    copyShares(mask_shares, tmp_mask);
    freeShares(tmp_mask);
}
#endif

static void mpc_xor2(uint32_t* output, shares_t* output_masks, const uint32_t* x,
                     const shares_t* x_masks,  const uint32_t* y, const shares_t* y_masks, const picnic_instance_t* params)
{
    xor_array(output, x, y, (params->input_size / 4));
    mpc_xor_masks(output_masks, x_masks, y_masks);
}

static void mpc_xor2_nl(uint32_t* output, shares_t* output_masks, const uint32_t* x,
                     const shares_t* x_masks, const uint32_t* y, const shares_t* y_masks,
                     size_t index, size_t num)
{
    xor_array_RC((uint8_t*)output, (uint8_t*)x, (uint8_t*)&y[index/32], 4);
    //xor masks
    mpc_xor_masks_nl(output_masks, x_masks, y_masks, index, num);
}

#if 0
/* Helper function when debugging MPC function that operate on masked values */
static void print_unmasked(char* label, uint32_t* state, shares_t* mask_shares, const picnic_instance_t* params)
{
    uint32_t tmp[LOWMC_MAX_STATE_SIZE];

    reconstructShares(tmp, mask_shares);
    xor_array(tmp, tmp, state, (params->input_size / 4));
    printHex(label, (uint8_t*)tmp, params->input_size);
}
#endif



static int contains(uint16_t* list, size_t len, size_t value)
{
    for (size_t i = 0; i < len; i++) {
        if (list[i] == value) {
            return 1;
        }
    }
    return 0;
}

static int indexOf(uint16_t* list, size_t len, size_t value)
{
    for (size_t i = 0; i < len; i++) {
        if (list[i] == value) {
            return i;
        }
    }
    assert(!"indexOf called on list where value is not found. (caller bug)");
    return -1;
}

static void getAuxBits(uint8_t* output, randomTape_t* tapes, const picnic_instance_t* params)
{
    size_t firstAuxIndex = params->lowmc->n + 1;
    size_t last = params->num_MPC_parties - 1;
    size_t pos = 0;

    memset(output, 0, params->view_size);
    size_t andSizeBits = 3 * params->lowmc->r * params->lowmc->m;
    for (size_t i = 0; i < andSizeBits * 2; i += 2) {
        uint8_t auxBit = getBit(tapes->tape[last], firstAuxIndex + i);
        setBit(output, pos, auxBit);
        pos++;
    }
}

static void setAuxBits(randomTape_t* tapes, uint8_t* input, const picnic_instance_t* params)
{
    size_t firstAuxIndex = params->lowmc->n + 1;
    size_t last = params->num_MPC_parties - 1;
    size_t pos = 0;

    for (size_t i = 0; i < params->view_size * 2 * 8; i += 2) {
        uint8_t auxBit = getBit(input, pos);
        setBit(tapes->tape[last], firstAuxIndex + i, auxBit);
        pos++;
    }
}

static int simulateOnline(uint32_t* maskedKey, shares_t* mask_shares, randomTape_t*
                          tapes, msgs_t* msgs, const uint32_t* plaintext, const uint32_t* pubKey, const picnic_instance_t* params)
{
    int ret = 0;
    uint32_t* roundKey = malloc(params->input_size);
    uint32_t* state = malloc(params->input_size);
    uint32_t* state2 = malloc(params->input_size);
    uint32_t* nl_part = malloc(params->lowmc->r * sizeof(uint32_t));
    shares_t* nl_part_masks = allocateShares(params->lowmc->r * 32);
    shares_t* key_masks = allocateShares(mask_shares->numWords);    // Make a copy to use when computing each round key
    shares_t* mask2_shares = allocateShares(mask_shares->numWords);
    uint8_t* unopened_msgs = NULL;

    if(msgs->unopened >= 0) { //We are in verify, save the unopenend parties msgs
        unopened_msgs = malloc(params->view_size + params->input_size);
        memcpy(unopened_msgs, msgs->msgs[msgs->unopened], params->view_size + params->input_size);
    }

    copyShares(key_masks, mask_shares);

#if defined(REDUCED_ROUND_KEY_COMPUTATION)
    mpc_matrix_mul(state, maskedKey, params->lowmc->k0_matrix->w64, mask_shares, params);       // roundKey = maskedKey * KMatrix[0]
    xor_array(state, state, plaintext, (params->input_size / 4));                                // state = plaintext + roundKey
    xor_array_RC((uint8_t*)state, (uint8_t*)state, (uint8_t*)params->lowmc->precomputed_constant_linear, params->input_size);  // state = state + precomp_const
    mpc_matrix_mul_nl_part(nl_part, maskedKey, params->lowmc->precomputed_non_linear_part_matrix->w64,
            params->lowmc->precomputed_constant_non_linear->w64, nl_part_masks, key_masks, params);
#if defined(OPTIMIZED_LINEAR_LAYER_EVALUATION)
    for (uint32_t r = 0; r < params->lowmc->r-1; r++) {
        mpc_sbox(state, mask_shares, tapes, msgs, unopened_msgs, params);
        mpc_xor2_nl(state, mask_shares, state, mask_shares, nl_part, nl_part_masks, r*32+2, 30);    // state += roundKey
        //mpc_matrix_mul(state, state, params->lowmc->rounds[r].l_matrix->w64, mask_shares, params);              // state = state * LMatrix (r-1)
        mpc_matrix_mul_z(state2, state, mask2_shares, mask_shares, params->lowmc->rounds[r].z_matrix->w64, params);
        mpc_shuffle((uint8_t*)state, mask_shares, params->lowmc->rounds[r].r_mask);
        mpc_matrix_addmul_r(state2, state, mask2_shares, mask_shares, params->lowmc->rounds[r].r_matrix->w64, params);
        for(uint32_t i = 0; i < params->lowmc->m*3; i++) {
            mask_shares->shares[i] = 0;
            setBit((uint8_t*)state, i, 0);
        }
        mpc_xor2(state, mask_shares, state, mask_shares, state2, mask2_shares, params);
    }
    mpc_sbox(state, mask_shares, tapes, msgs, unopened_msgs, params);
    mpc_xor2_nl(state, mask_shares, state, mask_shares, nl_part, nl_part_masks, (params->lowmc->r-1)*32+2, 30);    // state += roundKey
    mpc_matrix_mul(state, state, params->lowmc->zr_matrix->w64, mask_shares, params);              // state = state * LMatrix (r-1)
#else
    for (uint32_t r = 0; r < params->lowmc->r; r++) {
        mpc_sbox(state, mask_shares, tapes, msgs, params);
        mpc_xor2_nl(state, mask_shares, state, mask_shares, nl_part, nl_part_masks, r*32+2, 30, params);    // state += roundKey
        mpc_matrix_mul(state, state, params->lowmc->rounds[r].l_matrix->w64, mask_shares, params);              // state = state * LMatrix (r-1)
    }
#endif
#else
    mpc_matrix_mul(roundKey, maskedKey, params->lowmc->k0_matrix->w64, mask_shares, params);       // roundKey = maskedKey * KMatrix[0]
    xor_array(state, roundKey, plaintext, (params->input_size / 4));                      // state = plaintext + roundKey

    shares_t* round_key_masks = allocateShares(mask_shares->numWords);
    for (uint32_t r = 0; r < params->lowmc->r; r++) {
        copyShares(round_key_masks, key_masks);
        mpc_matrix_mul(roundKey, maskedKey, params->lowmc->rounds[r].k_matrix->w64, round_key_masks, params);

        mpc_sbox(state, mask_shares, tapes, msgs, params);
        mpc_matrix_mul(state, state, params->lowmc->rounds[r].l_matrix->w64, mask_shares, params);              // state = state * LMatrix (r-1)
        xor_array_RC(state, state, (const uint8_t*)params->lowmc->rounds[r].constant->w64, params->input_size);              // state += RConstant
        mpc_xor2(state, mask_shares, roundKey, round_key_masks, state, mask_shares, params);    // state += roundKey
    }
    freeShares(round_key_masks);
#endif

    /* Unmask the output, and check that it's correct */
    if (msgs->unopened >= 0) {
        /* During signature verification we have the shares of the output for
         * the unopened party already in msgs, but not in mask_shares. */
        for (size_t i = 0; i < params->lowmc->n; i++) {
            uint8_t share = getBit(unopened_msgs, msgs->pos + i);
            setBit((uint8_t*)&mask_shares->shares[i],  msgs->unopened, share);
        }

    }
    uint32_t output[LOWMC_MAX_STATE_SIZE];
    reconstructShares(output, mask_shares);
    xor_array(output, output, state, (params->input_size / 4));

    if (memcmp(output, pubKey, params->input_size) != 0) {
        printf("%s: output does not match pubKey\n", __func__);
        printHex("pubKey", (uint8_t*)pubKey, params->input_size);
        printHex("output", (uint8_t*)output, params->input_size);
        ret = -1;
        goto Exit;
    }

    broadcast(mask_shares, msgs);
    msgsTranspose(msgs);

    free(unopened_msgs);
    free(state);
    free(state2);
    free(roundKey);
    free(nl_part);
    freeShares(key_masks);
    freeShares(mask2_shares);
    freeShares(nl_part_masks);

Exit:
    return ret;
}

static size_t bitsToChunks(size_t chunkLenBits, const uint8_t* input, size_t inputLen, uint16_t* chunks)
{
    if (chunkLenBits > inputLen * 8) {
        assert(!"Invalid input to bitsToChunks: not enough input");
        return 0;
    }
    size_t chunkCount = ((inputLen * 8) / chunkLenBits);

    for (size_t i = 0; i < chunkCount; i++) {
        chunks[i] = 0;
        for (size_t j = 0; j < chunkLenBits; j++) {
            chunks[i] += getBit(input, i * chunkLenBits + j) << j;
            assert(chunks[i] < (1 << chunkLenBits));
        }
        chunks[i] = le16toh(chunks[i]);
    }

    return chunkCount;
}

static size_t appendUnique(uint16_t* list, uint16_t value, size_t position)
{
    if (position == 0) {
        list[position] = value;
        return position + 1;
    }

    for (size_t i = 0; i < position; i++) {
        if (list[i] == value) {
            return position;
        }
    }
    list[position] = value;
    return position + 1;
}

static void HCP(uint16_t* challengeC, uint16_t* challengeP, commitments_t* Ch,
                uint8_t* hCv, uint8_t* salt, const uint32_t* pubKey, const uint32_t* plaintext, const uint8_t* message,
                size_t messageByteLength, const picnic_instance_t* params)
{
    Keccak_HashInstance ctx;
    uint8_t h[MAX_DIGEST_SIZE] = { 0 };

    assert(params->num_opened_rounds < params->num_rounds);

#if 0  // Print out inputs when debugging
    printf("\n");
    for (size_t t = 0; t < params->num_rounds; t++) {
        printf("%s Ch[%lu]", __func__, t);
        printHex("", Ch->hashes[t], params->digest_size);

    }
    printHex("hCv", hCv, params->digest_size);

    printf("%s salt", __func__);
    printHex("", salt, params->seed_size);
    printf("%s pubKey", __func__);
    printHex("", (uint8_t*)pubKey, params->input_size);
    printf("%s plaintext", __func__);
    printHex("", (uint8_t*)plaintext, params->input_size);

#endif

    hash_init(&ctx, params);
    for (size_t t = 0; t < params->num_rounds; t++) {
        hash_update(&ctx, Ch->hashes[t], params->digest_size);
    }

    hash_update(&ctx, hCv, params->digest_size);
    hash_update(&ctx, salt, params->seed_size);
    hash_update(&ctx, (uint8_t*)pubKey, params->input_size);
    hash_update(&ctx, (uint8_t*)plaintext, params->input_size);
    hash_update(&ctx, message, messageByteLength);
    hash_final(&ctx);
    hash_squeeze(&ctx, h, params->digest_size);

    // Populate C
    uint32_t bitsPerChunkC = ceil_log2(params->num_rounds);
    uint32_t bitsPerChunkP = ceil_log2(params->num_MPC_parties);
    uint16_t* chunks = calloc(params->digest_size * 8 / bitsPerChunkP, sizeof(uint16_t));

    size_t countC = 0;
    while (countC < params->num_opened_rounds) {
        size_t numChunks = bitsToChunks(bitsPerChunkC, h, params->digest_size, chunks);
        for (size_t i = 0; i < numChunks; i++) {
            if (chunks[i] < params->num_rounds) {
                countC = appendUnique(challengeC, chunks[i], countC);
            }
            if (countC == params->num_opened_rounds) {
                break;
            }
        }

        hash_init_prefix(&ctx, params, HASH_PREFIX_1);
        hash_update(&ctx, h, params->digest_size);
        hash_final(&ctx);
        hash_squeeze(&ctx, h, params->digest_size);
    }

    // Note that we always compute h = H(h) after setting C
    size_t countP = 0;

    while (countP < params->num_opened_rounds) {
        size_t numChunks = bitsToChunks(bitsPerChunkP, h, params->digest_size, chunks);
        for (size_t i = 0; i < numChunks; i++) {
            if (chunks[i] < params->num_MPC_parties) {
                challengeP[countP] = chunks[i];
                countP++;
            }
            if (countP == params->num_opened_rounds) {
                break;
            }
        }

        hash_init_prefix(&ctx, params, HASH_PREFIX_1);
        hash_update(&ctx, h, params->digest_size);
        hash_final(&ctx);
        hash_squeeze(&ctx, h, params->digest_size);
    }

#if 0   // Print challenge when debugging
    printf("C = ");
    for (size_t i = 0; i < countC; i++) {
        printf("%u, ", challengeC[i]);
    }
    printf("\n");

    printf("P = ");
    for (size_t i = 0; i < countP; i++) {
        printf("%u, ", challengeP[i]);
    }
    printf("\n");
#endif

    free(chunks);

}

static uint16_t* getMissingLeavesList(uint16_t* challengeC, const picnic_instance_t* params)
{
    size_t missingLeavesSize = params->num_rounds - params->num_opened_rounds;
    uint16_t* missingLeaves = calloc(missingLeavesSize, sizeof(uint16_t));
    size_t pos = 0;

    for (size_t i = 0; i < params->num_rounds; i++) {
        if (!contains(challengeC, params->num_opened_rounds, i)) {
            missingLeaves[pos] = i;
            pos++;
        }
    }

    return missingLeaves;
}

int verify_picnic2(signature2_t* sig, const uint32_t* pubKey, const uint32_t* plaintext, const uint8_t* message, size_t messageByteLength,
                   const picnic_instance_t* params)
{
    commitments_t* C = allocateCommitments(params, 0);
    commitments_t Ch = { 0 };
    commitments_t Cv = { 0 };
    msgs_t* msgs = allocateMsgs(params);
    tree_t* treeCv = createTree(params->num_rounds, params->digest_size);
    size_t challengeSizeBytes = params->num_opened_rounds * sizeof(uint16_t);
    uint16_t* challengeC = malloc(challengeSizeBytes);
    uint16_t* challengeP = malloc(challengeSizeBytes);
    tree_t** seeds = calloc(params->num_rounds, sizeof(tree_t*));
    randomTape_t* tapes = malloc(params->num_rounds * sizeof(randomTape_t));
    tree_t* iSeedsTree = createTree(params->num_rounds, params->seed_size);
    int ret = reconstructSeeds(iSeedsTree, sig->challengeC, params->num_opened_rounds, sig->iSeedInfo, sig->iSeedInfoLen, sig->salt, 0, params);

    if (ret != 0) {
        ret = -1;
        goto Exit;
    }

    /* Populate seeds with values from the signature */
    for (size_t t = 0; t < params->num_rounds; t++) {
        if (!contains(sig->challengeC, params->num_opened_rounds, t)) {
            /* Expand iSeed[t] to seeds for each parties, using a seed tree */
            seeds[t] = generateSeeds(params->num_MPC_parties, getLeaf(iSeedsTree, t), sig->salt, t, params);
        }
        else {
            /* We don't have the initial seed for the round, but instead a seed
             * for each unopened party */
            seeds[t] = createTree(params->num_MPC_parties, params->seed_size);
            size_t P_index = indexOf(sig->challengeC, params->num_opened_rounds, t);
            uint16_t hideList[1];
            hideList[0] = sig->challengeP[P_index];
            ret = reconstructSeeds(seeds[t], hideList, 1,
                                   sig->proofs[t].seedInfo, sig->proofs[t].seedInfoLen,
                                   sig->salt, t, params);
            if (ret != 0) {
                printf("Failed to reconstruct seeds for round %lu\n", t);
                ret = -1;
                goto Exit;
            }
        }
    }

    /* Commit */
    size_t last = params->num_MPC_parties - 1;
    uint8_t auxBits[MAX_AUX_BYTES];
    for (size_t t = 0; t < params->num_rounds; t++) {
        /* Compute random tapes for all parties.  One party for each repitition
         * challengeC will have a bogus seed; but we won't use that party's
         * random tape. */
        createRandomTapes(&tapes[t], getLeaves(seeds[t]), sig->salt, t, params);

        if (!contains(sig->challengeC, params->num_opened_rounds, t)) {
            /* We're given iSeed, have expanded the seeds, compute aux from scratch so we can comnpte Com[t] */
            computeAuxTape(&tapes[t], params);
            for (size_t j = 0; j < last; j++) {
                commit(C[t].hashes[j], getLeaf(seeds[t], j), NULL, sig->salt, t, j, params);
            }
            getAuxBits(auxBits, &tapes[t], params);
            commit(C[t].hashes[last], getLeaf(seeds[t], last), auxBits, sig->salt, t, last, params);
        }
        else {
            /* We're given all seeds and aux bits, execpt for the unopened 
             * party, we get their commitment */
            size_t unopened = sig->challengeP[indexOf(sig->challengeC, params->num_opened_rounds, t)];
            for (size_t j = 0; j < last; j++) {
                if (j != unopened) {
                    commit(C[t].hashes[j], getLeaf(seeds[t], j), NULL, sig->salt, t, j, params);
                }
            }
            if (last != unopened) {
                commit(C[t].hashes[last], getLeaf(seeds[t], last), sig->proofs[t].aux, sig->salt, t, last, params);
            }

            memcpy(C[t].hashes[unopened], sig->proofs[t].C, params->digest_size);
        }

    }


    /* Commit to the commitments */
    allocateCommitments2(&Ch, params, params->num_rounds);
    for (size_t t = 0; t < params->num_rounds; t++) {
        commit_h(Ch.hashes[t], &C[t], params);
    }

    /* Commit to the views */
    allocateCommitments2(&Cv, params, params->num_rounds);
    shares_t* mask_shares = allocateShares(params->lowmc->n);
    for (size_t t = 0; t < params->num_rounds; t++) {
        if (contains(sig->challengeC, params->num_opened_rounds, t)) {
            /* 2. When t is in C, we have everything we need to re-compute the view, as an honest signer would.
             * We simulate the MPC with one fewer party; the unopned party's values are all set to zero. */
            size_t unopened = sig->challengeP[indexOf(sig->challengeC, params->num_opened_rounds, t)];
            size_t tapeLengthBytes = 2 * params->view_size + params->input_size;
            setAuxBits(&tapes[t], sig->proofs[t].aux, params);
            memset(tapes[t].tape[unopened], 0, tapeLengthBytes);
            memcpy(msgs[t].msgs[unopened], sig->proofs[t].msgs, params->view_size + params->input_size );
            msgs[t].unopened = unopened;

            tapesToWords(mask_shares, &tapes[t]);
            int ret = simulateOnline((uint32_t*)sig->proofs[t].input, mask_shares, &tapes[t], &msgs[t], plaintext, pubKey, params);
            if (ret != 0) {
                printf("MPC simulation failed for round %lu, signature invalid\n", t);
                ret = -1;
                freeShares(mask_shares);
                goto Exit;
            }
            commit_v(Cv.hashes[t], sig->proofs[t].input, &msgs[t], params);
        }
        else {
            Cv.hashes[t] = NULL;
        }
    }
    freeShares(mask_shares);

    size_t missingLeavesSize = params->num_rounds - params->num_opened_rounds;
    uint16_t* missingLeaves = getMissingLeavesList(sig->challengeC, params);
    ret = addMerkleNodes(treeCv, missingLeaves, missingLeavesSize, sig->cvInfo, sig->cvInfoLen);
    free(missingLeaves);
    if (ret != 0) {
        ret = -1;
        goto Exit;
    }

    ret = verifyMerkleTree(treeCv, Cv.hashes, sig->salt, params);
    if (ret != 0) {
        ret = -1;
        goto Exit;
    }

    /* Compute the challenge; two lists of integers */
    HCP(challengeC, challengeP, &Ch, treeCv->nodes[0], sig->salt, pubKey, plaintext, message, messageByteLength, params);

    /* Compare to challenge from signature */
    if ( memcmp(sig->challengeC, challengeC, challengeSizeBytes) != 0 ||
         memcmp(sig->challengeP, challengeP, challengeSizeBytes) != 0 ) {
        printf("Challenge does not match, signature invalid\n");
        ret = -1;
        goto Exit;
    }

    ret = EXIT_SUCCESS;

Exit:

    free(challengeC);
    free(challengeP);
    freeCommitments(C);
    freeCommitments2(&Cv);
    freeCommitments2(&Ch);
    freeMsgs(msgs);
    freeTree(treeCv);
    freeTree(iSeedsTree);
    for (size_t t = 0; t < params->num_rounds; t++) {
        freeRandomTape(&tapes[t]);
        freeTree(seeds[t]);
    }
    free(seeds);
    free(tapes);

    return ret;
}

static void computeSaltAndRootSeed(uint8_t* saltAndRoot, size_t saltAndRootLength, uint32_t* privateKey, uint32_t* pubKey,
                                   uint32_t* plaintext, const uint8_t* message, size_t messageByteLength, const picnic_instance_t* params)
{
    Keccak_HashInstance ctx;

    hash_init(&ctx, params);
    hash_update(&ctx, (uint8_t*)privateKey, params->input_size);
    hash_update(&ctx, message, messageByteLength);
    hash_update(&ctx, (uint8_t*)pubKey, params->input_size);
    hash_update(&ctx, (uint8_t*)plaintext, params->input_size);
    uint16_t stateSizeLE = htole16((uint16_t) params->lowmc->n);
    hash_update(&ctx, (uint8_t*)&stateSizeLE, sizeof(uint16_t));
    hash_final(&ctx);
    hash_squeeze(&ctx, saltAndRoot, saltAndRootLength);
}

int sign_picnic2(uint32_t* privateKey, uint32_t* pubKey, uint32_t* plaintext, const uint8_t* message,
                 size_t messageByteLength, signature2_t* sig, const picnic_instance_t* params)
{
    int ret = 0;
    uint8_t* saltAndRoot = malloc(params->seed_size * 2);

    computeSaltAndRootSeed(saltAndRoot, params->seed_size * 2, privateKey, pubKey, plaintext, message, messageByteLength, params);
    memcpy(sig->salt, saltAndRoot, params->seed_size);
    tree_t* iSeedsTree = generateSeeds(params->num_rounds, saltAndRoot + params->seed_size, sig->salt, 0, params);
    uint8_t** iSeeds = getLeaves(iSeedsTree);
    free(saltAndRoot);

    randomTape_t* tapes = malloc(params->num_rounds * sizeof(randomTape_t));
    tree_t** seeds = malloc(params->num_rounds * sizeof(tree_t*));
    for (size_t t = 0; t < params->num_rounds; t++) {
        seeds[t] = generateSeeds(params->num_MPC_parties, iSeeds[t], sig->salt, t, params);
        createRandomTapes(&tapes[t], getLeaves(seeds[t]), sig->salt, t, params);
    }

    /* Preprocessing; compute aux tape for the N-th player, for each parallel rep */
    uint8_t auxBits[MAX_AUX_BYTES];
    for (size_t t = 0; t < params->num_rounds; t++) {
        computeAuxTape(&tapes[t], params);
    }

    /* Commit to seeds and aux bits */
    commitments_t* C = allocateCommitments(params, 0);
    for (size_t t = 0; t < params->num_rounds; t++) {
        for (size_t j = 0; j < params->num_MPC_parties - 1; j++) {
            commit(C[t].hashes[j], getLeaf(seeds[t], j), NULL, sig->salt, t, j, params);
        }
        size_t last = params->num_MPC_parties - 1;
        getAuxBits(auxBits, &tapes[t], params);
        commit(C[t].hashes[last], getLeaf(seeds[t], last), auxBits, sig->salt, t, last, params);
    }

    /* Simulate the online phase of the MPC */
    inputs_t inputs = allocateInputs(params);
    msgs_t* msgs = allocateMsgs(params);
    shares_t* mask_shares = allocateShares(params->lowmc->n);
    for (size_t t = 0; t < params->num_rounds; t++) {
        uint32_t* maskedKey = (uint32_t*)inputs[t];

        tapesToWords(mask_shares, &tapes[t]);
        reconstructShares(maskedKey, mask_shares);                                      // maskedKey = masks
        xor_array(maskedKey, maskedKey, privateKey, (params->input_size / 4));            // maskedKey += privateKey

        int rv = simulateOnline(maskedKey, mask_shares, &tapes[t], &msgs[t], plaintext, pubKey, params);
        if (rv != 0) {
            printf("MPC simulation failed, aborting signature\n");
            ret = -1;
        }
    }
    freeShares(mask_shares);

    /* Commit to the commitments and views */
    commitments_t Ch;
    allocateCommitments2(&Ch, params, params->num_rounds);
    commitments_t Cv;
    allocateCommitments2(&Cv, params, params->num_rounds);
    for (size_t t = 0; t < params->num_rounds; t++) {
        commit_h(Ch.hashes[t], &C[t], params);
        commit_v(Cv.hashes[t], inputs[t], &msgs[t], params);
    }

    /* Create a Merkle tree with Cv as the leaves */
    tree_t* treeCv = createTree(params->num_rounds, params->digest_size);
    buildMerkleTree(treeCv, Cv.hashes, sig->salt, params);

    /* Compute the challenge; two lists of integers */
    uint16_t* challengeC = sig->challengeC;
    uint16_t* challengeP = sig->challengeP;
    HCP(challengeC, challengeP, &Ch, treeCv->nodes[0], sig->salt, pubKey, plaintext, message, messageByteLength, params);

    /* Send information required for checking commitments with Merkle tree.
     * The commitments the verifier will be missing are those not in challengeC. */
    size_t missingLeavesSize = params->num_rounds - params->num_opened_rounds;
    uint16_t* missingLeaves = getMissingLeavesList(challengeC, params);
    size_t cvInfoLen = 0;
    uint8_t* cvInfo = openMerkleTree(treeCv, missingLeaves, missingLeavesSize, &cvInfoLen);
    sig->cvInfo = cvInfo;
    sig->cvInfoLen = cvInfoLen;
    free(missingLeaves);

    /* Reveal iSeeds for unopned rounds, those in {0..T-1} \ ChallengeC. */
    sig->iSeedInfo = malloc(params->num_rounds * params->seed_size);
    sig->iSeedInfoLen = revealSeeds(iSeedsTree, challengeC, params->num_opened_rounds,
                                    sig->iSeedInfo, params->num_rounds * params->seed_size, params);
    sig->iSeedInfo = realloc(sig->iSeedInfo, sig->iSeedInfoLen);

    /* Assemble the proof */
    proof2_t* proofs = sig->proofs;
    for (size_t t = 0; t < params->num_rounds; t++) {
        if (contains(challengeC, params->num_opened_rounds, t)) {
            allocateProof2(&proofs[t], params);
            size_t P_index = indexOf(challengeC, params->num_opened_rounds, t);
            proofs[t].unOpenedIndex = challengeP[P_index];

            uint16_t hideList[1];
            hideList[0] = challengeP[P_index];
            proofs[t].seedInfo = malloc(params->num_MPC_parties * params->seed_size);
            proofs[t].seedInfoLen = revealSeeds(seeds[t], hideList, 1, proofs[t].seedInfo, params->num_MPC_parties * params->seed_size, params);
            proofs[t].seedInfo = realloc(proofs[t].seedInfo, proofs[t].seedInfoLen);

            size_t last = params->num_MPC_parties - 1;
            if (challengeP[P_index] != last) {
                getAuxBits(proofs[t].aux, &tapes[t], params);
            }

            memcpy(proofs[t].input, inputs[t], params->input_size);
            memcpy(proofs[t].msgs, msgs[t].msgs[challengeP[P_index]], params->view_size + params->input_size );
            memcpy(proofs[t].C, C[t].hashes[proofs[t].unOpenedIndex], params->digest_size);
        }
    }

    sig->proofs = proofs;

#if 0
    printf("\n-----------------\n\nSelf-Test, trying to verify signature:\n");
    int ret = verify_picnic2(sig, pubKey, plaintext, message, messageByteLength, params);
    if (ret != 0) {
        printf("Verification failed; signature invalid\n");
        ret = -1;
    }
    else {
        printf("Verification succeeded\n\n");
    }
    printf("-----------------\n\nSelf-Test complete\n");

#endif

    for (size_t t = 0; t < params->num_rounds; t++) {
        freeRandomTape(&tapes[t]);
        freeTree(seeds[t]);
    }
    free(tapes);
    free(seeds);
    freeTree(iSeedsTree);
    freeTree(treeCv);

    freeCommitments(C);
    freeCommitments2(&Ch);
    freeCommitments2(&Cv);
    freeInputs(inputs);
    freeMsgs(msgs);

    return ret;

}


static int inRange(uint16_t* list, size_t len, size_t low, size_t high)
{
    for (size_t i = 0; i < len; i++) {
        if (list[i] > high || list[i] < low) {
            return 0;
        }
    }
    return 1;
}

static int unique(uint16_t* list, size_t len)
{
    for (size_t i = 0; i < len; i++) {
        for (size_t j = 0; j < len; j++) {
            if (j != i && list[i] == list[j]) {
                return 0;
            }
        }
    }
    return 1;
}

static int arePaddingBitsZero(uint8_t* data, size_t byteLength, size_t bitLength)
{
    for (size_t i = bitLength; i < byteLength * 8; i++) {
        uint8_t bit_i = getBit(data, i);
        if (bit_i != 0) {
            return 0;
        }
    }
    return 1;
}

int deserializeSignature2(signature2_t* sig, const uint8_t* sigBytes, size_t sigBytesLen, const picnic_instance_t* params)
{
    /* Read the challenge and salt */
    size_t bytesRequired = 4 * params->num_opened_rounds + params->seed_size;

    if (sigBytesLen < bytesRequired) {
        return EXIT_FAILURE;
    }

    memcpy(sig->challengeC, sigBytes, 2 * params->num_opened_rounds);
    sigBytes += 2 * params->num_opened_rounds;
    memcpy(sig->challengeP, sigBytes, 2 * params->num_opened_rounds);
    sigBytes += 2 * params->num_opened_rounds;
    memcpy(sig->salt, sigBytes, params->seed_size);
    sigBytes += params->seed_size;

    for (size_t i = 0; i < params->num_opened_rounds; i++) {
        sig->challengeC[i] = le16toh(sig->challengeC[i]);
        sig->challengeP[i] = le16toh(sig->challengeP[i]);
    }

    if (!inRange(sig->challengeC, params->num_opened_rounds, 0, params->num_rounds - 1)) {
        return EXIT_FAILURE;
    }
    if (!unique(sig->challengeC, params->num_opened_rounds)) {
        return EXIT_FAILURE;
    }
    if (!inRange(sig->challengeP, params->num_opened_rounds, 0, params->num_MPC_parties - 1)) {
        return EXIT_FAILURE;
    }

    /* Add size of iSeeds tree data */
    sig->iSeedInfoLen = revealSeedsSize(params->num_rounds, sig->challengeC, params->num_opened_rounds, params);
    bytesRequired += sig->iSeedInfoLen;

    /* Add the size of the Cv Merkle tree data */
    size_t missingLeavesSize = params->num_rounds - params->num_opened_rounds;
    uint16_t* missingLeaves = getMissingLeavesList(sig->challengeC, params);
    sig->cvInfoLen = openMerkleTreeSize(params->num_rounds, missingLeaves, missingLeavesSize, params);
    bytesRequired += sig->cvInfoLen;
    free(missingLeaves);

    /* Compute the number of bytes required for the proofs */
    uint16_t hideList[1] = { 0 };
    size_t seedInfoLen = revealSeedsSize(params->num_MPC_parties, hideList, 1, params);
    for (size_t t = 0; t < params->num_rounds; t++) {
        if (contains(sig->challengeC, params->num_opened_rounds, t)) {
            size_t P_t = sig->challengeP[indexOf(sig->challengeC, params->num_opened_rounds, t)];
            if (P_t != (params->num_MPC_parties - 1)) {
                bytesRequired += params->view_size;
            }
            bytesRequired += params->digest_size;
            bytesRequired += params->input_size;
            bytesRequired += params->input_size + params->view_size;
            bytesRequired += seedInfoLen;
        }
    }

    /* Fail if the signature does not have the exact number of bytes we expect */
    if (sigBytesLen != bytesRequired) {
        printf("%s: sigBytesLen = %lu, expected bytesRequired = %lu\n", __func__, sigBytesLen, bytesRequired);
        return EXIT_FAILURE;
    }

    sig->iSeedInfo = malloc(sig->iSeedInfoLen);
    memcpy(sig->iSeedInfo, sigBytes, sig->iSeedInfoLen);
    sigBytes += sig->iSeedInfoLen;

    sig->cvInfo = malloc(sig->cvInfoLen);
    memcpy(sig->cvInfo, sigBytes, sig->cvInfoLen);
    sigBytes += sig->cvInfoLen;

    /* Read the proofs */
    for (size_t t = 0; t < params->num_rounds; t++) {
        if (contains(sig->challengeC, params->num_opened_rounds, t)) {
            allocateProof2(&sig->proofs[t], params);
            sig->proofs[t].seedInfoLen = seedInfoLen;
            sig->proofs[t].seedInfo = malloc(sig->proofs[t].seedInfoLen);
            memcpy(sig->proofs[t].seedInfo, sigBytes, sig->proofs[t].seedInfoLen);
            sigBytes += sig->proofs[t].seedInfoLen;

            size_t P_t = sig->challengeP[indexOf(sig->challengeC, params->num_opened_rounds, t)];
            if (P_t != (params->num_MPC_parties - 1) ) {
                memcpy(sig->proofs[t].aux, sigBytes, params->view_size);
                sigBytes += params->view_size;
                if (!arePaddingBitsZero(sig->proofs[t].aux, params->view_size, 3 * params->lowmc->r * params->lowmc->m)) {
                    printf("%s: failed while deserializing aux bits\n", __func__);
                    return -1;
                }
            }

            memcpy(sig->proofs[t].input, sigBytes, params->seed_size);
            sigBytes += params->input_size;

            size_t msgsByteLength = params->input_size + params->view_size;
            memcpy(sig->proofs[t].msgs, sigBytes, msgsByteLength);
            sigBytes += msgsByteLength;
            size_t msgsBitLength = params->lowmc->n + 3 * params->lowmc->r * params->lowmc->m;
            if (!arePaddingBitsZero(sig->proofs[t].msgs, msgsByteLength, msgsBitLength)) {
                printf("%s: failed while deserializing msgs bits\n", __func__);
                return -1;
            }

            memcpy(sig->proofs[t].C, sigBytes, params->digest_size);
            sigBytes += params->digest_size;
        }
    }

    return EXIT_SUCCESS;
}

int serializeSignature2(const signature2_t* sig, uint8_t* sigBytes, size_t sigBytesLen, const picnic_instance_t* params)
{
    uint8_t* sigBytesBase = sigBytes;

    /* Compute the number of bytes required for the signature */
    size_t bytesRequired = 4 * params->num_opened_rounds + params->seed_size; /* challenge and salt */

    bytesRequired += sig->iSeedInfoLen;                                         /* Encode only iSeedInfo, the length will be recomputed by deserialize */
    bytesRequired += sig->cvInfoLen;

    for (size_t t = 0; t < params->num_rounds; t++) {   /* proofs */
        if (contains(sig->challengeC, params->num_opened_rounds, t)) {
            size_t P_t = sig->challengeP[indexOf(sig->challengeC, params->num_opened_rounds, t)];
            bytesRequired += sig->proofs[t].seedInfoLen;
            if (P_t != (params->num_MPC_parties - 1)) {
                bytesRequired += params->view_size;
            }
            bytesRequired += params->digest_size;
            bytesRequired += params->input_size;
            bytesRequired += params->input_size + params->view_size;
        }
    }

    if (sigBytesLen < bytesRequired) {
        return -1;
    }

    memcpy(sigBytes, sig->challengeC, 2 * params->num_opened_rounds);
    uint16_t* challengeC = (uint16_t*)sigBytes;
    sigBytes += 2 * params->num_opened_rounds;
    memcpy(sigBytes, sig->challengeP, 2 * params->num_opened_rounds);
    uint16_t* challengeP = (uint16_t*)sigBytes;
    sigBytes += 2 * params->num_opened_rounds;
    memcpy(sigBytes, sig->salt, params->seed_size);
    sigBytes += params->seed_size;

    for (size_t i = 0; i < params->num_opened_rounds; i++) {
        challengeC[i] = le16toh(sig->challengeC[i]);
        challengeP[i] = le16toh(sig->challengeP[i]);
    }

    memcpy(sigBytes, sig->iSeedInfo, sig->iSeedInfoLen);
    sigBytes += sig->iSeedInfoLen;
    memcpy(sigBytes, sig->cvInfo, sig->cvInfoLen);
    sigBytes += sig->cvInfoLen;

    /* Write the proofs */
    for (size_t t = 0; t < params->num_rounds; t++) {
        if (contains(sig->challengeC, params->num_opened_rounds, t)) {
            memcpy(sigBytes, sig->proofs[t].seedInfo,  sig->proofs[t].seedInfoLen);
            sigBytes += sig->proofs[t].seedInfoLen;

            size_t P_t = sig->challengeP[indexOf(sig->challengeC, params->num_opened_rounds, t)];

            if (P_t != (params->num_MPC_parties - 1) ) {
                memcpy(sigBytes, sig->proofs[t].aux, params->view_size);
                sigBytes += params->view_size;
            }

            memcpy(sigBytes, sig->proofs[t].input, params->seed_size);
            sigBytes += params->input_size;

            memcpy(sigBytes, sig->proofs[t].msgs, params->input_size + params->view_size);
            sigBytes += params->input_size + params->view_size;

            memcpy(sigBytes, sig->proofs[t].C, params->digest_size);
            sigBytes += params->digest_size;
        }
    }

    return (int)(sigBytes - sigBytesBase);
}

int impl_sign_picnic2(const picnic_instance_t* instance, const uint8_t* plaintext, const uint8_t* private_key,
                      const uint8_t* public_key, const uint8_t* msg, size_t msglen, uint8_t* signature,
                      size_t* signature_len) {
    int ret;
    signature2_t* sig = (signature2_t*)malloc(sizeof(signature2_t));
    allocateSignature2(sig, instance);
    if (sig == NULL) {
        return -1;
    }
    ret = sign_picnic2((uint32_t*)private_key, (uint32_t*)public_key, (uint32_t*)plaintext, msg,
                       msglen, sig, instance);
    if (ret != EXIT_SUCCESS) {
        fprintf(stderr, "Failed to create signature\n");
        fflush(stderr);
        freeSignature2(sig, instance);
        free(sig);
        return -1;
    }
    ret = serializeSignature2(sig, signature, *signature_len, instance);
    if (ret == -1) {
        fprintf(stderr, "Failed to serialize signature\n");
        fflush(stderr);
        freeSignature2(sig, instance);
        free(sig);
        return -1;
    }
    *signature_len = ret;

    freeSignature2(sig, instance);
    free(sig);
    return 0;
}

int impl_verify_picnic2(const picnic_instance_t* instance, const uint8_t* plaintext, const uint8_t* public_key,
                const uint8_t* msg, size_t msglen, const uint8_t* signature, size_t signature_len) {
    int ret;
    signature2_t* sig = (signature2_t*)malloc(sizeof(signature2_t));
    allocateSignature2(sig, instance);
    if (sig == NULL) {
        return -1;
    }

    ret = deserializeSignature2(sig, signature, signature_len, instance);
    if (ret != EXIT_SUCCESS) {
        fprintf(stderr, "Failed to deserialize signature\n");
        fflush(stderr);
        freeSignature2(sig, instance);
        free(sig);
        return -1;
    }

    ret = verify_picnic2(sig, (uint32_t*)public_key,
                         (uint32_t*)plaintext, msg, msglen, instance);
    if (ret != EXIT_SUCCESS) {
        /* Signature is invalid, or verify function failed */
        freeSignature2(sig, instance);
        free(sig);
        return -1;
    }

    freeSignature2(sig, instance);
    free(sig);
    return 0;
}
