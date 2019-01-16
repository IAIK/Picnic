/*
Implementation by the Keccak Team, namely, Guido Bertoni, Joan Daemen,
Michaël Peeters, Gilles Van Assche and Ronny Van Keer,
hereby denoted as "the implementer".

For more information, feedback or questions, please refer to our website:
https://keccak.team/

To the extent possible under law, the implementer has waived all copyright
and related or neighboring rights to the source code in this file.
http://creativecommons.org/publicdomain/zero/1.0/
*/

#include <string.h>
#include "KeccakHashtimes4.h"

/* ---------------------------------------------------------------- */

HashReturn Keccak_HashInitializetimes4(Keccak_HashInstancetimes4 *instance, unsigned int rate, unsigned int capacity, unsigned int hashbitlen, unsigned char delimitedSuffix)
{
    HashReturn result;

    if (delimitedSuffix == 0)
        return FAIL;
    result = (HashReturn)KeccakWidth1600times4_SpongeInitialize(&instance->sponge, rate, capacity);
    if (result != SUCCESS)
        return result;
    instance->fixedOutputLength = hashbitlen;
    instance->delimitedSuffix = delimitedSuffix;
    return SUCCESS;
}

/* ---------------------------------------------------------------- */

HashReturn Keccak_HashUpdatetimes4(Keccak_HashInstancetimes4 *instance, const BitSequence **data, BitLength databitlen)
{
    if ((databitlen % 8) == 0)
        return (HashReturn)KeccakWidth1600times4_SpongeAbsorb(&instance->sponge, data, databitlen/8);
    else {
        //TODO: do we need this codepath for non-full bytes?
        return FAIL;
        /*
        HashReturn ret = (HashReturn)KeccakWidth1600times4_SpongeAbsorb(&instance->sponge, data, databitlen/8);
        if (ret == SUCCESS) {
            // The last partial byte is assumed to be aligned on the least significant bits
            unsigned char lastByte[4];
            unsigned short delimitedLastBytes[4];
            unsigned char oneByte[4][1];
            unsigned char* oneBytePtr[4];
            for(unsigned int instanceIndex = 0; instanceIndex < 4; instanceIndex++) {
                lastByte[instanceIndex] = data[instanceIndex][databitlen / 8];
                // Concatenate the last few bits provided here with those of the suffix
                delimitedLastBytes[instanceIndex] = (unsigned short) ((unsigned short) lastByte[instanceIndex] |
                                                                      ((unsigned short) instance->delimitedSuffix[instanceIndex]
                                                                              << (databitlen % 8)));
                if ((delimitedLastBytes[instanceIndex] & 0xFF00) == 0x0000) {
                    instance->delimitedSuffix[instanceIndex] = delimitedLastBytes[instanceIndex] & 0xFF;
                } else {
                    oneByte[instanceIndex][0] = delimitedLastBytes[instanceIndex] & 0xFF;
                    instance->delimitedSuffix[instanceIndex] = (delimitedLastBytes[instanceIndex] >> 8) & 0xFF;
                }
            }
            //TODO: fix with if/else from above...
            ret = (HashReturn) KeccakWidth1600times4_SpongeAbsorb(&instance->sponge, oneByte, 1);
        }
        return ret;
        */
    }
}

/* ---------------------------------------------------------------- */

HashReturn Keccak_HashFinaltimes4(Keccak_HashInstancetimes4 *instance, BitSequence **hashval)
{
    HashReturn ret = (HashReturn)KeccakWidth1600times4_SpongeAbsorbLastFewBits(&instance->sponge, instance->delimitedSuffix);
    if (ret == SUCCESS)
        return (HashReturn)KeccakWidth1600times4_SpongeSqueeze(&instance->sponge, hashval, instance->fixedOutputLength/8);
    else
        return ret;
}

/* ---------------------------------------------------------------- */

HashReturn Keccak_HashSqueezetimes4(Keccak_HashInstancetimes4 *instance, BitSequence **data, BitLength databitlen)
{
    if ((databitlen % 8) != 0)
        return FAIL;
    return (HashReturn)KeccakWidth1600times4_SpongeSqueeze(&instance->sponge, data, databitlen/8);
}
