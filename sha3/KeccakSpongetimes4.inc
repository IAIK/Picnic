/*
Implementation by the Keccak, Keyak and Ketje Teams, namely, Guido Bertoni,
Joan Daemen, Michaël Peeters, Gilles Van Assche and Ronny Van Keer, hereby
denoted as "the implementer".

For more information, feedback or questions, please refer to our websites:
http://keccak.noekeon.org/
http://keyak.noekeon.org/
http://ketje.noekeon.org/

To the extent possible under law, the implementer has waived all copyright
and related or neighboring rights to the source code in this file.
http://creativecommons.org/publicdomain/zero/1.0/
*/

#define JOIN0(a, b)                     a ## b
#define JOIN(a, b)                      JOIN0(a, b)

#define Sponge                          JOIN(prefix, _Sponge)
#define SpongeInstance                  JOIN(prefix, _SpongeInstance)
#define SpongeInitialize                JOIN(prefix, _SpongeInitialize)
#define SpongeAbsorb                    JOIN(prefix, _SpongeAbsorb)
#define SpongeAbsorbLastFewBits         JOIN(prefix, _SpongeAbsorbLastFewBits)
#define SpongeSqueeze                   JOIN(prefix, _SpongeSqueeze)

#define PlSnP_statesSizeInBytes           JOIN(PlSnP, _statesSizeInBytes)
#define PlSnP_statesAlignment             JOIN(PlSnP, _statesAlignment)
#define PlSnP_StaticInitialize            JOIN(PlSnP, _StaticInitialize)
#define PlSnP_InitializeAll               JOIN(PlSnP, _InitializeAll)
#define PlSnP_AddByte                     JOIN(PlSnP, _AddByte)
#define PlSnP_AddBytes                    JOIN(PlSnP, _AddBytes)
#define PlSnP_ExtractBytes                JOIN(PlSnP, _ExtractBytes)

int Sponge(unsigned int rate, unsigned int capacity, const unsigned char **input, size_t inputByteLen, unsigned char suffix, unsigned char *output, size_t outputByteLen)
{
//    ALIGN(PlSnP_statesAlignment) unsigned char state[PlSnP_statesSizeInBytes];
//    unsigned int partialBlock;
//    const unsigned char *curInput = input;
//    unsigned char *curOutput = output;
//    unsigned int rateInBytes = rate/8;
//
//    if (rate+capacity != PlSnP_width)
//        return 1;
//    if ((rate <= 0) || (rate > PlSnP_width) || ((rate % 8) != 0))
//        return 1;
//    if (suffix == 0)
//        return 1;
//
//    /* Initialize the state */
//    PlSnP_StaticInitialize();
//    PlSnP_InitializeAll(state);
//
//    /* First, absorb whole blocks */
//#ifdef PlSnP_FastLoop_Absorb
//    if (((rateInBytes % (PlSnP_width/200)) == 0) && (inputByteLen >= rateInBytes)) {
//        /* fast lane: whole lane rate */
//        size_t j;
//        j = PlSnP_FastLoop_Absorb(state, rateInBytes/(PlSnP_width/200), curInput, inputByteLen);
//        curInput += j;
//        inputByteLen -= j;
//    }
//#endif
//    while(inputByteLen >= (size_t)rateInBytes) {
//        #ifdef KeccakReference
//        displayBytes(1, "Block to be absorbed", curInput, rateInBytes);
//        #endif
//        PlSnP_AddBytes(state, curInput, 0, rateInBytes);
//        PlSnP_Permute(state);
//        curInput += rateInBytes;
//        inputByteLen -= rateInBytes;
//    }
//
//    /* Then, absorb what remains */
//    partialBlock = (unsigned int)inputByteLen;
//    #ifdef KeccakReference
//    displayBytes(1, "Block to be absorbed (part)", curInput, partialBlock);
//    #endif
//    PlSnP_AddBytes(state, curInput, 0, partialBlock);
//
//    /* Finally, absorb the suffix */
//    #ifdef KeccakReference
//    {
//        unsigned char delimitedData1[1];
//        delimitedData1[0] = suffix;
//        displayBytes(1, "Block to be absorbed (last few bits + first bit of padding)", delimitedData1, 1);
//    }
//    #endif
//    /* Last few bits, whose delimiter coincides with first bit of padding */
//    PlSnP_AddByte(state, suffix, partialBlock);
//    /* If the first bit of padding is at position rate-1, we need a whole new block for the second bit of padding */
//    if ((suffix >= 0x80) && (partialBlock == (rateInBytes-1)))
//        PlSnP_Permute(state);
//    /* Second bit of padding */
//    PlSnP_AddByte(state, 0x80, rateInBytes-1);
//    #ifdef KeccakReference
//    {
//        unsigned char block[PlSnP_width/8];
//        memset(block, 0, PlSnP_width/8);
//        block[rateInBytes-1] = 0x80;
//        displayBytes(1, "Second bit of padding", block, rateInBytes);
//    }
//    #endif
//    PlSnP_Permute(state);
//    #ifdef KeccakReference
//    displayText(1, "--- Switching to squeezing phase ---");
//    #endif
//
//    /* First, output whole blocks */
//    while(outputByteLen > (size_t)rateInBytes) {
//        PlSnP_ExtractBytes(state, curOutput, 0, rateInBytes);
//        PlSnP_Permute(state);
//        #ifdef KeccakReference
//        displayBytes(1, "Squeezed block", curOutput, rateInBytes);
//        #endif
//        curOutput += rateInBytes;
//        outputByteLen -= rateInBytes;
//    }
//
//    /* Finally, output what remains */
//    partialBlock = (unsigned int)outputByteLen;
//    PlSnP_ExtractBytes(state, curOutput, 0, partialBlock);
//    #ifdef KeccakReference
//    displayBytes(1, "Squeezed block (part)", curOutput, partialBlock);
//    #endif
//
    return 0;
}

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */

int SpongeInitialize(SpongeInstance *instance, unsigned int rate, unsigned int capacity)
{
    if (rate+capacity != PlSnP_width)
        return 1;
    if ((rate <= 0) || (rate > PlSnP_width) || ((rate % 8) != 0))
        return 1;
    PlSnP_StaticInitialize();
    PlSnP_InitializeAll(instance->state);
    instance->rate = rate;
    instance->byteIOIndex = 0;
    instance->squeezing = 0;

    return 0;
}

/* ---------------------------------------------------------------- */

int SpongeAbsorb(SpongeInstance *instance, const unsigned char **data, size_t dataByteLen)
{
    size_t i, j;
    unsigned int partialBlock;
    const unsigned char *curData[4];
    unsigned int rateInBytes = instance->rate/8;

    if (instance->squeezing)
        return 1; /* Too late for additional input */

    i = 0;
    if(dataByteLen > 0) {
        for (unsigned int instanceIndex = 0; instanceIndex < 4; instanceIndex++) {
            curData[instanceIndex] = data[instanceIndex];
        }
    }
    while(i < dataByteLen) {
        if ((instance->byteIOIndex == 0) && (dataByteLen >= (i + rateInBytes))) {
#ifdef PlSnP_FastLoop_Absorb
            /* processing full blocks first */
            if ((rateInBytes % (PlSnP_width/200)) == 0) {
                /* fast lane: whole lane rate */
                for(unsigned int instanceIndex = 0; instanceIndex < 4; instanceIndex++) {
                    j = PlSnP_FastLoop_Absorb(instance->state, rateInBytes/(PlSnP_width/200), curData[instanceIndex], dataByteLen - i);
                    curData[instanceIndex] += j;
                }
                i += j;
            }
            else {
#endif
                for(j=dataByteLen-i; j>=rateInBytes; j-=rateInBytes) {
                    for(unsigned int instanceIndex = 0; instanceIndex < 4; instanceIndex++) {
                        PlSnP_AddBytes(instance->state, instanceIndex, curData[instanceIndex], 0, rateInBytes);
                        curData[instanceIndex]+=rateInBytes;
                    }
                    PlSnP_Permute(instance->state);
                }
                i = dataByteLen - j;
#ifdef PlSnP_FastLoop_Absorb
            }
#endif
        }
        else {
            /* normal lane: using the message queue */
            partialBlock = (unsigned int)(dataByteLen - i);
            if (partialBlock+instance->byteIOIndex > rateInBytes)
                partialBlock = rateInBytes-instance->byteIOIndex;
            i += partialBlock;

            for(unsigned int instanceIndex = 0; instanceIndex < 4; instanceIndex++) {
                PlSnP_AddBytes(instance->state, instanceIndex, curData[instanceIndex], instance->byteIOIndex, partialBlock);
                curData[instanceIndex] += partialBlock;
            }
            instance->byteIOIndex += partialBlock;
            if (instance->byteIOIndex == rateInBytes) {
                PlSnP_Permute(instance->state);
                instance->byteIOIndex = 0;
            }
        }
    }
    return 0;
}

/* ---------------------------------------------------------------- */

int SpongeAbsorbLastFewBits(SpongeInstance *instance, unsigned char delimitedData)
{
    unsigned int rateInBytes = instance->rate/8;

    if (delimitedData == 0)
        return 1;
    if (instance->squeezing)
        return 1; /* Too late for additional input */

    /* Last few bits, whose delimiter coincides with first bit of padding */
    for(unsigned int instanceIndex = 0; instanceIndex < 4; instanceIndex++) {
        PlSnP_AddByte(instance->state, instanceIndex, delimitedData, instance->byteIOIndex);
    }

    /* If the first bit of padding is at position rate-1, we need a whole new block for the second bit of padding */
    if ((delimitedData >= 0x80) && (instance->byteIOIndex == (rateInBytes-1)))
        PlSnP_Permute(instance->state);
    /* Second bit of padding */
    for(unsigned int instanceIndex = 0; instanceIndex < 4; instanceIndex++) {
        PlSnP_AddByte(instance->state, instanceIndex, 0x80, rateInBytes - 1);
    }
    PlSnP_Permute(instance->state);
    instance->byteIOIndex = 0;
    instance->squeezing = 1;
    return 0;
}

/* ---------------------------------------------------------------- */

int SpongeSqueeze(SpongeInstance *instance, unsigned char **data, size_t dataByteLen)
{
    size_t i, j;
    unsigned int partialBlock;
    unsigned int rateInBytes = instance->rate/8;
    unsigned char *curData[4] = { NULL, NULL, NULL, NULL};

    if (!instance->squeezing)
        SpongeAbsorbLastFewBits(instance, 0x01);

    i = 0;
    if(dataByteLen > 0) {
        for (unsigned int instanceIndex = 0; instanceIndex < 4; instanceIndex++) {
            curData[instanceIndex] = data[instanceIndex];
        }
    }
    while(i < dataByteLen) {
        if ((instance->byteIOIndex == rateInBytes) && (dataByteLen >= (i + rateInBytes))) {
            for(j=dataByteLen-i; j>=rateInBytes; j-=rateInBytes) {
                PlSnP_Permute(instance->state);
                for(unsigned int instanceIndex = 0; instanceIndex < 4; instanceIndex++) {
                    PlSnP_ExtractBytes(instance->state, instanceIndex, curData[instanceIndex], 0, rateInBytes);
                    curData[instanceIndex]+=rateInBytes;
                }
            }
            i = dataByteLen - j;
        }
        else {
            /* normal lane: using the message queue */
            if (instance->byteIOIndex == rateInBytes) {
                PlSnP_Permute(instance->state);
                instance->byteIOIndex = 0;
            }
            partialBlock = (unsigned int)(dataByteLen - i);
            if (partialBlock+instance->byteIOIndex > rateInBytes)
                partialBlock = rateInBytes-instance->byteIOIndex;
            i += partialBlock;

            for(unsigned int instanceIndex = 0; instanceIndex < 4; instanceIndex++) {
                PlSnP_ExtractBytes(instance->state, instanceIndex, curData[instanceIndex], instance->byteIOIndex, partialBlock);
                curData[instanceIndex] += partialBlock;
            }
            instance->byteIOIndex += partialBlock;
        }
    }
    return 0;
}

/* ---------------------------------------------------------------- */

#undef Sponge
#undef SpongeInstance
#undef SpongeInitialize
#undef SpongeAbsorb
#undef SpongeAbsorbLastFewBits
#undef SpongeSqueeze
#undef PlSnP_statesSizeInBytes
#undef PlSnP_statesAlignment
#undef PlSnP_StaticInitialize
#undef PlSnP_InitializeAll
#undef PlSnP_AddByte
#undef PlSnP_AddBytes
#undef PlSnP_ExtractBytes
