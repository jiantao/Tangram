/*
 * =====================================================================================
 *
 *       Filename:  TGM_FragLenHist.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  02/08/2012 02:24:40 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include "khash.h"
#include "TGM_Error.h"
#include "TGM_Utilities.h"
#include "TGM_FragLenHist.h"

// index of the mode count for invalid pair mode 
#define INVALID_PAIR_MODE_SET_INDEX 1

#define DEFAULT_NUM_HIST_ELMNT 200

// fragment length hash
KHASH_MAP_INIT_INT(fragLen, uint64_t);


static inline int CompareFragLenBin(const void* a, const void* b)
{
    const uint32_t* first = a;
    const uint32_t* second = b;

    if (*first < *second)
        return -1;
    else if (*first > *second)
        return 1;
    else
        return 0;
}

static void TGM_FragLenHistToMature(TGM_FragLenHist* pHist)
{
    khash_t(fragLen)* pRawHist = pHist->rawHist;

    pHist->size = kh_size(pRawHist);

    if (pHist->size > pHist->capacity)
    {
        free(pHist->fragLen);
        free(pHist->freq);

        pHist->capacity = 2 * pHist->size;

        pHist->fragLen = (uint32_t*) malloc(pHist->capacity * sizeof(uint32_t));
        if(pHist->fragLen == NULL)
            TGM_ErrQuit("ERROR: Not enough memory for the storage of the fragment length array in the fragment length histogram object.\n");

        pHist->freq = (uint64_t*) malloc(pHist->capacity * sizeof(uint64_t));
        if(pHist->freq == NULL)
            TGM_ErrQuit("ERROR: Not enough memory for the storage of the frequency array in the fragment length histogram object.\n");
    }

    unsigned int i = 0;
    for (khiter_t khIter = kh_begin(pRawHist); khIter != kh_end(pRawHist); ++khIter)
    {
        if (kh_exist(pRawHist, khIter))
        {
            pHist->fragLen[i] = kh_key(pRawHist, khIter);
            ++i;
        }
    }

    qsort(pHist->fragLen, pHist->size, sizeof(uint32_t), CompareFragLenBin);

    double cumFreq = 0.0;
    double totalFragLen = 0.0;
    uint64_t totalFreq = pHist->modeCount[0];
    double cdf = 0;
    uint32_t fragLenQual = 0;

    TGM_Bool foundMedian = FALSE;
    for (unsigned int j = 0; j != pHist->size; ++j)
    {
        khiter_t khIter = kh_get(fragLen, pRawHist, pHist->fragLen[j]);
        if (khIter == kh_end(pRawHist))
            TGM_ErrQuit("ERROR: Cannot find the fragment length frequency from the hash table.\n");

        pHist->freq[j] = kh_value(pRawHist, khIter);

        totalFragLen += pHist->fragLen[j] * pHist->freq[j];
        cumFreq += pHist->freq[j];
        cdf = cumFreq / totalFreq;

        if (!foundMedian && cdf >= 0.5)
        {
            pHist->median = pHist->fragLen[j];
            foundMedian = TRUE;
        }

        cdf = cdf > 0.5 ? 1.0 - cdf : cdf;

        fragLenQual = DoubleRoundToInt(-10.0 * log10(cdf));
        kh_value(pRawHist, khIter) = fragLenQual;
    }

    pHist->mean = totalFragLen / totalFreq;

    pHist->stdev = 0.0;
    for (unsigned int j = 0; j != pHist->size; ++j)
    {
        pHist->stdev += (double) pHist->freq[j] * pow(pHist->mean - pHist->fragLen[j], 2);
    }

    if (totalFreq != 1)
        pHist->stdev = sqrt(pHist->stdev / (double) (totalFreq - 1));
}

TGM_FragLenHistArray* TGM_FragLenHistArrayAlloc(unsigned int capacity)
{
    TGM_FragLenHistArray* pHistArray = NULL;
    TGM_ARRAY_ALLOC(pHistArray, capacity, TGM_FragLenHistArray, TGM_FragLenHist);

    return pHistArray;
}

void TGM_FragLenHistArrayFree(TGM_FragLenHistArray* pHistArray)
{
    if (pHistArray != NULL)
    {
        for (unsigned int i = 0; i != pHistArray->capacity; ++i)
        {
            free(pHistArray->data[i].fragLen);
            free(pHistArray->data[i].freq);

            kh_destroy(fragLen, pHistArray->data[i].rawHist);
        }

        free(pHistArray->data);
        free(pHistArray);
    }
}

TGM_FragLenHistLite* TGM_FragLenHistLiteAlloc(unsigned int capacity)
{
    TGM_FragLenHistLite* pHistLite = (TGM_FragLenHistLite*) malloc(sizeof(TGM_FragLenHistLite));
    if (pHistLite == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for a fragment length histogram object.\n");

    pHistLite->size = 0;
    pHistLite->capacity = 0;

    pHistLite->fragLen = NULL;
    pHistLite->freq = NULL;

    TGM_FragLenHistLiteInit(pHistLite, (capacity + 1) / 2);

    return pHistLite;
}

void TGM_FragLenHistLiteFree(TGM_FragLenHistLite* pHistLite)
{
    if (pHistLite != NULL)
    {
        free(pHistLite->fragLen);
        free(pHistLite->freq);

        free(pHistLite);
    }
}

void TGM_FragLenHistLiteArrayFree(TGM_FragLenHistLiteArray* pHistLiteArray)
{
    if (pHistLiteArray != NULL)
    {
        for (unsigned int i = 0; i != pHistLiteArray->size; ++i)
        {
            free(pHistLiteArray->data[i].fragLen);
            free(pHistLiteArray->data[i].freq);
        }

        free(pHistLiteArray->data);
        free(pHistLiteArray);
    }
}

void TGM_FragLenHistArrayClear(TGM_FragLenHistArray* pHistArray)
{
    for (unsigned int i = 0; i != pHistArray->size; ++i)
    {
        pHistArray->data[i].size = 0;
        pHistArray->data[i].mean = 0.0;
        pHistArray->data[i].median = 0.0;
        pHistArray->data[i].stdev = 0.0;
        pHistArray->data[i].modeCount[0] = 0;
        pHistArray->data[i].modeCount[1] = 0;

        kh_clear(fragLen, pHistArray->data[i].rawHist);
    }
}

void TGM_FragLenHistArrayInit(TGM_FragLenHistArray* pHistArray, unsigned int newSize)
{
    TGM_FragLenHistArrayClear(pHistArray);

    if (newSize > pHistArray->capacity)
    {
        TGM_ARRAY_RESIZE(pHistArray, newSize * 2, TGM_FragLenHist);
        memset(pHistArray->data + pHistArray->size, 0, (newSize * 2 - pHistArray->size) * sizeof(TGM_FragLenHist));
    }

    for (unsigned int i = 0; i != newSize; ++i)
    {
        if (pHistArray->data[i].rawHist == NULL)
        {
            pHistArray->data[i].rawHist = kh_init(fragLen);
            kh_resize(fragLen, pHistArray->data[i].rawHist, DEFAULT_NUM_HIST_ELMNT);
        }
    }

    pHistArray->size = newSize;
}

TGM_Status TGM_FragLenHistArrayUpdate(TGM_FragLenHistArray* pHistArray, unsigned int backHistIndex, uint32_t fragLen)
{
    if (backHistIndex > pHistArray->size)
        return TGM_ERR;

    TGM_FragLenHist* pCurrHist = pHistArray->data + (pHistArray->size - backHistIndex);

    // if the pair mode is not valid
    // we only updated the count of the invalid pair and return
    if (fragLen == 0)
    {
        ++(pCurrHist->modeCount[INVALID_PAIR_MODE_SET_INDEX]);
        return TGM_OK;
    }

    // because we can have up to 2 different pair mode sets (4 differen pair modes)
    // we should choose which histogram we should update
    // the first one (with index 0 or 1) or the second one(2, 3)
    khash_t(fragLen)* pCurrHash = pCurrHist->rawHist;

    if (pCurrHash == NULL)
    {
        pCurrHash = kh_init(fragLen);
        kh_resize(fragLen, pCurrHash, 20);
    }

    int ret = 0;
    khiter_t khIter = kh_put(fragLen, pCurrHash, fragLen, &ret);

    if (ret == 0)
        kh_value(pCurrHash, khIter) += 1;
    else
        kh_value(pCurrHash, khIter) = 1;

    ++(pCurrHist->modeCount[0]);
    pCurrHist->rawHist = pCurrHash;

    return TGM_OK;
}

int TGM_FragLenHistLiteArrayGetFragLenQual(const TGM_FragLenHistLiteArray* pHistArray, uint32_t readGrpID, uint32_t fragLen, uint32_t median)
{
    TGM_FragLenHistLite* pCurrHist = pHistArray->data + readGrpID;

    double totalFreq = pCurrHist->freq[pCurrHist->size - 1];
    uint32_t* pTarget = bsearch(&fragLen, pCurrHist->fragLen, pCurrHist->size, sizeof(uint32_t), TGM_CompareUint);

    int fragLenQual = -1;
    if (pTarget != NULL)
    {
        unsigned int index = pTarget - pCurrHist->fragLen;
        double cdf = (double) pCurrHist->freq[index] / totalFreq;

        if (fragLen > median)
            cdf = 1.0 - cdf;

        fragLenQual = DoubleRoundToInt(-10.0 * log10(cdf));
    }

    return fragLenQual;
}

void TGM_FragLenHistArrayFinalize(TGM_FragLenHistArray* pHistArray)
{
    for (unsigned int i = 0; i != pHistArray->size; ++i)
        TGM_FragLenHistToMature(&(pHistArray->data[i]));
}

void TGM_FragLenHistArrayWriteHeader(uint32_t size, FILE* output)
{
    int ret = fseeko(output, 0, SEEK_SET);
    if (ret != 0)
        TGM_ErrQuit("ERROR: Cannot seek the histogram file.\n");

    unsigned int writeSize = fwrite(&size, sizeof(uint32_t), 1, output);
    if (writeSize != 1)
        TGM_ErrQuit("ERROR: Cannot write the number of histogram into the file.\n");
}

void TGM_FragLenHistArrayWrite(const TGM_FragLenHistArray* pHistArray, FILE* output)
{
    unsigned int writeSize = 0;
    for (unsigned int i = 0; i != pHistArray->size; ++i)
    {
        // ignore those read groups that have too few read pairs
        uint32_t histSize = pHistArray->data[i].size >= MIN_FRAGLEN_HIST_SIZE ? pHistArray->data[i].size : 0;

        writeSize = fwrite(&(histSize), sizeof(uint32_t), 1, output);
        if (writeSize != 1)
            TGM_ErrQuit("ERROR: Cannot write the size of the histogram into the file.\n");

        if (histSize != 0)
        {
            writeSize = fwrite(pHistArray->data[i].fragLen, sizeof(uint32_t), histSize, output);
            if (writeSize != histSize)
                TGM_ErrQuit("ERROR: Cannot write the fragment length of the histogram into the file.\n");

            writeSize = fwrite(pHistArray->data[i].freq, sizeof(uint64_t), histSize, output);
            if (writeSize != histSize)
                TGM_ErrQuit("ERROR: Cannot write the frequency of the histogram into the file.\n");
        }
    }

    fflush(output);
}

TGM_FragLenHistLiteArray* TGM_FragLenHistLiteArrayRead(FILE* pHistArrayInput)
{
    uint32_t numHist = 0;
    unsigned int readSize = fread(&numHist, sizeof(uint32_t), 1, pHistArrayInput);
    if (readSize != 1)
        TGM_ErrQuit("ERROR: Cannot read the number of hist arrays.\n");

    TGM_FragLenHistLiteArray* pHistArray = NULL;

    TGM_ARRAY_ALLOC(pHistArray, numHist, TGM_FragLenHistLiteArray, TGM_FragLenHistLite);
    for (unsigned int i = 0; i != numHist; ++i)
    {
        TGM_FragLenHistLite* pHistLite = pHistArray->data + i;
        readSize = fread(&(pHistLite->size), sizeof(uint32_t), 1, pHistArrayInput);
        if (readSize != 1)
            TGM_ErrQuit("ERROR: Cannot read the size of the fragment length histogram.\n");

        pHistLite->capacity = pHistLite->size;
        if (pHistLite->size != 0)
        {
            pHistLite->fragLen = (uint32_t*) malloc(sizeof(uint32_t) * pHistLite->capacity);
            if (pHistLite->fragLen == NULL)
                TGM_ErrQuit("ERROR: Not enough memory for the fragment length array.\n");

            pHistLite->freq = (uint64_t*) malloc(sizeof(uint64_t) * pHistLite->capacity);
            if (pHistLite->freq == NULL)
                TGM_ErrQuit("ERROR: Not enough memory for the frequecy array.\n");

            readSize = fread(pHistLite->fragLen, sizeof(uint32_t), pHistLite->size, pHistArrayInput);
            if (readSize != pHistLite->size)
                TGM_ErrQuit("ERROR: Cannot read the fragment length from the histogram file.\n");

            readSize = fread(pHistLite->freq, sizeof(uint64_t), pHistLite->size, pHistArrayInput);
            if (readSize != pHistLite->size)
                TGM_ErrQuit("ERROR: Cannot read the fragment length frequency from the histogram file.\n");

            for (unsigned int j = 1; j < pHistLite->size; ++j)
                pHistLite->freq[j] += pHistLite->freq[j - 1];
        }
        else
        {
            pHistLite->fragLen = NULL;
            pHistLite->freq = NULL;
        }
    }

    pHistArray->size = numHist;

    return pHistArray;
}


void TGM_FragLenHistLiteInit(TGM_FragLenHistLite* pHistLite, uint32_t newSize)
{
    if (newSize > pHistLite->capacity)
    {
        pHistLite->capacity = 2 * newSize;
        pHistLite->fragLen = (uint32_t*) realloc(pHistLite->fragLen, sizeof(uint32_t) * pHistLite->capacity);
        if (pHistLite->fragLen == NULL)
            TGM_ErrQuit("ERROR: Not enough memory for the fragment length array.\n");

        pHistLite->freq = (uint64_t*) realloc(pHistLite->freq, sizeof(uint64_t) * pHistLite->capacity);
        if (pHistLite->freq == NULL)
            TGM_ErrQuit("ERROR: Not enough memory for the frequecy array.\n");
    }
}

void TGM_FragLenHistLiteRead(TGM_FragLenHistLite* pHistLite, FILE* input)
{
    unsigned int readSize = fread(&(pHistLite->size), sizeof(uint32_t), 1, input);
    if (readSize != 1)
        TGM_ErrQuit("ERROR: Cannot read the size of histogram from the fragment length histogram file.\n");

    if (pHistLite->size != 0)
    {
        TGM_FragLenHistLiteInit(pHistLite, pHistLite->size);

        readSize = fread(pHistLite->fragLen, sizeof(uint32_t), pHistLite->size, input);
        if (readSize != pHistLite->size)
            TGM_ErrQuit("ERROR: Cannot read the fragment length array from the file.\n");

        readSize = fread(pHistLite->freq, sizeof(uint64_t), pHistLite->size, input);
        if (readSize != pHistLite->size)
            TGM_ErrQuit("ERROR: Cannot read the frequency array from the file.\n");
    }
}

void TGM_FragLenHistLiteWrite(const TGM_FragLenHistLite* pHistLite, FILE* output)
{
    unsigned int writeSize = fwrite(&(pHistLite->size), sizeof(uint32_t), 1, output);
    if (writeSize != 1)
        TGM_ErrQuit("ERROR: Cannot write the size of histogram.\n");

    if (pHistLite->size != 0)
    {
        writeSize = fwrite(pHistLite->fragLen, sizeof(uint32_t), pHistLite->size, output);
        if (writeSize != pHistLite->size)
            TGM_ErrQuit("ERROR: Cannot write the fragment length array into the file.\n");

        writeSize = fwrite(pHistLite->freq, sizeof(uint64_t), pHistLite->size, output);
        if (writeSize != pHistLite->size)
            TGM_ErrQuit("ERROR: Cannot write the fragment length array into the file.\n");
    }
}

void TGM_FragLenHistTransfer(TGM_FragLenHistLite* pHistLite, FILE* input, FILE* output)
{
    uint32_t numHist = 0;
    unsigned int readSize = fread(&numHist, sizeof(uint32_t), 1, input);
    if (readSize != 1)
        TGM_ErrQuit("ERROR: Cannot read the number of histograms from the fragment length histogram file.\n");

    for (unsigned int i = 0; i != numHist; ++i)
    {
        TGM_FragLenHistLiteRead(pHistLite, input);
        TGM_FragLenHistLiteWrite(pHistLite, output);
    }

    fflush(output);
}
