/*
 * =====================================================================================
 *
 *       Filename:  TGM_FragLenTable.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  05/09/2012 04:08:49 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "TGM_Error.h"
#include "TGM_Utilities.h"
#include "TGM_FragLenTable.h"

using namespace Tangram;

FragLenTable::FragLenTable()
{

}

FragLenTable::~FragLenTable()
{

}

void FragLenTable::Read(FILE* fpHistInput)
{
    uint32_t numHist = 0;
    
    unsigned int readSize = fread(&numHist, sizeof(uint32_t), 1, fpHistInput);
    if (readSize != 1)
        TGM_ErrQuit("ERROR: Cannot read the number of hist arrays.\n");

    fragLenTable.Init(numHist);
    fragLenTable.SetSize(numHist);

    for (unsigned int i = 0; i != numHist; ++i)
    {
        uint32_t numElmnt = 0;
        readSize = fread(&(numElmnt), sizeof(uint32_t), 1, fpHistInput);
        if (readSize != 1)
            TGM_ErrQuit("ERROR: Cannot read the size of the fragment length histogram.\n");

        if (numElmnt > 0)
        {
            fragLenTable[i].size = numElmnt;

            fragLenTable[i].fragLen = (uint32_t*) calloc(sizeof(uint32_t), numElmnt);
            if (fragLenTable[i].fragLen == NULL)
                TGM_ErrQuit("ERROR: Not enough memory for the fragment length array.\n");

            readSize = fread(fragLenTable[i].fragLen, sizeof(uint32_t), numElmnt, fpHistInput);
            if (readSize != numElmnt)
                TGM_ErrQuit("ERROR: Cannot read the fragment length from the histogram file.\n");

            fragLenTable[i].freq = (uint64_t*) calloc(sizeof(uint64_t), numElmnt);
            if (fragLenTable[i].freq == NULL)
                TGM_ErrQuit("ERROR: Not enough memory for the frequency array.\n");

            readSize = fread(fragLenTable[i].freq , sizeof(uint64_t), numElmnt, fpHistInput);
            if (readSize != numElmnt)
                TGM_ErrQuit("ERROR: Cannot read the fragment length frequency from the histogram file.\n");

            for (unsigned int j = 1; j < numElmnt; ++j)
                fragLenTable[i].freq[j] += fragLenTable[i].freq[j - 1];
        }
        else
        {
            fragLenTable[i].fragLen = NULL;
            fragLenTable[i].freq = NULL;
            fragLenTable[i].size = 0;
        }
    }
}

int FragLenTable::GetQuality(uint32_t readGrpID, uint32_t targetFragLen, uint32_t median) const 
{
    const uint32_t* pFragLen = fragLenTable[readGrpID].fragLen;
    const uint64_t* pFreq = fragLenTable[readGrpID].freq;
    unsigned int size = fragLenTable[readGrpID].size;

    const uint32_t* fragLenPos = (const uint32_t*) bsearch(&targetFragLen, pFragLen, size, sizeof(uint32_t), CompareFragLen);

    int quality = -1;
    if (fragLenPos != NULL)
    {
        int index = fragLenPos - pFragLen;
        uint64_t totalFreq = pFreq[size - 1];
        double cdf = (double) pFreq[index] / totalFreq;

        if (targetFragLen > median)
            cdf = 1.0 - cdf;

        quality = DoubleRoundToInt(-10.0 * log10(cdf));
    }

    return quality;
}

void FragLenTable::Destory(void)
{
    unsigned int numHist = fragLenTable.Size();
    for (unsigned int i = 0; i != numHist; ++i)
    {
        free(fragLenTable[i].fragLen);
        free(fragLenTable[i].freq);
    }
}
