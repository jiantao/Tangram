/*
 * =====================================================================================
 *
 *       Filename:  TGM_LibTable.cpp
 *
 *    Description:  Library information table
 *
 *        Created:  05/08/2012 10:36:35 AM
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Jiantao Wu (), 
 *   Inistitution:  Boston College
 *
 * =====================================================================================
 */

#include <cstdio>

#include "khash.h"
#include "TGM_Error.h"
#include "TGM_LibTable.h"

using namespace std;
using namespace Tangram;

#define MAX_LIB_DIFF 350

KHASH_MAP_INIT_STR(name, uint32_t);

LibTable::LibTable()
{
    specialRefHash = NULL;
    readGrpHash = NULL;
}

LibTable::~LibTable()
{
    khash_t(name)* pHash = (khash_t(name)*) readGrpHash;
    kh_destroy(name, pHash);

    pHash = (khash_t(name)*) specialRefHash;
    kh_destroy(name, pHash);

    unsigned int size = readGrpNames.Size();
    for (unsigned int i = 0; i != size; ++i)
        free(readGrpNames[i]);

    size = sampleNames.Size();
    for (unsigned int i = 0; i != size; ++i)
        free(sampleNames[i]);

    size = anchorNames.Size();
    for (unsigned int i = 0; i != size; ++i)
        free(anchorNames[i]);

    size = specialRefNames.Size();
    for (unsigned int i = 0; i != size; ++i)
        free(specialRefNames[i]);
}

bool LibTable::Read(FILE* fpLibInput, int maxFragDiff)
{
    unsigned int readSize = 0;
    uint32_t sizeAC = 0;
    uint32_t sizeSM = 0;
    uint32_t sizeRG = 0;
    uint32_t sizeSP = 0;

    readSize = fread(&sizeAC, sizeof(uint32_t), 1, fpLibInput);
    if (readSize != 1)
        TGM_ErrQuit("ERROR: Cannot read the number of anchors from the library file.\n");

    readSize = fread(&sizeSM, sizeof(uint32_t), 1, fpLibInput);
    if (readSize != 1)
        TGM_ErrQuit("ERROR: Cannot read the number of sample from the library file.\n");

    readSize = fread(&sizeRG, sizeof(uint32_t), 1, fpLibInput);
    if (readSize != 1)
        TGM_ErrQuit("ERROR: Cannot read the number of read group from the library file.\n");

    ReadAnchors(sizeAC, fpLibInput);
    ReadSamples(sizeSM, fpLibInput);
    ReadReadGrps(sizeRG, fpLibInput);

    if (maxFragDiff > 0)
        CheckReadGrps(maxFragDiff);

    readSize = fread(&sizeSP, sizeof(uint32_t), 1, fpLibInput);
    if (readSize != 1)
        TGM_ErrQuit("ERROR: Cannot read the number of special references from the library file.\n");

    if (sizeSP > 0)
    {
        ReadSpecialRef(sizeSP, fpLibInput);
    }

    return true;
}

bool LibTable::GetReadGrpID(uint32_t& readGrpID, const char* readGrpName) const
{
    const khash_t(name)* pHash = (const khash_t(name)*) readGrpHash;
    khiter_t khIter = kh_get(name, pHash, readGrpName);

    if (khIter != kh_end(pHash))
    {
        readGrpID = kh_value(pHash, khIter);
        return true;
    }
    else
        return false;
}

bool LibTable::GetSampleID(uint32_t& sampleID, uint32_t readGrpID) const
{
    if (readGrpID <= readGrpNames.Size())
    {
        sampleID = sampleMap[readGrpID];
        return true;
    }
    else
        return false;
}

bool LibTable::GetSpecialRefID(uint32_t& specialRefID, const char* specialRefName) const
{
    const khash_t(name)* pHash = (const khash_t(name)*) specialRefHash;
    khiter_t khIter = kh_get(name, pHash, specialRefName);

    if (khIter != kh_end(pHash))
    {
        specialRefID = kh_value(pHash, khIter);
        return true;
    }
    else
        return false;
}

void LibTable::ReadAnchors(uint32_t sizeAC, FILE* fpLibInput)
{
    unsigned int md5Len = sizeAC * MD5_STR_LEN;

    anchorNames.Init(sizeAC);
    anchorLength.Init(sizeAC);
    md5.Init(md5Len);

    anchorNames.SetSize(sizeAC);
    anchorLength.SetSize(sizeAC);
    md5.SetSize(md5Len);

    unsigned int readSize = 0;

    int32_t* pAnchorLength = anchorLength.GetPointer(0);
    readSize = fread(pAnchorLength, sizeof(int32_t), sizeAC, fpLibInput);
    if (readSize != sizeAC)
        TGM_ErrQuit("ERROR: Cannot read the length of anchors from the library file.\n");

    char* pMd5 = md5.GetPointer(0);
    readSize = fread(pMd5, sizeof(char), md5Len, fpLibInput);
    if (readSize != md5Len)
        TGM_ErrQuit("ERROR: Cannot read the md5 strings from the library file.\n");

    uint32_t nameLen = 0;
    for (unsigned int i = 0; i != sizeAC; ++i)
    {
        readSize = fread(&nameLen, sizeof(uint32_t), 1, fpLibInput);
        if (readSize != 1)
            TGM_ErrQuit("ERROR: Cannot read the length of anchor name from the library file.\n");

        anchorNames[i] = (char*) malloc((nameLen + 1) * sizeof(char));
        if (anchorNames[i] == NULL)
            TGM_ErrQuit("ERROR: Not enough memory for the storage of the anchor name.\n");

        anchorNames[i][nameLen] = '\0';

        readSize = fread(anchorNames[i], sizeof(char), nameLen, fpLibInput);
        if (readSize != nameLen)
            TGM_ErrQuit("ERROR: Cannot read the anchor name from the library file.\n");
    }

    uint32_t specialPrefixLen = 0;
    readSize = fread(&(specialPrefixLen), sizeof(uint32_t), 1, fpLibInput);
    if (readSize != 1)
        TGM_ErrQuit("ERROR: Cannot read the length of the special prefix from the library file.\n");

    if (specialPrefixLen > 0)
    {
        char prefix[specialPrefixLen + 1];
        prefix[specialPrefixLen] = '\0';

        readSize = fread(prefix, sizeof(char), specialPrefixLen, fpLibInput);
        if (readSize != specialPrefixLen)
            TGM_ErrQuit("ERROR: Cannot read the special prefix from the library file.\n");

        specialPrefix = prefix;
    }
}

void LibTable::ReadSamples(uint32_t sizeSM, FILE* fpLibInput)
{
    unsigned int readSize = 0;
    uint32_t nameLen = 0;

    sampleNames.Init(sizeSM);
    sampleNames.SetSize(sizeSM);

    for (unsigned int i = 0; i != sizeSM; ++i)
    {
        readSize = fread(&nameLen, sizeof(uint32_t), 1, fpLibInput);
        if (readSize != 1)
            TGM_ErrQuit("ERROR: Cannot read the length of sample name from the library file.\n");

        sampleNames[i] = (char*) malloc((nameLen + 1) * sizeof(char));
        if (sampleNames[i] == NULL)
            TGM_ErrQuit("ERROR: Not enough memory for the storage of the sample name.\n");

        sampleNames[i][nameLen] = '\0';
        
        readSize = fread(sampleNames[i], sizeof(char), nameLen, fpLibInput);
        if (readSize != nameLen)
            TGM_ErrQuit("ERROR: Cannot read the sample name from the library file.\n");
    }
}

void LibTable::ReadReadGrps(uint32_t sizeRG, FILE* fpLibInput)
{
    uint32_t nameLen = 0;
    unsigned int readSize = 0;

    sampleMap.Init(sizeRG);
    sampleMap.SetSize(sizeRG);

    int32_t* pSampleMap = sampleMap.GetPointer(0);
    readSize = fread(pSampleMap, sizeof(int32_t), sizeRG, fpLibInput);
    if (readSize != sizeRG)
        TGM_ErrQuit("ERROR: Cannot read the sample map from the library file.\n");

    readGrpNames.Init(sizeRG);
    readGrpNames.SetSize(sizeRG);

    readGrpHash = kh_init(name);

    for (unsigned int i = 0; i != sizeRG; ++i)
    {
        readSize = fread(&nameLen, sizeof(uint32_t), 1, fpLibInput);
        if (readSize != 1)
            TGM_ErrQuit("ERROR: Cannot read the length of read group name from the library file.\n");

        readGrpNames[i] = (char*) malloc((nameLen + 1) * sizeof(char));
        if (readGrpNames[i] == NULL)
            TGM_ErrQuit("ERROR: Not enough memory for the storage of the sample name.\n");

        readGrpNames[i][nameLen] = '\0';

        readSize = fread(readGrpNames[i], sizeof(char), nameLen, fpLibInput);
        if (readSize != nameLen)
            TGM_ErrQuit("ERROR: Cannot read the read group name from the library file.\n");

        UpdateReadGrpHash(i, readGrpNames[i]);
    }

    readSize = fread(&(fragLenMax), sizeof(uint32_t), 1, fpLibInput);
    if (readSize != 1)
        TGM_ErrQuit("ERROR: Cannot read the maximum fragment length from the library file.\n");

    if (fragLenMax > 0)
    {
        readSize = fread(&(cutoff), sizeof(double), 1, fpLibInput);
        if (readSize != 1)
            TGM_ErrQuit("ERROR: Cannot read the cutoff from the library file.\n");

        readSize = fread(&(trimRate), sizeof(double), 1, fpLibInput);
        if (readSize != 1)
            TGM_ErrQuit("ERROR: Cannot read the trim rate from the library file.\n");

        seqTech.Init(sizeRG);
        seqTech.SetSize(sizeRG);

        int8_t* pSeqTech = seqTech.GetPointer(0);
        readSize = fread(pSeqTech, sizeof(int8_t), sizeRG, fpLibInput);
        if (readSize != sizeRG)
            TGM_ErrQuit("ERROR: Cannot read the sequencing technology from the library file.\n");

        libInfo.Init(sizeRG);
        libInfo.SetSize(sizeRG);

        LibInfo* pLibInfo = libInfo.GetPointer(0);
        readSize = fread(pLibInfo, sizeof(LibInfo), sizeRG, fpLibInput);
        if (readSize != sizeRG)
            TGM_ErrQuit("ERROR: Cannot read the library information from the library file.\n");
    }
}

// get rid of those bad libraries
void LibTable::CheckReadGrps(int maxFragDiff)
{
    unsigned int size = libInfo.Size();
    int32_t newMaxFrag = 0;
    for (unsigned int i = 0; i != size; ++i)
    {
        if (seqTech[i] == ST_ILLUMINA)
        {
            if (libInfo[i].fragLenHigh - libInfo[i].fragLenLow > maxFragDiff)
            {
                libInfo[i].fragLenMedian = 0;
                libInfo[i].fragLenHigh = 0;
            }
            else
            {
                if (libInfo[i].fragLenHigh > newMaxFrag)
                    newMaxFrag = libInfo[i].fragLenHigh;
            }
        }
    }

    if (newMaxFrag == 0)
        TGM_ErrQuit("ERROR: No valid read groups after filtering.\n");

    fragLenMax = newMaxFrag;
}

void LibTable::ReadSpecialRef(uint32_t sizeSP, FILE* fpLibInput)
{
    specialRefNames.Init(sizeSP);
    specialRefNames.SetSize(sizeSP);

    specialRefHash = kh_init(name);

    unsigned int readSize = 0;
    for (unsigned int i = 0; i != sizeSP; ++i)
    {
        specialRefNames[i] = (char*) calloc(sizeof(char), 3);
        readSize = fread(specialRefNames[i], sizeof(char), 2, fpLibInput);
        if (readSize != 2)
            TGM_ErrQuit("ERROR: Cannot read the special reference ID into the library information file.\n");

        UpdateSpecialRefHash(i, specialRefNames[i]);
    }
}

void LibTable::UpdateReadGrpHash(uint32_t readGrpID, const char* readGrpName)
{
    khash_t(name)* pHash = (khash_t(name)*) readGrpHash;

    int ret = 0;
    khiter_t khIter = kh_put(name, pHash, readGrpName, &ret);

    if (ret != 0)
        kh_value(pHash, khIter) = readGrpID;
}

void LibTable::UpdateSpecialRefHash(uint32_t specialRefID, const char* specialRefName)
{
    khash_t(name)* pHash = (khash_t(name)*) specialRefHash;

    int ret = 0;
    khiter_t khIter = kh_put(name, pHash, specialRefName, &ret);

    if (ret != 0)
        kh_value(pHash, khIter) = specialRefID;
}
