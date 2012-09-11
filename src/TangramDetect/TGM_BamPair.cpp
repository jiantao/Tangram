/*
 * =====================================================================================
 *
 *       Filename:  TGM_BamPair.cpp
 *
 *    Description:  Table used to store bam pairs
 *
 *        Created:  05/09/2012 06:52:08 PM
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Jiantao Wu (), 
 *   Inistitution:  Boston College
 *
 * =====================================================================================
 */

#ifndef __STDC_CONSTANT_MACROS
#define __STDC_CONSTANT_MACROS
#endif

#ifndef __STDC_LIMIT_MACROS
#define __STDC_LIMIT_MACROS
#endif

#include <stdlib.h>
#include <stdint.h>
#include <cmath>
#include <cstring>

#include "khash.h"
#include "TGM_Utilities.h"
#include "TGM_BamPair.h"

using namespace std;
using namespace BamTools;
using namespace Tangram;

// KHASH_MAP_INIT_STR(name, ObjIndex);

BamPairTable::BamPairTable(const LibTable& inLibTable, const FragLenTable& inFraglenTable, uint16_t minSoftSize, uint8_t minMQ, uint8_t minSpMQ)
              : libTable(inLibTable), fragLenTable(inFraglenTable), minSoftSize(minSoftSize), minMQ(minMQ), minSpMQ(minSpMQ)
{
    pairStat.spRef[0][2] = '\0';
    pairStat.spRef[1][2] = '\0';

    // readNameHash = kh_init(name);
    longPairs.Init(10);
    invertedPairs.Init(10);
    specialPairs.Init(10);
    orphanPairs.Init(10);
    softPairs.Init(10);

    numInverted3 = 0;
}

BamPairTable::~BamPairTable()
{
    // kh_destroy(name, (khash_t(name)*) readNameHash);
    unsigned int size = orphanPairs.Size();
    for (unsigned int i = 0; i != size; ++i)
        TGM_SeqClean(&(orphanPairs[i].read));

    size = softPairs.Size();
    for (unsigned int i = 0; i != size; ++i)
    {
        free(softPairs[i].cigar);
        TGM_SeqClean(&(softPairs[i].read));
    }
}

void BamPairTable::Update(const BamAlignment& alignment)
{
    pAlignment = &alignment;

    if (BamPairFilter())
    {
        if (IsOrphanPair())
        {
            UpdateOrphanPair();
        }
        else if (SetPairStat())
        {
            switch (pairStat.readPairType)
            {
                case PT_NORMAL:
                case PT_UNKNOWN:
                    break;
                case PT_LONG:
                    UpdateLocalPair(longPairs);
                    break;
                case PT_SHORT:
                    // UpdateLocalPair(shortPairs);
                    break;
                case PT_REVERSED:
                    // UpdateLocalPair(reversedPairs);
                    break;
                case PT_INVERTED3:
                case PT_INVERTED5:
                    UpdateInvertedPair();
                    break;
                case PT_SPECIAL3:
                case PT_SPECIAL5:
                    UpdateSpecialPair();
                    break;
                case PT_CROSS:
                    UpdateCrossPair();
                    break;
                case PT_SOFT3:
                case PT_SOFT5:
                    UpdateSoftPair();
                    break;
                default:
                    break;
            }
        }
    }
}

bool BamPairTable::BamPairFilter(void) const
{
    // filter out those alignment that:
    // 1. is not paired
    // 2. both mates unmapped
    // 3. either of the mate reference ID is negative
    // 4. failed the quality control
    // 5. is duplcated marked
    // 6. is not primary alignment
    // 7. name is not valid
    if (!(*pAlignment).IsPaired() 
        || (!pAlignment->IsMapped() && !pAlignment->IsMateMapped()) 
        || (*pAlignment).RefID < 0 || pAlignment->MateRefID < 0
        || pAlignment->IsFailedQC() || pAlignment->IsDuplicate() || !pAlignment->IsPrimaryAlignment()
        || pAlignment->Name == "0" || pAlignment->Name == "*") 
    {
        return false;
    }

    return true;
}

bool BamPairTable::SetPairStat(void)
{
    // initialize the special reference name
    pairStat.spRef[0][0] = ' ';
    pairStat.spRef[1][0] = ' ';

    // initialize the soft size
    memset(pairStat.softSize, 0, sizeof(int) * 4);

    // boolean variable to indicate if this is the up mate (mate at upstream)
    bool isUpMate = true;

    if (pAlignment->RefID == pAlignment->MateRefID)
    {
        if (pAlignment->Position > pAlignment->MatePosition)
            isUpMate = false;
    }
    else if (pAlignment->RefID > pAlignment->MateRefID)
        isUpMate = false;

    // parse the ZA string
    if (pAlignment->GetTag("ZA", zaStr))
    {
        ParseZAstr(isUpMate);
    }
    else
    {
        // FIXME: Need to add some support for no za tag
        return false;
    }

    if ((isUpMate && pairStat.spRef[0][0] != ' ') || (!isUpMate && pairStat.spRef[1][0] != ' '))
        return false;

    // get the pair orientation
    pairStat.orient = GetPairOrient(isUpMate);
    if (pairStat.orient == TGM_BAD_PAIR_MODE)
        return false;

    if (pAlignment->GetTag("RG", readGrpName))
    {
        if (!libTable.GetReadGrpID(pairStat.readGrpID, readGrpName.c_str()))
            return false;

        // filter out those pairs whose libraries are not good
        unsigned int fragLenMedian = libTable.GetFragLenMedian(pairStat.readGrpID);
        if (fragLenMedian == 0)
            return false;

        // get the sequencing technology
        SeqTech seqTech = libTable.GetSeqTech(pairStat.readGrpID);
        pairStat.orient = (PairOrient) TGM_SeqTechMap[seqTech][pairStat.orient];
    }
    else
        return false;

    if (pAlignment->RefID == pAlignment->MateRefID)
        pairStat.fragLen = abs(pAlignment->InsertSize);
    else
        pairStat.fragLen = -1;

    pairStat.readPairType = GetPairType(isUpMate);

    return true;
}

void BamPairTable::ParseZAstr(bool isUpMate)
{
    unsigned char whichMate = 0;
    int numMappings = 0;
    const char* currFieldPos = zaStr.c_str() + 1;

    if ((*currFieldPos == '&' && isUpMate) || (*currFieldPos == '@' && !isUpMate))
        whichMate = 1;

    currFieldPos += 2;
    unsigned int currFieldNum = 1;
    const char* currFieldEnd = NULL;

    const char* cigarStr = NULL;
    unsigned int cigarLen = 0;

    const char* mdStr = NULL;
    unsigned int mdLen = 0;

    do
    {
        ++currFieldNum;
        currFieldEnd = strpbrk(currFieldPos, ";>");

        if (currFieldEnd - currFieldPos == 0)
        {
            ++currFieldPos;
            continue;
        }

        switch (currFieldNum)
        {
            case 2:
                pairStat.bestMQ[whichMate] = atoi(currFieldPos);
                break;
            case 3:
                pairStat.secMQ[whichMate] = atoi(currFieldPos);
                break;
            case 4:
                pairStat.spRef[whichMate][0] = *currFieldPos;
                pairStat.spRef[whichMate][1] = *(currFieldPos + 1);
                break;
            case 5:
                numMappings = atoi(currFieldPos);
                pairStat.numMappings[whichMate] = numMappings;
                if (numMappings > UINT16_MAX)
                    pairStat.numMappings[whichMate] = UINT16_MAX;
                break;
            case 6:
                cigarStr = currFieldPos;
                cigarLen = currFieldEnd - currFieldPos;
                break;
            case 7:
                mdStr = currFieldPos;
                mdLen = currFieldEnd - currFieldPos;
                break;
            default:
                break;
        }

        currFieldPos = currFieldEnd + 1;

    }while(*currFieldEnd != '>');

    if (cigarStr != NULL && mdStr != NULL)
    {
        pairStat.numMM[whichMate] = GetNumMismatchFromZA(&(pairStat.end[whichMate]), pairStat.softSize[whichMate], cigarStr, cigarLen, mdStr, mdLen);
        pairStat.end[whichMate] = pairStat.end[whichMate] + pAlignment->MatePosition - 1;
    }
    else
    {
        pairStat.numMM[whichMate] = GetNumMismatchFromBam(pairStat.softSize[whichMate]);
        pairStat.end[whichMate] = pAlignment->GetEndPosition(false,  true);
    }


    // handle the sencond ZA tag
    currFieldPos += 3;
    currFieldNum = 1;
    currFieldEnd = NULL;

    cigarStr = NULL;
    cigarLen = 0;

    mdStr = NULL;
    mdLen = 0;

    whichMate ^= 1;

    do
    {
        ++currFieldNum;
        currFieldEnd = strpbrk(currFieldPos, ";>");

        if (currFieldEnd - currFieldPos == 0)
        {
            ++currFieldPos;
            continue;
        }

        switch (currFieldNum)
        {
            case 2:
                pairStat.bestMQ[whichMate] = atoi(currFieldPos);
                break;
            case 3:
                pairStat.secMQ[whichMate] = atoi(currFieldPos);
                break;
            case 4:
                pairStat.spRef[whichMate][0] = *currFieldPos;
                pairStat.spRef[whichMate][1] = *(currFieldPos + 1);
                break;
            case 5:
                numMappings = atoi(currFieldPos);
                pairStat.numMappings[whichMate] = numMappings;
                if (numMappings > UINT16_MAX)
                    pairStat.numMappings[whichMate] = UINT16_MAX;
                break;
            case 6:
                cigarStr = currFieldPos;
                cigarLen = currFieldEnd - currFieldPos;
                break;
            case 7:
                mdStr = currFieldPos;
                mdLen = currFieldEnd - currFieldPos;
            default:
                break;
        }

        currFieldPos = currFieldEnd + 1;

    }while(*currFieldEnd != '>');

    if (cigarStr != NULL && mdStr != NULL)
    {
        pairStat.numMM[whichMate] = GetNumMismatchFromZA(&(pairStat.end[whichMate]), pairStat.softSize[whichMate], cigarStr, cigarLen, mdStr, mdLen);
        pairStat.end[whichMate] = pairStat.end[whichMate] + pAlignment->MatePosition - 1;
    }
    else
    {
        pairStat.numMM[whichMate] = GetNumMismatchFromBam(pairStat.softSize[whichMate]);
        pairStat.end[whichMate] = pAlignment->GetEndPosition(false, true);
    }
}

int BamPairTable::GetNumMismatchFromBam(int softSize[2])
{
    int numMM = 0;
    unsigned int numCigarField = pAlignment->CigarData.size();
    for (unsigned i = 0; i != numCigarField; ++i)
    {
        if (pAlignment->CigarData[i].Type == 'I' || pAlignment->CigarData[i].Type == 'D')
            numMM += (pAlignment->CigarData[i].Length);

        if (pAlignment->CigarData[i].Type == 'S')
        {
            if (i == 0)
                softSize[1] = pAlignment->CigarData[i].Length;
            else
                softSize[0] = pAlignment->CigarData[i].Length;
        }
    }

    string mdStr;
    if (pAlignment->GetTag("MD", mdStr))
    {
        const char* mdFieldPos = mdStr.c_str();
        while (mdFieldPos != NULL && *mdFieldPos != '\0')
        {
            if (isdigit(*mdFieldPos))
            {
                ++mdFieldPos;
                continue;
            }

            const char* mdFieldEnd = mdFieldPos + 1;
            while (!isdigit(*mdFieldEnd) && *mdFieldEnd != '\0')
                ++mdFieldEnd;

            if (*mdFieldPos != '^')
                numMM += mdFieldEnd - mdFieldPos;

            mdFieldPos = mdFieldEnd;
        }
    }

    return numMM;
}

int BamPairTable::GetNumMismatchFromZA(int32_t* pLen, int softSize[2], const char* cigarStr, unsigned int cigarLen, const char* mdStr, unsigned int mdLen)
{
    int numMM = 0;
    int len = 0;

    if (cigarStr == NULL)
        return numMM;

    for (const char* cigarEnd = strpbrk(cigarStr, "DIMS"); cigarEnd != NULL && cigarEnd - cigarStr < cigarLen; cigarEnd = strpbrk(cigarEnd + 1, "DIMS"))
    {
        const char* currPos = cigarEnd - 1;
        while (isdigit(*currPos))
            --currPos;

        int size = atoi(currPos + 1);

        if (*cigarEnd == 'S')
        {
            if (cigarEnd == cigarStr + cigarLen - 1)
                softSize[0] = size;
            else
                softSize[1] = size;
        }
        else if (*cigarEnd == 'M')
        {
            len += size;
        }
        else
        {
            numMM += size;
            
            if (*cigarEnd == 'D')
                len += size;
        }
    }

    const char* mdFieldPos = mdStr;
    while (mdFieldPos != NULL && mdFieldPos - mdStr < mdLen)
    {
        if (isdigit(*mdFieldPos))
        {
            ++mdFieldPos;
            continue;
        }

        const char* mdFieldEnd = mdFieldPos + 1;
        while (!isdigit(*mdFieldEnd) && *mdFieldEnd != '\0')
            ++mdFieldEnd;

        if (*mdFieldPos != '^')
            numMM += mdFieldEnd - mdFieldPos;

        mdFieldPos = mdFieldEnd;
    }

    *pLen = len;
    
    return numMM;
}

PairOrient BamPairTable::GetPairOrient(bool isUpMate)
{
    int8_t upMode = 0;
    int8_t downMode = 0;

    if (pAlignment->IsReverseStrand())
        upMode |= 1;

    if (pAlignment->IsMateReverseStrand())
        downMode |= 1;

    if (pAlignment->IsFirstMate() && pAlignment->IsSecondMate())
        return TGM_BAD_PAIR_MODE;
    else if (pAlignment->IsSecondMate())
        upMode |= (1 << 1);
    else if (pAlignment->IsFirstMate())
        downMode |= (1 << 1);
    else
        return TGM_BAD_PAIR_MODE;

    if (!isUpMate)
        TGM_SWAP(upMode, downMode, unsigned int);

    return (PairOrient) TGM_PairModeMap[((upMode << 2) | downMode)];
}

PairType BamPairTable::GetPairType(bool isUpMate)
{
    if (pairStat.spRef[0][0] != ' ' && pairStat.spRef[1][0] != ' ')
        return PT_UNKNOWN;
    else if (pairStat.spRef[0][0] == ' ' && pairStat.spRef[1][0] != ' ')
    {
        if (pairStat.bestMQ[0] > minSpMQ)
        {
            switch (pairStat.orient)
            {
                case TGM_1F2F:
                case TGM_1F2R:
                case TGM_2F1F:
                case TGM_2F1R:
                    return PT_SPECIAL3;
                    break;
                case TGM_1R2F:
                case TGM_1R2R:
                case TGM_2R1F:
                case TGM_2R1R:
                    return PT_SPECIAL5;
                    break;
                default:
                    return PT_UNKNOWN;
                    break;
            }
        }
        else
            return PT_UNKNOWN;
    }
    else if (pairStat.spRef[0][0] != ' ' && pairStat.spRef[1][0] == ' ')
    {
        if (pairStat.bestMQ[1] >= minSpMQ)
        {
            switch(pairStat.orient)
            {
                case TGM_2F1F:
                case TGM_2R1F:
                case TGM_1F2F:
                case TGM_1R2F:
                    return PT_SPECIAL3;
                    break;
                case TGM_2F1R:
                case TGM_2R1R:
                case TGM_1F2R:
                case TGM_1R2R:
                    return PT_SPECIAL5;
                    break;
                default:
                    return PT_UNKNOWN;
                    break;
            }
        }
        else
            return PT_UNKNOWN;
    }

    if (pairStat.fragLen == -1)
        return PT_CROSS;

    int32_t fragLenHigh = libTable.GetFragLenHigh(pairStat.readGrpID);
    int32_t fragLenLow = libTable.GetFragLenLow(pairStat.readGrpID);

    PairType type1 = (PairType) PairTypeMap[0][pairStat.orient];
    PairType type2 = (PairType) PairTypeMap[1][pairStat.orient];

    if (type1 == PT_NORMAL || type2 == PT_NORMAL)
    {
        if (pairStat.fragLen > fragLenHigh)
            return PT_LONG;
        else if (pairStat.fragLen < fragLenLow)
            return PT_SHORT;
        else
        {
            unsigned int upSoftSize = pairStat.softSize[0][0] + pairStat.softSize[0][1];
            unsigned int downSoftSize = pairStat.softSize[1][0] + pairStat.softSize[1][1];

            if (upSoftSize < minSoftSize && downSoftSize < minSoftSize)
            {
                return PT_NORMAL;
            }
            else if ( upSoftSize >= minSoftSize && downSoftSize >= minSoftSize)
            {
                return PT_UNKNOWN;
            }
            else
            {
                // we need to skip the anchor mate in a soft pair
                if ((isUpMate && downSoftSize >= minSoftSize) || (!isUpMate && upSoftSize >= minSoftSize))
                    return PT_UNKNOWN;

                if (pairStat.softSize[0][0] >= minSoftSize || pairStat.softSize[1][0] >= minSoftSize)
                    return PT_SOFT3;
                else if (pairStat.softSize[0][1] >= minSoftSize || pairStat.softSize[1][1] >= minSoftSize)
                    return PT_SOFT5;
                else
                    return PT_UNKNOWN;
            }
        }
    }
    else if ((type1 == PT_UNKNOWN && type2 == PT_UNKNOWN) || (type1 != PT_UNKNOWN && type2 != PT_UNKNOWN))
    {
        return PT_UNKNOWN;
    }
    else
    {
        return (type1 == PT_UNKNOWN ? type2 : type1);
    }
}

void BamPairTable::UpdateOrphanPair(void)
{
    /*
    khash_t(name)* pHash = (khash_t(name)*) readNameHash;
    khiter_t khIter = 0;
    int ret = 0;

    khIter = kh_get(name, pHash, pAlignment->Name.c_str());
    if (khIter != kh_end(pHash))
    {
        objIndex = kh_value(pHash, khIter);
        readNames.Return(objIndex.poolIdx);
        kh_del(name, pHash, khIter);
    }
    else
    {
        if (orphanPairs.size() <= numOrphanPairs)
            orphanPairs.resize(2 * numOrphanPairs);

        objIndex.idx = numOrphanPairs;
        ++numOrphanPairs;

        string& readName = readNames.Alloc(objIndex.poolIdx);
        readName.assign(pAlignment->Name);
        khIter = kh_put(name, pHash, readName.c_str(), &ret);

        kh_value(pHash, khIter) = objIndex;
    }
    */

    if (!pAlignment->IsMapped())
    {
        if (!pAlignment->GetTag("RG", readGrpName))
            return;

        uint32_t readGrpID = 0;
        if (!libTable.GetReadGrpID(readGrpID, readGrpName.c_str()))
            return;

        bool isUpMate = false;
        if (pAlignment->GetTag("ZA", zaStr))
            ParseZAstr(isUpMate);
        else
            return;

        if (pairStat.bestMQ[0] < minMQ)
            return;

        if (orphanPairs.IsFull())
        {
            orphanPairs.Resize(orphanPairs.Size() * 2);
            orphanPairs.InitToEnd();
        }

        OrphanPair& newOrphanPair = orphanPairs.End();

        newOrphanPair.bestMQ = pairStat.bestMQ[0];
        newOrphanPair.readGrpID = readGrpID;

        newOrphanPair.anchorPos = pAlignment->Position;
        newOrphanPair.anchorEnd = pairStat.end[0];
        newOrphanPair.refID = pAlignment->RefID;

        if (pAlignment->IsSecondMate())
        {
            if (pAlignment->IsMateReverseStrand())
                newOrphanPair.aOrient = TGM_1R;
            else
                newOrphanPair.aOrient = TGM_1F;
        }
        else
        {
            if (pAlignment->IsMateReverseStrand())
                newOrphanPair.aOrient = TGM_2R;
            else
                newOrphanPair.aOrient = TGM_2F;
        }

        TGM_SeqCpy(&(newOrphanPair.read), pAlignment->QueryBases.c_str(), pAlignment->QueryBases.size());

        orphanPairs.Increment();
    }
}

void BamPairTable::UpdateSpecialPair(void)
{
    // check which mate is the anchor
    int anchorIndex = 0;
    int specialIndex = 0;

    if (pairStat.spRef[0][0] == ' ')
    {
        anchorIndex = 0;
        specialIndex = 1;
    }
    else
    {
        anchorIndex = 1;
        specialIndex = 0;
    }

    if (specialPairs.IsFull())
        specialPairs.Resize(specialPairs.Size() * 2);

    SpecialPair& newSpecialPair = specialPairs.End();

    if (pAlignment->RefID == pAlignment->MateRefID)
    {
        int32_t fragLenHigh = libTable.GetFragLenHigh(pairStat.readGrpID);
        int32_t fragLenLow = libTable.GetFragLenLow(pairStat.readGrpID);
        int32_t fragLenMedian = libTable.GetFragLenMedian(pairStat.readGrpID);

        // widen the fragment length window to clean up the clusters
        int32_t spFragLenHigh = fragLenHigh + (fragLenHigh - fragLenMedian);
        int32_t spFragLenLow = fragLenLow - (fragLenMedian - fragLenLow);

        int32_t mapDist = 0;
        if (pairStat.readPairType == PT_SPECIAL3)                               // anchor is on the forward strand
            mapDist = pairStat.end[specialIndex] - pAlignment->Position;
        else                                                                    // anchor is on the reversed strand
            mapDist = pairStat.end[anchorIndex] - pAlignment->MatePosition;

        // we have to filter out those special pairs with normal fragment length
        if (mapDist >= spFragLenLow && mapDist <= spFragLenHigh)
            return;
    }

    newSpecialPair.readGrpID = pairStat.readGrpID;
    newSpecialPair.refID[0] = pAlignment->RefID;
    newSpecialPair.refID[1] = pAlignment->MateRefID;

    newSpecialPair.pos[0] = pAlignment->Position;
    newSpecialPair.pos[1] = pAlignment->MatePosition;

    newSpecialPair.end[0] = pairStat.end[anchorIndex];
    newSpecialPair.end[1] = pairStat.end[specialIndex];

    newSpecialPair.numMM[0] = pairStat.numMM[anchorIndex];
    newSpecialPair.numMM[1] = pairStat.numMM[specialIndex];

    newSpecialPair.bestMQ[0] = pairStat.bestMQ[anchorIndex];
    newSpecialPair.bestMQ[1] = pairStat.bestMQ[specialIndex];

    //  get the special reference ID
    uint32_t spRefID = 0;
    libTable.GetSpecialRefID(spRefID, pairStat.spRef[specialIndex]);

    newSpecialPair.numSpeicalHits = pairStat.numMappings[specialIndex];
    newSpecialPair.orient = pairStat.orient;
    newSpecialPair.pairType = pairStat.readPairType;
    newSpecialPair.specialID = spRefID;

    if (pAlignment->IsReverseStrand())
        newSpecialPair.aStrand = 1;
    else
        newSpecialPair.aStrand = 0;

    if (pAlignment->IsMateReverseStrand())
        newSpecialPair.sStrand = 1;
    else
        newSpecialPair.sStrand = 0;

    specialPairs.Increment();
}


void BamPairTable::UpdateLocalPair(Array<LocalPair>& localPairs)
{
    if (pAlignment->Position > pAlignment->MatePosition)
        return;

    if (pairStat.bestMQ[0] < minMQ || pairStat.bestMQ[1] < minMQ)
        return;

    if (localPairs.IsFull())
        localPairs.Resize(localPairs.Size() * 2);

    LocalPair& newLocalPair = localPairs.End();

    newLocalPair.readGrpID = pairStat.readGrpID;
    newLocalPair.refID = pAlignment->RefID;

    newLocalPair.pos[0] = pAlignment->Position;
    newLocalPair.pos[1] = pAlignment->MatePosition;

    newLocalPair.end[0] = pairStat.end[0];
    newLocalPair.end[1] = pairStat.end[1];

    newLocalPair.numMM[0] = pairStat.numMM[0];
    newLocalPair.numMM[1] = pairStat.numMM[1];

    newLocalPair.bestMQ[0] = pairStat.bestMQ[0];
    newLocalPair.bestMQ[1] = pairStat.bestMQ[1];

    uint32_t fragLenMedian = libTable.GetFragLenMedian(pairStat.readGrpID);

    if (pairStat.readPairType != PT_REVERSED)
        newLocalPair.fragLenQual = fragLenTable.GetQuality(pairStat.readGrpID, pairStat.fragLen, fragLenMedian);
    else
        newLocalPair.fragLenQual = INVALID_FRAG_LEN_QUAL;

    newLocalPair.readPairType = PT_LONG;
    newLocalPair.orient = pairStat.orient;

    localPairs.Increment();
}

void BamPairTable::UpdateInvertedPair(void)
{
    if (pAlignment->Position > pAlignment->MatePosition)
        return;

    if (pairStat.bestMQ[0] < minMQ || pairStat.bestMQ[1] < minMQ)
        return;

    if (invertedPairs.IsFull())
        invertedPairs.Resize(invertedPairs.Size() * 2);

    LocalPair& newInvertedPair = invertedPairs.End();

    newInvertedPair.readGrpID = pairStat.readGrpID;
    newInvertedPair.refID = pAlignment->RefID;

    if (pairStat.readPairType == PT_INVERTED3)
    {
        newInvertedPair.pos[0] = pAlignment->Position;
        newInvertedPair.pos[1] = pAlignment->MatePosition;

        newInvertedPair.end[0] = pairStat.end[0];
        newInvertedPair.end[1] = pairStat.end[1];

        newInvertedPair.numMM[0] = pairStat.numMM[0];
        newInvertedPair.numMM[1] = pairStat.numMM[1];

        newInvertedPair.bestMQ[0] = pairStat.bestMQ[0];
        newInvertedPair.bestMQ[1] = pairStat.bestMQ[1];

        ++numInverted3;
    }
    else
    {
        newInvertedPair.pos[0] = pAlignment->MatePosition;
        newInvertedPair.pos[1] = pAlignment->Position;

        newInvertedPair.end[0] = pairStat.end[1];
        newInvertedPair.end[1] = pairStat.end[0];

        newInvertedPair.numMM[0] = pairStat.numMM[1];
        newInvertedPair.numMM[1] = pairStat.numMM[0];

        newInvertedPair.bestMQ[0] = pairStat.bestMQ[1];
        newInvertedPair.bestMQ[1] = pairStat.bestMQ[0];
    }

    newInvertedPair.fragLenQual = INVALID_FRAG_LEN_QUAL;

    newInvertedPair.readPairType = pairStat.readPairType;
    newInvertedPair.orient = pairStat.orient;

    invertedPairs.Increment();
}

void BamPairTable::UpdateCrossPair(void)
{

}

void BamPairTable::UpdateSoftPair(void)
{
    if (softPairs.IsFull())
    {
        softPairs.Resize(softPairs.Size() * 2);
        softPairs.InitToEnd();
    }

    SoftPair& newSoftPair = softPairs.End();

    newSoftPair.readGrpID = pairStat.readGrpID;
    newSoftPair.refID = pAlignment->RefID;

    newSoftPair.pos = pAlignment->Position;
    newSoftPair.matePos = pAlignment->MatePosition;

    if (newSoftPair.pos <= newSoftPair.matePos)
    {
        if (pairStat.bestMQ[1] < minMQ)
            return;

        newSoftPair.end = pairStat.end[0];
    }
    else
    {
        if (pairStat.bestMQ[0] < minMQ)
            return;

        newSoftPair.end = pairStat.end[1];
    }

    newSoftPair.readPairType = pairStat.readPairType;
    newSoftPair.orient = pairStat.orient;

    if (pAlignment->IsReverseStrand())
        newSoftPair.sStrand = 1;
    else
        newSoftPair.sStrand = 0;

    if (pAlignment->IsMateReverseStrand())
        newSoftPair.aStrand = 1;
    else
        newSoftPair.aStrand = 0;

    newSoftPair.cigarLen = pAlignment->CigarData.size();
    newSoftPair.cigar = TransferCigar(pAlignment->CigarData);

    TGM_SeqCpy(&(newSoftPair.read), pAlignment->QueryBases.c_str(), pAlignment->QueryBases.size());
    softPairs.Increment();
}

uint32_t* BamPairTable::TransferCigar(const std::vector<BamTools::CigarOp>& cigarData)
{
    uint32_t* cigar = (uint32_t*) malloc(sizeof(uint32_t) * cigarData.size());
    if (cigar == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for the cigar string.\n");

    for (unsigned int i = 0; i != cigarData.size(); ++i)
    {
        cigar[i] = cigarData[i].Length;
        cigar[i] = cigar[i] << BAM_CIGAR_SHIFT;

        switch(cigarData[i].Type)
        {
            case 'D':
                cigar[i] |= BAM_CDEL;
                break;
            case 'H':
                cigar[i] |= BAM_CHARD_CLIP;
                break;
            case 'I':
                cigar[i] |= BAM_CINS;
                break;
            case 'M':
                cigar[i] |= BAM_CMATCH;
                break;
            case 'N':
                cigar[i] |= BAM_CREF_SKIP;
                break;
            case 'P':
                cigar[i] |= BAM_CPAD;
                break;
            case 'S':
                cigar[i] |= BAM_CSOFT_CLIP;
                break;
            default:
                break;
        }
    }

    return cigar;
}
