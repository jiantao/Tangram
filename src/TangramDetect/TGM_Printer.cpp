/*
 * =====================================================================================
 *
 *       Filename:  TGM_Printer.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  08/16/2012 03:35:16 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include "TGM_Printer.h"

using namespace Tangram;


Printer::Printer(const Detector* pDetector, const Aligner* pAligner, const Reference* pRef, const LibTable& libTable, const BamPairTable& bamPairTable)
                : pDetector(pDetector), pAligner(pAligner), pRef(pRef), libTable(libTable), bamPairTable(bamPairTable)
{
    InitFeatures();
}

Printer::~Printer()
{

}

void Printer::Print(void)
{
    for (unsigned int i = 0; i != NUM_SV_TYPES; ++i)
    {
        switch (i)
        {
            case SV_SPECIAL:
                PrintSpecial();
                break;
            default:
                break;
        }
    }
}

void Printer::PrintSpecial(void)
{
    unsigned int numSp = libTable.GetNumSpecialRef();

    for (unsigned int i = 0; i != numSp; ++i)
    {
        unsigned int numEvents = pDetector->pSpecialEventsTable[i].Size();
        const SpecialEvent* pRpSpecials = NULL;

        if (numEvents > 0)
            pRpSpecials = pDetector->pSpecialEventsTable[i].GetPointer(0);

        const SplitEvent* pSplitSpecials = NULL;
        unsigned int splitLen = 0;

        if (pAligner != NULL)
            pSplitSpecials = pAligner->GetSpecialStartFromZA(splitLen, i);

        InitPrintIdx(numEvents, splitLen);

        while (GetNextSpecial(pRpSpecials, pSplitSpecials))
        {
            sampleMap.clear();
            InitFeatures();

            if (printIdx.pRpSpecial != NULL)
                SetSampleInfoRpSpecial(*(printIdx.pRpSpecial));

            if (printIdx.pSplitEvent != NULL)
                SetSampleInfoSplit(*(printIdx.pSplitEvent));

            SetSampleString();
            SetSpecialFeatures(i);

            printf("chr%s\t%d\t%d\t%c\t%d\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%s\n", features.anchorName, features.pos, features.pos + features.len + 1, features.strand, 
                    features.rpFrag[0], features.rpFrag[1], features.splitFrag[0], features.splitFrag[1], features.spRefName, features.pos5[0], features.pos5[1], 
                    features.pos3[0], features.pos3[1], formatted.str().c_str());
        }
    }

    if (pAligner != NULL)
    {
        unsigned int numFamily = pRef->familyName.size();

        for (unsigned int i = 0; i != numFamily; ++i)
        {
            int zaID = pRef->familyToZA[i];
            if (zaID < 0)
            {
                unsigned int len = 0;
                const SplitEvent* pSplitSpecials = pAligner->GetSpecialStartFromFamily(len, i);

                for (unsigned int j = 0; j != len; ++j)
                {
                    sampleMap.clear();
                    InitFeatures();

                    SetSampleInfoSplit(pSplitSpecials[j]);

                    SetSampleString();
                    SetSpecialFeaturesFromSplit(pSplitSpecials[j]);

                    printf("chr%s\t%d\t%d\t%c\t%d\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%s\n", features.anchorName, features.pos, features.pos + features.len + 1, features.strand, 
                            features.rpFrag[0], features.rpFrag[1], features.splitFrag[0], features.splitFrag[1], features.spRefName, features.pos5[0], features.pos5[1], 
                            features.pos3[0], features.pos3[1], formatted.str().c_str());
                }
            }
        }
    }
}

bool Printer::GetNextSpecial(const SpecialEvent* pRpSpecials, const SplitEvent* pSplitSpecials)
{
    if (printIdx.rpIdx == printIdx.rpSize)
        printIdx.pRpSpecial = NULL;
    else
        printIdx.pRpSpecial = pRpSpecials + printIdx.rpIdx;

    if (printIdx.pRpSpecial != NULL)
    {
        printIdx.pRpSpecial = NULL;
        while (printIdx.rpIdx < printIdx.rpSize)
        {
            if (pRpSpecials[printIdx.rpIdx].splitIdx >= 0)
                ++(printIdx.rpIdx);
            else
            {
                printIdx.pRpSpecial = pRpSpecials + printIdx.rpIdx;
                break;
            }
        }
    }

    if (printIdx.splitIdx == printIdx.splitSize)
        printIdx.pSplitEvent = NULL;
    else
        printIdx.pSplitEvent = pSplitSpecials + printIdx.splitIdx;

    if (printIdx.pRpSpecial == NULL && printIdx.pSplitEvent == NULL)
        return false;
    if (printIdx.pRpSpecial != NULL && printIdx.pSplitEvent != NULL)
    {
        if (printIdx.pRpSpecial->pos <= (unsigned int) printIdx.pSplitEvent->pos)
            printIdx.pSplitEvent = NULL;
        else
            printIdx.pRpSpecial = NULL;
    }

    if (printIdx.pRpSpecial != NULL)
        ++(printIdx.rpIdx);

    if (printIdx.pSplitEvent != NULL)
    {
        int rpIdx = printIdx.pSplitEvent->rpIdx;
        if (rpIdx >= 0)
        {
            printIdx.pRpSpecial = pRpSpecials + rpIdx;
        }

        ++(printIdx.splitIdx);
    }

    return true;
}

void Printer::SetSampleInfoRpSpecial(const SpecialEvent& rpSpecial)
{
    for (unsigned int k = 0; k != rpSpecial.numFrag[0]; ++k)
    {
        const SpecialPair& specialPair = bamPairTable.specialPairs[rpSpecial.origIndex[0][k]];
        unsigned int readGrpID = specialPair.readGrpID;
        unsigned int sampleID = 0;

        if (libTable.GetSampleID(sampleID, readGrpID))
            ++sampleMap[sampleID];
    }

    for (unsigned int k = 0; k != rpSpecial.numFrag[1]; ++k)
    {
        const SpecialPair& specialPair = bamPairTable.specialPairs[rpSpecial.origIndex[1][k]];
        unsigned int readGrpID = specialPair.readGrpID;
        unsigned int sampleID = 0;

        if (libTable.GetSampleID(sampleID, readGrpID))
            ++sampleMap[sampleID];
    }
}

void Printer::SetSampleInfoSplit(const SplitEvent& splitEvent)
{
    for (unsigned int i = 0; i != splitEvent.size3; ++i)
    {
        if (splitEvent.first3[i].isMajor)
        {
            unsigned int readGrpID = 0;
            unsigned int sampleID = 0;
            unsigned int origIdx = splitEvent.first3[i].origIdx;

            if (splitEvent.first3[i].isSoft)
                readGrpID = bamPairTable.softPairs[origIdx].readGrpID;
            else
                readGrpID = bamPairTable.orphanPairs[origIdx].readGrpID;

            if (libTable.GetSampleID(sampleID, readGrpID))
                ++sampleMap[sampleID];

            ++(features.splitFrag[0]);
        }
    }

    for (unsigned int i = 0; i != splitEvent.size5; ++i)
    {
        if (splitEvent.first5[i].isMajor)
        {
            unsigned int readGrpID = 0;
            unsigned int sampleID = 0;
            unsigned int origIdx = splitEvent.first5[i].origIdx;

            if (splitEvent.first5[i].isSoft)
                readGrpID = bamPairTable.softPairs[origIdx].readGrpID;
            else
                readGrpID = bamPairTable.orphanPairs[origIdx].readGrpID;

            if (libTable.GetSampleID(sampleID, readGrpID))
                ++sampleMap[sampleID];

            ++(features.splitFrag[1]);
        }
    }
}

void Printer::SetSampleString(void)
{
    const Array<char*>& sampleNames = libTable.GetSampleNames();

    formatted.clear();
    formatted.str("");

    for (std::map<unsigned int, unsigned int>::const_iterator itor = sampleMap.begin(); itor != sampleMap.end(); ++itor)
    {
        unsigned int sampleID = itor->first;
        const char* name = sampleNames[sampleID];
        formatted << name << "_" << itor->second << "_";
    }
}

void Printer::SetSpecialFeatures(unsigned int zaSpRefID)
{
    if (printIdx.pSplitEvent != NULL)
    {
        const SplitEvent& splitEvent = *(printIdx.pSplitEvent);
        SetSpecialFeaturesFromSplit(splitEvent);
    }

    if (printIdx.pRpSpecial != NULL)
    {
        const SpecialEvent& rpSpecial = *(printIdx.pRpSpecial);
        SetSpecialFeaturesFromRp(rpSpecial, zaSpRefID);
    }
}

void Printer::SetSpecialFeaturesFromSplit(const SplitEvent& splitEvent)
{
    const Array<char*>* pSpecialRefs = libTable.GetSpecialRefNames();

    features.pos = splitEvent.pos;
    features.len = splitEvent.len;

    features.strand = splitEvent.strand == 0 ? '+' : '-';

    int zaID = pRef->familyToZA[splitEvent.pSpecialData->familyID];
    if (zaID < 0)
        features.spRefName = pRef->familyName[splitEvent.pSpecialData->familyID].c_str();
    else
        features.spRefName = (*pSpecialRefs)[zaID];

    features.anchorName = pRef->refHeader.names[splitEvent.refID];

    features.pos5[0] = splitEvent.pos5[0];
    features.pos5[1] = splitEvent.pos5[1];

    features.pos3[0] = splitEvent.pos3[0];
    features.pos3[1] = splitEvent.pos3[1];
}

void Printer::SetSpecialFeaturesFromRp(const SpecialEvent& rpSpecial, unsigned int zaSpRefID)
{
    const Array<char*>& anchorNames = libTable.GetAnchorNames();
    const Array<char*>* pSpecialRefs = libTable.GetSpecialRefNames();

    features.rpFrag[0] = rpSpecial.numFrag[0];
    features.rpFrag[1] = rpSpecial.numFrag[1];

    if (features.anchorName == NULL)
    {
        features.anchorName = anchorNames[rpSpecial.refID];
        features.spRefName = (*pSpecialRefs)[zaSpRefID];
        
        features.pos = rpSpecial.pos;
        features.len = rpSpecial.length;

        features.pos5[0] = rpSpecial.pos5[0];
        features.pos5[1] = rpSpecial.pos5[1];

        features.pos3[0] = rpSpecial.pos3[0];
        features.pos3[1] = rpSpecial.pos3[1];
    }
}
