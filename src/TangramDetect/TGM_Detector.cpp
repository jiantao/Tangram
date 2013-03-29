/*
 * =====================================================================================
 *
 *       Filename:  TGM_Detector.cpp
 *
 *    Description:  Detect SV using read pair method
 *
 *        Created:  05/14/2012 01:04:16 PM
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

#include <limits.h>

#include "TGM_Utilities.h"
#include "TGM_Detector.h"

using namespace Tangram;

Detector::Detector(const DetectPars& detectPars, const LibTable& libTable, const BamPairTable& bamPairTable)
                  : pairAttrbtTable(libTable, bamPairTable), detectPars(detectPars), libTable(libTable), bamPairTable(bamPairTable)
{

}

Detector::~Detector()
{
    unsigned int numSp = libTable.GetNumSpecialRef();
    for (unsigned int i = 0; i != numSp; ++i)
    {
        for (unsigned int j = 0; j != pSpecialEventsTable[i].Size(); ++j)
        {
            if (pSpecialEventsTable[i][j].numFrag[0] > 0)
                free(pSpecialEventsTable[i][j].origIndex[0]);

            if (pSpecialEventsTable[i][j].numFrag[1] > 0)
                free(pSpecialEventsTable[i][j].origIndex[1]);
        }
    }

    if (pSpecialEventsTable != NULL)
        delete [] pSpecialEventsTable;
}

void Detector::Init(void)
{
    int numSp = libTable.GetNumSpecialRef();
    if (numSp > 0)
    {
        pSpecialEventsTable = new (std::nothrow) Array<SpecialEvent>[numSp];
        if (pSpecialEventsTable == NULL)
            TGM_ErrQuit("ERROR: Not enough memory for the special events array.\n");
    }
    else
        pSpecialEventsTable = NULL;

    // initialize the pair attributes table
    pairAttrbtTable.Init();
}

void Detector::CallEvents(void)
{
    unsigned int detectSet = detectPars.detectSet;

    for (unsigned int i = SV_DELETION; i <= SV_INTER_CHR_TRNSLCTN; ++i)
    {
        unsigned int shiftSize = i - 1;
        switch(i)
        {
            case SV_DELETION:
                if ((detectSet & (1 << shiftSize)) != 0)
                    CallDeletion();
                break;
            case SV_TANDEM_DUP:
                if ((detectSet & (1 << shiftSize)) != 0)
                    CallTandemDup();
                break;
            case SV_INVERSION:
                if ((detectSet & (1 << shiftSize)) != 0)
                    CallInversion();
                break;
            case SV_SPECIAL:
                if ((detectSet & (1 << shiftSize)) != 0)
                    CallSpecial();
                break;
            case SV_INTER_CHR_TRNSLCTN:
                if ((detectSet & (1 << shiftSize)) != 0)
                    CallTranslocation();
                break;
            default:
                break;
        }
    }
}

void Detector::CallDeletion(void)
{

}

void Detector::CallTandemDup(void)
{

}

void Detector::CallInversion(void)
{
    pairAttrbtTable.MakeInversion();

    cluster3.Init(pairAttrbtTable.invertedAttrbts, detectPars.minClusterSize, DEFAULT_MIN_STD);
    cluster5.Init(pairAttrbtTable.invertedAttrbts + 1, detectPars.minClusterSize, DEFAULT_MIN_STD);

    cluster3.Make();
    cluster5.Make();

    MakeInversions();
    // MergeInversions();
}

void Detector::MakeInversions(void)
{
    const Array<LocalPair>& invPairs = bamPairTable.invertedPairs;

    const Array<ClusterElmnt>* pCluster3 = cluster3.GetCluterElmnts();
    const Array<ClusterElmnt>* pCluster5 = cluster5.GetCluterElmnts();

    const Array<unsigned int>& next3 = cluster3.GetNextArray();
    const Array<unsigned int>& next5 = cluster5.GetNextArray();

    unsigned int actualNum3 = cluster3.GetActualNumElmnts();
    unsigned int actualNum5 = cluster5.GetActualNumElmnts();
    unsigned int actualNum = actualNum3 + actualNum5;

    unsigned int numEvents3 = pCluster3->Size();
    unsigned int numEvents5 = pCluster5->Size();

    // merged events will be stored in this array
    // so we make the capacity of this array bigger
    invEvents[0].Init(actualNum);
    invEvents[1].Init(actualNum5);

    // loop over the 3' clusters (3' mate of the fragment is in the inverted region)
    for (unsigned int i = 0; i != numEvents3; ++i)
    {
        const ClusterElmnt* pClusterElmnt = pCluster3->GetPointer(i);
        
        // skip the bad cluster
        if (pClusterElmnt->numReadPair == 0)
            continue;

        unsigned int j = pClusterElmnt->startIndex;

        unsigned int posMin5 = UINT_MAX;
        unsigned int posMax5 = 0;
        unsigned int endMax5 = 0;

        unsigned int posMin3 = UINT_MAX;
        unsigned int posMax3 = 0;
        unsigned int endMax3 = 0;

        Inversion* pNewInv = &(invEvents[0].End());
        pNewInv->fragLenMax = 0;

        const LocalPair* pInvPair = NULL;

        // loop over the inverted pairs in the current cluster
        do
        {
            unsigned int origIndex = (*cluster3.pPairAttrbts)[j].origIndex;
            pInvPair = invPairs.GetPointer(origIndex);

            unsigned int fragLenHigh = libTable.GetFragLenHigh(pInvPair->readGrpID);

            if (fragLenHigh > pNewInv->fragLenMax)
                pNewInv->fragLenMax = fragLenHigh;

            if ((unsigned int) pInvPair->pos[0] < posMin5)
                posMin5 = pInvPair->pos[0];
            
            if ((unsigned int) pInvPair->pos[0] > posMax5)
                posMax5 = pInvPair->pos[0];

            if ((unsigned int) pInvPair->end[0] > endMax5)
                endMax5 = pInvPair->end[0];

            if ((unsigned int) pInvPair->pos[1] < posMin3)
                posMin3 = pInvPair->pos[1];

            if ((unsigned int) pInvPair->pos[1] > posMax3)
                posMax3 = pInvPair->pos[1];

            if ((unsigned int) pInvPair->end[1] > endMax3)
                endMax3 = pInvPair->end[1];

            j = next3[j];

        }while(j != pClusterElmnt->startIndex);

        int pos = endMax5 + 1;
        int len = endMax3 - endMax5;
        int posU = ((posMax5 - posMin5) + (posMax3 - posMin3)) / (2 * pClusterElmnt->numReadPair);

        // skip those events with negative length
        if (len < 1)
            continue;

        pNewInv->refID = pInvPair->refID;

        pNewInv->pos = pos;
        pNewInv->posU = posU;
        pNewInv->length = len;

        pNewInv->numFrag[0] = pClusterElmnt->numReadPair;
        pNewInv->numFrag[1] = 0;

        pNewInv->clusterID[0] = pClusterElmnt->startIndex;
        pNewInv->clusterID[1] = -1;

        pNewInv->pos5[0] = posMin5;
        pNewInv->pos5[1] = endMax5;

        pNewInv->pos3[0] = posMin3;
        pNewInv->pos3[1] = endMax3;

        // index pointing to the split-read signal
        pNewInv->splitIdx = -1;

        // merge the new event with the previous event if they overlap
        if (invEvents[0].Size() > 0 && IsInvOverlapped(invEvents[0].Last(), *pNewInv))
        {
            Inversion mergedInv;
            DoInversionMerge(mergedInv, &(invEvents[0].Last()), pNewInv);

            if (IsSignificantInv(mergedInv))
                invEvents[0].Last() = mergedInv;
        }
        else
        {
            if (IsSignificantInv(*pNewInv))
                invEvents[0].Increment();
        }
    }

    invEvents[0].Sort(CompareInversion);

    // loop over the 5' clusters (5' of the fragment is in the inverted region)
    for (unsigned int i = 0; i != numEvents5; ++i)
    {
        const ClusterElmnt* pClusterElmnt = pCluster5->GetPointer(i);
        
        // skip the bad cluster
        if (pClusterElmnt->numReadPair == 0)
            continue;

        unsigned int j = pClusterElmnt->startIndex;

        unsigned int posMin5 = UINT_MAX;
        unsigned int posMax5 = 0;
        unsigned int endMax5 = 0;

        unsigned int posMin3 = UINT_MAX;
        unsigned int posMax3 = 0;
        unsigned int endMax3 = 0;

        // the new event will be stored in the element after the last one
        Inversion* pNewInv = &(invEvents[1].End());
        pNewInv->fragLenMax = 0;

        const LocalPair* pInvPair = NULL;

        // loop over the inverted pairs in the current cluster
        do
        {
            unsigned int origIndex = (*cluster5.pPairAttrbts)[j].origIndex;
            pInvPair = invPairs.GetPointer(origIndex);

            unsigned int fragLenHigh = libTable.GetFragLenHigh(pInvPair->readGrpID);

            if (fragLenHigh > pNewInv->fragLenMax)
                pNewInv->fragLenMax = fragLenHigh;

            if ((unsigned int) pInvPair->pos[1] < posMin5)
                posMin5 = pInvPair->pos[1];

            if ((unsigned int) pInvPair->pos[1] > posMax5)
                posMax5 = pInvPair->pos[1];

            if ((unsigned int) pInvPair->end[1] > endMax5)
                endMax5 = pInvPair->end[1];

            if ((unsigned int) pInvPair->pos[0] < posMin3)
                posMin3 = pInvPair->pos[0];

            if ((unsigned int) pInvPair->pos[0] > posMax3)
                posMax3 = pInvPair->pos[0];

            if ((unsigned int) pInvPair->end[0] > endMax3)
                endMax3 = pInvPair->end[0];

            j = next5[j];

        }while(j != pClusterElmnt->startIndex);

        int pos = posMin5;
        int len = posMin3 - posMin5;

        // position uncertainty
        int posU = ((posMax5 - posMin5) + (posMax3 - posMin3)) / (2 * pClusterElmnt->numReadPair);

        // skip those events with negative length
        if (len < 1)
            continue;

        pNewInv->refID = pInvPair->refID;

        pNewInv->pos = pos;
        pNewInv->length = len;
        pNewInv->posU = posU;

        pNewInv->numFrag[0] = 0;
        pNewInv->numFrag[1] = pClusterElmnt->numReadPair;

        pNewInv->clusterID[0] = -1;
        pNewInv->clusterID[1] = pClusterElmnt->startIndex;

        pNewInv->pos5[0] = posMin5;
        pNewInv->pos5[1] = endMax5;

        pNewInv->pos3[0] = posMin3;
        pNewInv->pos3[1] = endMax3;

        // index pointing to the split-read signal
        pNewInv->splitIdx = -1;

        // merge the new event with the previous event if they overlap
        if (invEvents[1].Size() > 0 && IsInvOverlapped(invEvents[1].Last(), *pNewInv))
        {
            Inversion mergedInv;
            DoInversionMerge(mergedInv, &(invEvents[1].Last()), pNewInv);

            if (IsSignificantInv(mergedInv))
                invEvents[1].Last() = mergedInv;
        }
        else
        {
            if (IsSignificantInv(*pNewInv))
                invEvents[1].Increment();
        }
    }

    invEvents[1].Sort(CompareInversion);
}

bool Detector::IsInvOverlapped(const Inversion& prevInv, const Inversion& newInv)
{
    unsigned int prevEnd = prevInv.pos + prevInv.length;
    unsigned int posDiff = abs((int) prevInv.pos - (int) newInv.pos);

    if (prevEnd > newInv.pos && posDiff < newInv.fragLenMax)
        return true;

    return false;
}

bool Detector::IsSignificantInv(const Inversion& inv)
{
    if (inv.length < (unsigned int) detectPars.minEventLength)
        return false;

    return true;
}

void Detector::DoInversionMerge(Inversion& mergedInv, const Inversion* pHeadInv, const Inversion* pTailInv)
{
    if (pHeadInv->numFrag[0] == 0 && pTailInv->numFrag[0] > 0)
        TGM_SWAP(pHeadInv, pTailInv, const Inversion*);

    unsigned int headCount = pHeadInv->numFrag[0] + pHeadInv->numFrag[1];
    unsigned int tailCount = pTailInv->numFrag[0] + pTailInv->numFrag[1];
    unsigned int totalCount = headCount + tailCount;

    mergedInv = *pHeadInv;
    mergedInv.pos3[0] = (pHeadInv->pos3[0] < pTailInv->pos3[0] ? pHeadInv->pos3[0] : pTailInv->pos3[0]);
    mergedInv.pos3[1] = (pHeadInv->pos3[1] > pTailInv->pos3[1] ? pHeadInv->pos3[1] : pTailInv->pos3[1]);
    mergedInv.pos5[0] = (pHeadInv->pos5[0] < pTailInv->pos5[0] ? pHeadInv->pos5[0] : pTailInv->pos5[0]);
    mergedInv.pos5[1] = (pHeadInv->pos5[1] > pTailInv->pos5[1] ? pHeadInv->pos5[1] : pTailInv->pos5[1]);

    bool bothForward = (pHeadInv->numFrag[1] == 0) && (pTailInv->numFrag[1] == 0);
    bool bothReversed = (pHeadInv->numFrag[0] == 0) && (pTailInv->numFrag[0] == 0);

    if (bothForward)
    {
        mergedInv.pos = mergedInv.pos5[1] + 1;
        mergedInv.length = mergedInv.pos3[1] - mergedInv.pos5[1] - 1;

        cluster3.Connect(pHeadInv->clusterID[0], pTailInv->clusterID[0]);
        mergedInv.clusterID[0] = pHeadInv->clusterID[0];
    }
    else if (bothReversed)
    {
        mergedInv.pos = mergedInv.pos5[0];
        mergedInv.length = mergedInv.pos3[0] - mergedInv.pos5[0] - 1;

        cluster5.Connect(pHeadInv->clusterID[1], pTailInv->clusterID[1]);
        mergedInv.clusterID[1] = pHeadInv->clusterID[1];
    }
    else
    {
        mergedInv.pos5[0] = pHeadInv->pos5[0];
        mergedInv.pos5[1] = (pHeadInv->pos5[1] * headCount + pTailInv->pos5[0] * tailCount) / totalCount;

        mergedInv.pos3[0] = (pHeadInv->pos3[1] * headCount +  pTailInv->pos3[0] * tailCount) / totalCount;
        mergedInv.pos3[1] = pTailInv->pos5[1];

        mergedInv.pos = mergedInv.pos5[1];
        mergedInv.length = mergedInv.pos3[0] - mergedInv.pos + 1;

        if (pHeadInv->numFrag[0] > 0 && pTailInv->numFrag[0] > 0) 
        {
            cluster3.Connect(pHeadInv->clusterID[0], pTailInv->clusterID[0]);
            mergedInv.clusterID[0] = pHeadInv->clusterID[0];
        }
        else if (pHeadInv->numFrag[1] > 0 && pTailInv->numFrag[1] > 0) 
        {
            cluster5.Connect(pHeadInv->clusterID[1], pTailInv->clusterID[1]);
            mergedInv.clusterID[1] = pHeadInv->clusterID[1];
        }
        else
        {
            mergedInv.clusterID[1] = pTailInv->clusterID[1];
        }
    }

    mergedInv.posU = DoubleRoundToInt((double) totalCount / ((double) headCount / pHeadInv->posU + (double) tailCount / pTailInv->posU));
    mergedInv.numFrag[0] = pHeadInv->numFrag[0] + pTailInv->numFrag[0];
    mergedInv.numFrag[1] = pHeadInv->numFrag[1] + pTailInv->numFrag[1];

    mergedInv.fragLenMax = pHeadInv->fragLenMax > pTailInv->fragLenMax ? pHeadInv->fragLenMax : pTailInv->fragLenMax;
}

void Detector::MergeInversions(void)
{
    /*  
    unsigned int numInv5 = invEvents[0].Size();
    unsigned int numInv3 = invEvents[1].Size();

    for (unsigned int i = 0; i != numInv5; ++i)
    {
        Inversion& inv5 = invEvents[0][i];
    }
    */
}

void Detector::CallSpecial(void)
{
    if (!libTable.HasSpecialRefNames())
        TGM_ErrQuit("ERROR: No special information detected in the bam file.\n"
                            "MEI detection aborted.\n");

    pairAttrbtTable.MakeSpecial();

    for (unsigned int i = 0; i != pairAttrbtTable.specialSize; i += 2)
    {
        cluster3.Init(pairAttrbtTable.pSpecialAttrbts + i, detectPars.minClusterSize, DEFAULT_MIN_STD);
        cluster5.Init(pairAttrbtTable.pSpecialAttrbts + i + 1, detectPars.minClusterSize, DEFAULT_MIN_STD);

        cluster3.Make();
        cluster5.Make();

        unsigned int spRefID = i / 2;

        MakeSpecialEvents(pSpecialEventsTable[spRefID]);
        MergeSpecialEvents(pSpecialEventsTable[spRefID]);
        SetOrigIndices(pSpecialEventsTable[spRefID], spRefID);
    }
}

void Detector::CallTranslocation(void)
{

}


void Detector::MakeSpecialEvents(Array<SpecialEvent>& specialEvents)
{
    const Array<SpecialPair>& specialPairs = bamPairTable.specialPairs;

    const Array<ClusterElmnt>* pCluster3 = cluster3.GetCluterElmnts();
    const Array<ClusterElmnt>* pCluster5 = cluster5.GetCluterElmnts();

    const Array<unsigned int>& next3 = cluster3.GetNextArray();
    const Array<unsigned int>& next5 = cluster5.GetNextArray();

    unsigned int actualNum3 = cluster3.GetActualNumElmnts();
    unsigned int actualNum5 = cluster5.GetActualNumElmnts();
    unsigned int actualNum = actualNum3 + actualNum5;

    unsigned int numEvents3 = pCluster3->Size();
    unsigned int numEvents5 = pCluster5->Size();

    specialEvents.Init(actualNum);
    specialEvents.SetSize(actualNum);

    unsigned int spIdx = 0;

    for (unsigned int j = 0; j != numEvents3; ++j)
    {
        const ClusterElmnt* pClusterElmnt = pCluster3->GetPointer(j);

        if (pClusterElmnt->numReadPair == 0)
            continue;

        unsigned int i = pClusterElmnt->startIndex;

        unsigned int posMin5 = UINT_MAX;
        unsigned int endMax5 = 0;

        unsigned int posMin3 = UINT_MAX;
        unsigned int endMax3 = 0;

        SpecialEvent* pSpecialEvent = specialEvents.GetPointer(spIdx);
        pSpecialEvent->fragLenMax = 0;

        ++spIdx;

        const SpecialPair* pSpecialPair = NULL;

        do
        {
            unsigned int origIndex = (*cluster3.pPairAttrbts)[i].origIndex;
            pSpecialPair = specialPairs.GetPointer(origIndex);

            unsigned int fragLenMedian = libTable.GetFragLenMedian(pSpecialPair->readGrpID);
            unsigned int fragLenHigh = libTable.GetFragLenHigh(pSpecialPair->readGrpID);

            if (fragLenHigh > pSpecialEvent->fragLenMax)
                pSpecialEvent->fragLenMax = fragLenHigh;

            unsigned int posUniq = pSpecialPair->pos[0];
            unsigned int endUniq = pSpecialPair->end[0];

            unsigned int posMultiple = posUniq + fragLenMedian - (pSpecialPair->end[1] - pSpecialPair->pos[1] + 1);
            unsigned int endMultiple = posUniq + fragLenMedian - 1;

            if (posUniq < posMin5)
                posMin5 = posUniq;

            if (endUniq > endMax5)
                endMax5 = endUniq;

            if (posMultiple < posMin3)
                posMin3 = posMultiple;

            if (endMultiple > endMax3)
                endMax3 = endMultiple;

            i = next3[i];

        }while(i != pClusterElmnt->startIndex);

        pSpecialEvent->refID = pSpecialPair->refID[0];

        pSpecialEvent->pos = endMax5;
        pSpecialEvent->numFrag[0] = pClusterElmnt->numReadPair;
        pSpecialEvent->numFrag[1] = 0;
        pSpecialEvent->clusterID[0] = pClusterElmnt->startIndex;

        pSpecialEvent->length = 0;

        pSpecialEvent->pos5[0] = posMin5;
        pSpecialEvent->pos5[1] = endMax5;

        pSpecialEvent->pos3[0] = posMin3;
        pSpecialEvent->pos3[1] = endMax3;

        pSpecialEvent->posUncertainty = DoubleRoundToInt((double) ((endMax5 - posMin5) + (endMax3 - posMin3)) / (double) (2 * pClusterElmnt->numReadPair));
    }

    for (unsigned int j = 0; j != numEvents5; ++j)
    {
        const ClusterElmnt* pClusterElmnt = pCluster5->GetPointer(j);

        if (pClusterElmnt->numReadPair == 0)
            continue;

        unsigned int i = pClusterElmnt->startIndex;

        unsigned int posMin5 = UINT_MAX;
        unsigned int endMax5 = 0;

        unsigned int posMin3 = UINT_MAX;
        unsigned int endMax3 = 0;

        SpecialEvent* pSpecialEvent = specialEvents.GetPointer(spIdx);
        pSpecialEvent->fragLenMax = 0;

        ++spIdx;

        const SpecialPair* pSpecialPair = NULL;

        do
        {
            unsigned int origIndex = (*cluster5.pPairAttrbts)[i].origIndex;
            pSpecialPair = specialPairs.GetPointer(origIndex);

            unsigned int fragLenMedian = libTable.GetFragLenMedian(pSpecialPair->readGrpID);
            unsigned int fragLenHigh = libTable.GetFragLenHigh(pSpecialPair->readGrpID);

            if (fragLenHigh > pSpecialEvent->fragLenMax)
                pSpecialEvent->fragLenMax = fragLenHigh;

            unsigned int posUniq = pSpecialPair->pos[0];
            unsigned int endUniq = pSpecialPair->end[0];

            unsigned int posMultiple = endUniq + 1 - fragLenMedian;
            unsigned int endMultiple = endUniq + pSpecialPair->end[1] + 1 - fragLenMedian - pSpecialPair->pos[1];

            if (posMultiple < posMin5)
                posMin5 = posMultiple;

            if (endMultiple > endMax5)
                endMax5 = endMultiple;

            if (posUniq < posMin3)
                posMin3 = posUniq;

            if (endUniq > endMax3)
                endMax3 = endUniq;

            i = next5[i];

        }while(i != pClusterElmnt->startIndex);

        pSpecialEvent->refID = pSpecialPair->refID[0];

        pSpecialEvent->pos = posMin3;
        pSpecialEvent->numFrag[0] = 0;
        pSpecialEvent->numFrag[1] = pClusterElmnt->numReadPair;
        pSpecialEvent->clusterID[1] = pClusterElmnt->startIndex;

        pSpecialEvent->length = 0;

        pSpecialEvent->pos5[0] = posMin5;
        pSpecialEvent->pos5[1] = endMax5;

        pSpecialEvent->pos3[0] = posMin3;
        pSpecialEvent->pos3[1] = endMax3;

        pSpecialEvent->posUncertainty = DoubleRoundToInt((double) ((endMax5 - posMin5) + (endMax3 - posMin3)) / (double) (2 * pClusterElmnt->numReadPair));
    }

    specialEvents.Sort(CompareSpecialEvents);
}

void Detector::MergeSpecialEvents(Array<SpecialEvent>& specialEvents)
{
    unsigned int headIndex = 1;
    unsigned int tailIndex = 0;
    unsigned int newSize = specialEvents.Size();
    unsigned int numSpecialEvents = specialEvents.Size();

    while (headIndex < numSpecialEvents)
    {
        SpecialEvent mergedEvent;
        uint32_t fragLenMax = specialEvents[headIndex].fragLenMax;
        if (specialEvents[tailIndex].fragLenMax > fragLenMax)
            fragLenMax = specialEvents[tailIndex].fragLenMax;

        if (abs(specialEvents[headIndex].pos - specialEvents[tailIndex].pos) < fragLenMax)
        {
            SpecialEvent* pHeadEvent = specialEvents.GetPointer(headIndex);
            SpecialEvent* pTailEvent = specialEvents.GetPointer(tailIndex);

            DoSpecialMerge(&mergedEvent, pHeadEvent, pTailEvent);

            *pTailEvent = mergedEvent;

            --newSize;
        }
        else
        {
            ++tailIndex;
            if (headIndex != tailIndex)
            {
                SpecialEvent* pHeadEvent = specialEvents.GetPointer(headIndex);
                SpecialEvent* pTailEvent = specialEvents.GetPointer(tailIndex);
                *pTailEvent = *pHeadEvent;
            }
        }

        ++headIndex;
    }

    specialEvents.SetSize(newSize);
    specialEvents.Sort(CompareSpecialEvents);
}

void Detector::DoSpecialMerge(SpecialEvent* pMergedEvent, const SpecialEvent* pHeadEvent, const SpecialEvent* pTailEvent)
{
    if (pTailEvent->numFrag[0] == 0)
        TGM_SWAP(pHeadEvent, pTailEvent, const SpecialEvent*);

    *pMergedEvent = *pTailEvent;

    bool bothForward = (pHeadEvent->numFrag[0] > 0) && (pTailEvent->numFrag[0] > 0);
    bool bothReversed = (pHeadEvent->numFrag[1] > 0) && (pTailEvent->numFrag[1] > 0);

    if (bothForward)
    {
        pMergedEvent->pos5[0] = (pTailEvent->pos5[0] < pHeadEvent->pos5[0] ? pTailEvent->pos5[0] : pHeadEvent->pos5[0]);
        pMergedEvent->pos5[1] = (pTailEvent->pos5[1] > pHeadEvent->pos5[1] ? pTailEvent->pos5[1] : pHeadEvent->pos5[1]);
        cluster3.Connect(pHeadEvent->clusterID[0], pTailEvent->clusterID[0]);
        pMergedEvent->clusterID[0] = pHeadEvent->clusterID[0];
    }
    else if (bothReversed)
    {
        pMergedEvent->pos3[0] = (pTailEvent->pos3[0] < pHeadEvent->pos3[0] ? pTailEvent->pos3[0] : pHeadEvent->pos3[0]);
        pMergedEvent->pos3[1] = (pTailEvent->pos3[1] > pHeadEvent->pos3[1] ? pTailEvent->pos3[1] : pHeadEvent->pos3[1]);
        cluster5.Connect(pHeadEvent->clusterID[1], pTailEvent->clusterID[1]);
        pMergedEvent->clusterID[1] = pHeadEvent->clusterID[1];
    }
    else
    {
        pMergedEvent->pos3[0] = pHeadEvent->pos3[0];
        pMergedEvent->pos3[1] = pHeadEvent->pos3[1];
        pMergedEvent->clusterID[1] = pHeadEvent->clusterID[1];
    }

    int numFrag5 = pMergedEvent->numFrag[0] = pHeadEvent->numFrag[0] + pTailEvent->numFrag[0];
    pMergedEvent->numFrag[1] = pHeadEvent->numFrag[1] + pTailEvent->numFrag[1];

    int32_t posTail = (pTailEvent->numFrag[0] > 0 ? pTailEvent->pos5[1] : pTailEvent->pos3[0]);
    int32_t posHead = (pHeadEvent->numFrag[1] > 0 ? pHeadEvent->pos3[0] : pHeadEvent->pos5[1]);

    pMergedEvent->fragLenMax = (pHeadEvent->fragLenMax > pTailEvent->fragLenMax ? pHeadEvent->fragLenMax : pTailEvent->fragLenMax);

    pMergedEvent->pos = (posHead < posTail ? posHead : posTail);
    pMergedEvent->posUncertainty = DoubleRoundToInt((double) (pMergedEvent->pos5[1] - pMergedEvent->pos5[0]) / numFrag5);
}

void Detector::SetOrigIndices(Array<SpecialEvent>& specialEvents, unsigned int spRefID)
{
    SpecialEvent* pSpecialEvent = NULL;

    unsigned int size = specialEvents.Size();

    for (unsigned int i = 0; i != size; ++i)
    {
        pSpecialEvent = specialEvents.GetPointer(i);
        pSpecialEvent->familyID = spRefID;
        pSpecialEvent->sense = 0.0;
        pSpecialEvent->splitIdx = -1;

        if (pSpecialEvent->numFrag[0] > 0)
        {
            pSpecialEvent->origIndex[0] = (uint32_t*) malloc(pSpecialEvent->numFrag[0] * sizeof(uint32_t));
            if (pSpecialEvent->origIndex[0] == NULL)
                TGM_ErrQuit("ERROR: Not enough memory for the original index array.\n");

            int centerIndex = pSpecialEvent->clusterID[0];
            int memberIndex = centerIndex;

            unsigned int j = 0;

            do
            {
                unsigned int origIndex = (*cluster3.pPairAttrbts)[memberIndex].origIndex;
                pSpecialEvent->origIndex[0][j] = origIndex;

                const SpecialPair& specialPair = bamPairTable.specialPairs[origIndex];
                if (specialPair.sStrand == 1)
                    pSpecialEvent->sense += 1.0;
                else
                    pSpecialEvent->sense -= 1.0;

                memberIndex = cluster3.next[memberIndex];
                ++j;

            }while(memberIndex != centerIndex);
        }
        else
            pSpecialEvent->origIndex[0] = NULL;

        if (pSpecialEvent->numFrag[1] > 0)
        {
            pSpecialEvent->origIndex[1] = (uint32_t*) malloc(pSpecialEvent->numFrag[1] * sizeof(uint32_t));
            if (pSpecialEvent->origIndex[1] == NULL)
                TGM_ErrQuit("ERROR: Not enough memory for the original index array.\n");

            int centerIndex = pSpecialEvent->clusterID[1];
            int memberIndex = centerIndex;

            unsigned int j = 0;

            do
            {
                unsigned int origIndex = (*cluster5.pPairAttrbts)[memberIndex].origIndex;
                pSpecialEvent->origIndex[1][j] = origIndex;

                const SpecialPair& specialPair = bamPairTable.specialPairs[origIndex];
                if (specialPair.sStrand == 0)
                    pSpecialEvent->sense += 1.0;
                else
                    pSpecialEvent->sense -= 1.0;

                memberIndex = cluster5.next[memberIndex];
                ++j;

            }while(memberIndex != centerIndex);
        }
        else
            pSpecialEvent->origIndex[1] = NULL;
    }
}
