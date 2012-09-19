/*
 * =====================================================================================
 *
 *       Filename:  TGM_Aligner.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/10/2012 08:20:59 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
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
#include <pthread.h>
#include <string>
#include <cstdlib>

#include "TGM_Sequence.h"
#include "TGM_Aligner.h"
#include "TGM_FirstMapThread.h"

using namespace Tangram;

static int ComparePartial(const void* a, const void* b)
{
    const PrtlAlgnmnt* pPartialOne = (const PrtlAlgnmnt*) a;
    const PrtlAlgnmnt* pPartialTwo = (const PrtlAlgnmnt*) b;

    if (pPartialOne->refPos < pPartialTwo->refPos)
        return -1;
    else if (pPartialOne->refPos > pPartialTwo->refPos)
        return 1;
    else
        return 0;
}

static int CompareSplitEvent(const void* a, const void* b)
{
    const SplitEvent* pSplitOne = (const SplitEvent*) a;
    const SplitEvent* pSplitTwo = (const SplitEvent*) b;

    if (pSplitOne->svType < pSplitTwo->svType)
        return -1;
    else if (pSplitOne->svType > pSplitTwo->svType)
        return 1;
    else if (pSplitOne->svType == SV_NORMAL)
        return 0;
    else
    {
        if (pSplitOne->svType != SV_SPECIAL)
        {
            if (pSplitOne->refID < pSplitTwo->refID)
                return -1;
            else if (pSplitOne->refID > pSplitTwo->refID)
                return 1;
            else
            {
                if (pSplitOne->pos < pSplitTwo->pos)
                    return -1;
                else if (pSplitOne->pos > pSplitTwo->pos)
                    return 1;
                else
                {
                    if (pSplitOne->len < pSplitTwo->len)
                        return -1;
                    else if (pSplitOne->len > pSplitTwo->len)
                        return 1;
                    else
                        return 0;
                }
            }
        }
        else
        {
            if (pSplitOne->pSpecialData->familyID < pSplitTwo->pSpecialData->familyID)
                return -1;
            else if (pSplitOne->pSpecialData->familyID > pSplitTwo->pSpecialData->familyID)
                return 1;
            else
            {
                if (pSplitOne->refID < pSplitTwo->refID)
                    return -1;
                else if (pSplitOne->refID > pSplitTwo->refID)
                    return 1;
                else
                {
                    if (pSplitOne->pos < pSplitTwo->pos)
                        return -1;
                    else if (pSplitOne->pos > pSplitTwo->pos)
                        return 1;
                    else
                    {
                        if (pSplitOne->len < pSplitTwo->len)
                            return -1;
                        else if (pSplitOne->len > pSplitTwo->len)
                            return 1;
                        else
                            return 0;
                    }
                }
            }
        }
    }
}

Aligner::Aligner(Detector& detector, BamPairTable& bamPairTable, const AlignerPars& alignerPars, const Reference& ref, const LibTable& libTable)
                 : detector(detector), bamPairTable(bamPairTable), pars(alignerPars), ref(ref), libTable(libTable)
{

}

Aligner::~Aligner()
{
    unsigned int eventSize = splitEvents.Size();
    for (unsigned int i = 0; i != eventSize; ++i)
    {
        unsigned int size3 = splitEvents[i].size3;
        for (unsigned int j = 0; j != size3; ++j)
        {
            if (!splitEvents[i].first3[j].isSoft)
                free(splitEvents[i].first3[j].cigar);

            free(splitEvents[i].second5[j].cigar);
        }

        if (size3 > 0)
        {
            free(splitEvents[i].first3);
            free(splitEvents[i].second5);
        }

        unsigned int size5 = splitEvents[i].size5;
        for (unsigned int j = 0; j != size5; ++j)
        {
            if (!splitEvents[i].first5[j].isSoft)
                free(splitEvents[i].first5[j].cigar);

            free(splitEvents[i].second3[j].cigar);
        }

        if (size5 > 0)
        {
            free(splitEvents[i].first5);
            free(splitEvents[i].second3);
        }

        free(splitEvents[i].pSpecialData);
    }
}

void Aligner::Map(void)
{
    FirstMap();

    SecondMap();

    Merge();
}

void Aligner::FirstMap(void)
{
    Array<OrphanPair>& orphanPairs = bamPairTable.orphanPairs;
    unsigned int orphanSize = orphanPairs.Size();
    unsigned int softSize = bamPairTable.softPairs.Size();
    unsigned int partialSize = orphanSize + softSize;

    Array<PrtlAlgnmnt> firstPartials;
    Array<RefRegion> refRegions;

    firstPartials.Init(partialSize);
    refRegions.Init(orphanSize);

    firstPartials.SetSize(partialSize);
    refRegions.SetSize(orphanSize);

    // make the thread joinable
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    FirstMapData* firstMapData = (FirstMapData*) malloc(pars.numThread * sizeof(FirstMapData));
    if (firstMapData == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for the first split mapping data.\n");

    unsigned int orphanStart = 0;
    unsigned int shareSize = orphanSize / pars.numThread;
    unsigned int shareRemainder = orphanSize % pars.numThread;

    FirstMapThread firstMapThread(pars, libTable, bamPairTable, ref);

    for (int i = 0; i != pars.numThread; ++i)
    {
        firstMapData[i].idx = i;
        firstMapData[i].pFirstMapThread = &firstMapThread;
        firstMapData[i].orphanStart = orphanStart;
        firstMapData[i].firstPartials = firstPartials.GetPointer(orphanStart);
        firstMapData[i].refRegions = refRegions.GetPointer(orphanStart);
        firstMapData[i].orphanPairs = orphanPairs.GetPointer(orphanStart);
        firstMapData[i].orphanSize = shareSize;

        if (shareRemainder != 0)
        {
            firstMapData[i].orphanSize += 1;
            --shareRemainder;
        }

        orphanStart += firstMapData[i].orphanSize;

        int ret = pthread_create(&(firstMapData[i].thread), &attr, &FirstMapThread::StartThread, (void*) &(firstMapData[i]));
        if (ret != 0)
            TGM_ErrQuit("ERROR: Unable to create threads.\n");
    }

    pthread_attr_destroy(&attr);
    for (int i = 0; i != pars.numThread; ++i)
    {
        void* status;
        int ret = pthread_join(firstMapData[i].thread, &status);
        if (ret != 0)
            TGM_ErrQuit("ERROR: Unable to join threads.\n");
    }

    free(firstMapData);

    unsigned int j = orphanSize;
    for (unsigned int i = 0; i != softSize; ++i, ++j)
        InsertSoft(firstPartials, bamPairTable.softPairs[i], i, j);

    firstPartials.Sort(ComparePartial);
    CleanFirstPartials(firstPartials);

    if (firstPartials.Size() == 0)
        return;

    MergeFirstPartials(firstPartials);

#ifdef DEBUG
    PrintSplitEvent();
#endif
}

void Aligner::InsertFirstPartial(Array<PrtlAlgnmnt>& firstPartials, const s_align* pAlignment, const RefRegion& refRegion, 
                                 unsigned int idx, bool isUpStream, PartialType partialType)
{
    if (firstPartials.IsFull())
        firstPartials.Resize(firstPartials.Size() * 2);

    PrtlAlgnmnt& partial = firstPartials.End();

    partial.refPos = pAlignment->ref_begin1 + refRegion.start;
    partial.refEnd = pAlignment->ref_end1 + refRegion.start;
    partial.readPos = pAlignment->read_begin1;
    partial.readEnd = pAlignment->read_end1;
    partial.cigar = pAlignment->cigar;
    partial.cigarLen = pAlignment->cigarLen;

    partial.origIdx = idx;
    partial.isReversed = isUpStream ? 0 : 1;
    partial.isSoft = 0;
    partial.partialType = partialType;

    firstPartials.Increment();
}

void Aligner::InsertSoft(Array<PrtlAlgnmnt>& firstPartials, const SoftPair& softPair, unsigned int origIdx, unsigned int idx) const
{
    PrtlAlgnmnt& partial = firstPartials[idx];

    partial.refPos = softPair.pos;
    partial.refEnd = softPair.end;

    partial.readPos = 0;
    partial.readEnd = softPair.read.len - 1;
    if ((softPair.cigar[0] & BAM_CIGAR_MASK) == BAM_CSOFT_CLIP)
        partial.readPos += ((softPair.cigar[0] >> BAM_CIGAR_SHIFT));

    unsigned int cigarEnd = softPair.cigarLen - 1;
    if ((softPair.cigar[cigarEnd] & BAM_CIGAR_MASK) == BAM_CSOFT_CLIP)
        partial.readEnd -= ((softPair.cigar[cigarEnd] >> BAM_CIGAR_SHIFT));

    partial.cigar = softPair.cigar;
    partial.cigarLen = softPair.cigarLen ;

    partial.origIdx = origIdx;
    partial.isReversed = softPair.sStrand == 1 ? 1 : 0;
    partial.isSoft = 1;
    partial.partialType = softPair.readPairType == PT_SOFT3 ? PARTIAL_3 : PARTIAL_5;
}

void Aligner::CleanFirstPartials(Array<PrtlAlgnmnt>& firstPartials)
{
    int size = firstPartials.Size();

    for (int i = size - 1; i >= 0; --i)
    {
        if (firstPartials[i].cigar == NULL)
            --size;
        else
            break;
    }

    firstPartials.SetSize(size);
}

void Aligner::MergeFirstPartials(const Array<PrtlAlgnmnt>& firstPartials)
{
    SimpleRegion simpleRegion;

    unsigned int currStart = 0;
    unsigned int count3 = 0;
    unsigned int count5 = 0; 

    if (firstPartials[0].partialType == PARTIAL_3)
        ++count3;
    else
        ++count5;

    simpleRegion.pos = firstPartials[0].refPos;
    simpleRegion.end = firstPartials[0].refEnd;
    simpleRegion.count = 1;

    splitEvents.Init(DEFAULT_PARTIAL_SIZE);

    unsigned int firstSize = firstPartials.Size();
    for (unsigned int i = 1; i != firstSize; ++i)
    {
        if (firstPartials[i].refPos <= simpleRegion.end)
        {
            ++simpleRegion.count;
            if (firstPartials[i].refEnd > simpleRegion.end)
                simpleRegion.end = firstPartials[i].refEnd;

            if (firstPartials[i].partialType == PARTIAL_3)
                ++count3;
            else
                ++count5;
        }
        else
        {
            InsertNewSplitEvent(firstPartials, simpleRegion, count3, count5, currStart);

            count3 = 0;
            count5 = 0;

            if (firstPartials[i].partialType == PARTIAL_3)
                ++count3;
            else
                ++count5;

            simpleRegion.pos = firstPartials[i].refPos;
            simpleRegion.end = firstPartials[i].refEnd;
            simpleRegion.count = 1;

            currStart = i;
        }
    }

    InsertNewSplitEvent(firstPartials, simpleRegion, count3, count5, currStart);
}

void Aligner::InsertNewSplitEvent(const Array<PrtlAlgnmnt>& firstPartials, const SimpleRegion& simpleRegion, unsigned int count3, unsigned int count5, unsigned int currStart)
{
    if (splitEvents.IsFull())
        splitEvents.Resize(splitEvents.Size() * 2);

    SplitEvent& event = splitEvents.End();
    splitEvents.Increment();


    event.size3 = count3;
    event.size5 = count5;

    event.pos3[0] = INT32_MAX;
    event.pos3[1] = 0;

    event.pos5[0] = INT32_MAX;
    event.pos5[1] = 0;

    event.rpIdx = -1;

    event.refID = -1;
    event.pos = -1;
    event.len = 0;

    event.majorCount = 0;

    if (count3 > 0)
    {
        event.first3 = (PrtlAlgnmnt*) malloc(sizeof(PrtlAlgnmnt) * count3);
        if (event.first3 == NULL)
            TGM_ErrQuit("ERROR: Not enough memory for the first 3 prime partial alignments.\n");
    }

    if (count5 > 0)
    {
        event.first5 = (PrtlAlgnmnt*) malloc(sizeof(PrtlAlgnmnt) * count5);
        if (event.first5 == NULL)
            TGM_ErrQuit("ERROR: Not enough memory for the first 5 prime partial alignments.\n");
    }

    unsigned int idx3 = 0;
    unsigned int idx5 = 0;

    unsigned int currEnd = currStart + simpleRegion.count;
    for (unsigned int j = currStart; j != currEnd; ++j)
    {
        const OrphanPair& orphanPair = bamPairTable.orphanPairs[firstPartials[j].origIdx];
        int readLen = orphanPair.read.len;
        int alignedLen = firstPartials[j].readEnd - firstPartials[j].readPos + 1;

        if (firstPartials[j].partialType == PARTIAL_3)
        {
            event.first3[idx3] = firstPartials[j];
            // here the major partial alignments are those
            // have long enough unaligned sequence (>= trigger lenght)
            if (readLen - alignedLen >= pars.minTriggerLen)
            {
                event.first3[idx3].isMajor = 1;
                ++event.majorCount;
            }
            else
                event.first3[idx3].isMajor = 0;

            ++idx3;
        }
        else
        {
            event.first5[idx5] = firstPartials[j];
            // here the major partial alignments are those
            // have long enough unaligned sequence (>= trigger lenght)
            if (readLen - alignedLen >= pars.minTriggerLen)
            {
                event.first5[idx5].isMajor = 1;
                ++event.majorCount;
            }
            else
                event.first5[idx5].isMajor = 0;

            ++idx5;
        }
    }
}

#ifdef DEBUG

void Aligner::Print(const PrtlAlgnmnt& partialAlignment, const TGM_Sequence& read)
{
    int refPos = partialAlignment.refPos;
    int refEnd = partialAlignment.refEnd;
    int readPos = partialAlignment.readPos;
    int readEnd = partialAlignment.readEnd;
    printf("RefPos: %d\nRefEnd: %d\nReadLen: %lu\nReadPos: %d\nReadEnd: %d\n", refPos, refEnd, read.len, readPos, readEnd);

    for (int i = 0; i != partialAlignment.cigarLen; ++i)
    {
        int type = (partialAlignment.cigar[i] & BAM_CIGAR_MASK);
        int length = ((partialAlignment.cigar[i] >> BAM_CIGAR_SHIFT));

        switch(type)
        {
            case BAM_CMATCH:
            case BAM_CDEL:
            case BAM_CSOFT_CLIP:
                for (int j = 0; j != length; ++j)
                    printf("%c", char_table[ref.refSeq[refPos++]]);
                break;
            case BAM_CINS:
                for (int j = 0; j != length; ++j)
                    printf("_");
                break;
                break;
            default:
                break;
        }
    }

    printf("\n");

    for (int i = 0; i != partialAlignment.cigarLen; ++i)
    {
        int type = (partialAlignment.cigar[i] & BAM_CIGAR_MASK);
        int length = ((partialAlignment.cigar[i] >> BAM_CIGAR_SHIFT));

        switch(type)
        {
            case BAM_CDEL:
            case BAM_CSOFT_CLIP:
                for (int j = 0; j != length; ++j)
                    printf("_");
                break;
            case BAM_CMATCH:
            case BAM_CINS:
                for (int j = 0; j != length; ++j)
                    printf("%c", char_table[read.seq[readPos++]]);
                break;
                break;
            default:
                break;
        }
    }

    printf("\n\n");
}

void Aligner::PrintSplitEvent(void)
{
    for (unsigned int i = 0; i != splitEvents.Size(); ++i)
    {
        int refID = 0;
        int numFrag = splitEvents[i].size3 + splitEvents[i].size5;

        const PrtlAlgnmnt* pPartial = NULL;
        if (splitEvents[i].size3 > 0)
            pPartial = &(splitEvents[i].first3[0]);
        else
            pPartial = &(splitEvents[i].first5[0]);

        if (pPartial->isSoft)
            refID = bamPairTable.softPairs[pPartial->origIdx].refID;
        else
            refID = bamPairTable.orphanPairs[pPartial->origIdx].refID;

        // printf("chr%d\t%d\t%d\t%d\n", refID + 1, pPartial->refPos, pPartial->refEnd + 1, numFrag);
        printf("%d\t%d\t%d\n", splitEvents[i].refID, splitEvents[i].pos, splitEvents[i].len);
    }
}

#endif

void Aligner::SecondMap(void)
{
    SecondMapData mapData;
    SecondMapThread secondMapThread(splitEvents, pars, libTable, bamPairTable, ref);
    InitSecondMapData(mapData, secondMapThread);

    // make the thread joinable
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    for (int i = 0; i != pars.numThread; ++i)
    {
        int ret = pthread_create(&(mapData.pTags[i].thread), &attr, &SecondMapThread::StartThread, (void*) &(mapData));
        if (ret != 0)
            TGM_ErrQuit("ERROR: Unable to create threads.\n");
    }

    for (int i = 0; i != pars.numThread; ++i)
    {
        void* status;
        int ret = pthread_join(mapData.pTags[i].thread, &status);
        if (ret != 0)
            TGM_ErrQuit("ERROR: Unable to join threads.\n");
    }

    pthread_attr_destroy(&attr);
    DestroySecondMapData(mapData);
}

void Aligner::InitSecondMapData(SecondMapData& mapData, SecondMapThread& secondMapThread)
{
    mapData.currIdx = 0;
    mapData.totalWorks = splitEvents.Size();
    mapData.pSecondMap = &secondMapThread;

    mapData.pTags = (SecondMapTag*) malloc(pars.numThread * sizeof(SecondMapTag));
    if (mapData.pTags == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for the seoncd map tags.\n");

    for (int i = 0; i != pars.numThread; ++i)
        mapData.pTags[i].idx = i;

    int status = pthread_mutex_init(&(mapData.mutex), NULL);
    if (status != 0)
        TGM_ErrQuit("ERROR: Cannot initiate the mutex of second map data.\n");
}

void Aligner::DestroySecondMapData(SecondMapData& mapData)
{
    free(mapData.pTags);
    pthread_mutex_destroy(&(mapData.mutex));
}

void Aligner::Merge(void)
{
    CountSplitEvents();
    splitEvents.Sort(CompareSplitEvent);

    for (unsigned int i = SV_DELETION; i <= SV_INTER_CHR_TRNSLCTN; ++i)
    {
        switch (i)
        {
            case SV_SPECIAL:
                MergeSpecial();
                break;
            default:
                break;
        }
    }
}

void Aligner::CountSplitEvents(void)
{
    unsigned int numFamily = ref.familyName.size();

    svCount.Init(NUM_SV_TYPES);
    familyCount.Init(numFamily);

    svCount.SetSize(NUM_SV_TYPES);
    familyCount.SetSize(numFamily);

    unsigned int numSplit = splitEvents.Size();
    for (unsigned int i = 0; i != numSplit; ++i)
    {
        unsigned int svType = splitEvents[i].svType;
        ++(svCount[svType]);

        if (svType == SV_SPECIAL)
            ++(familyCount[splitEvents[i].pSpecialData->familyID]);
    }

    for (unsigned int i = 1; i < NUM_SV_TYPES; ++i)
        svCount[i] += svCount[i - 1];

    for (unsigned int i = 1; i < numFamily; ++i)
        familyCount[i] += familyCount[i - 1];
}

void Aligner::MergeSpecial(void)
{
    // unsigned int numSp = libTable.GetNumSpecialRef();
    unsigned int numFamily = ref.familyName.size();
    unsigned int start = svCount[SV_SPECIAL - 1];

    for (unsigned int i = 0; i != numFamily; ++i)
    {
        unsigned int spStart = i == 0 ? start : start + familyCount[i - 1];
        unsigned int spEnd = start + familyCount[i];
        
        int zaID = ref.familyToZA[i];
        if (zaID >= 0)
        {
            Array<SpecialEvent>& rpSpecials = detector.pSpecialEventsTable[zaID];
            if (rpSpecials.Size() == 0)
                continue;

            for (unsigned int j = spStart; j != spEnd; ++j)
            {
                DoMergeSpecial(rpSpecials, j);
            }
        }
    }
}

void Aligner::DoMergeSpecial(Array<SpecialEvent>& rpSpecials, unsigned int currSplitIdx)
{
    SplitEvent& splitEvent = splitEvents[currSplitIdx];
    SpecialEvent fakeSpecial;

    fakeSpecial.refID = splitEvent.refID;
    fakeSpecial.pos = splitEvent.pos;
    fakeSpecial.length = splitEvent.len;

    int rpSize = rpSpecials.Size();
    int rpIdx = rpSpecials.UpperBound(fakeSpecial, CompareSpecialEvents);
    if (rpIdx < 0) 
        rpIdx = rpSize - 1;


    int idx = rpIdx;

    int bestDist = -1;
    int bestRpIdx = -1;
    int lastBestDist = -1;

    while (idx >= 0)
    {
        int dist3 = abs(splitEvent.pos - rpSpecials[idx].pos3[0]);
        int dist5 = abs(splitEvent.pos - rpSpecials[idx].pos5[1]);

        int currBestDist = 0;
        int maxDist = rpSpecials[idx].fragLenMax;

        currBestDist = dist3 < dist5 ? dist3 : dist5;
        if (currBestDist > maxDist)
        {
            if (lastBestDist < 0 || currBestDist < lastBestDist)
                lastBestDist = currBestDist;
            else
                break;
        }
        else
        {
            if (bestDist < 0 || currBestDist < bestDist)
            {
                lastBestDist = currBestDist;
                bestDist = currBestDist;
                bestRpIdx = idx;
            }
            else
                break;
        }

        --idx;
    }

    idx = rpIdx + 1;
    lastBestDist = -1;
    while (idx < rpSize)
    {
        int dist3 = abs(splitEvent.pos - rpSpecials[idx].pos3[0]);
        int dist5 = abs(splitEvent.pos - rpSpecials[idx].pos5[1]);

        int currBestDist = 0;
        int maxDist = rpSpecials[idx].fragLenMax;

        currBestDist = dist3 < dist5 ? dist3 : dist5;
        if (currBestDist > maxDist)
            break;

        if (bestDist < 0 || currBestDist < bestDist)
        {
            bestDist = currBestDist;
            bestRpIdx = idx;
        }
        else
            break;

        ++idx;
    }

    if (bestDist < 0)
        return;

    int competitorIdx = rpSpecials[bestRpIdx].splitIdx;
    if (competitorIdx < 0)
    {
        rpSpecials[bestRpIdx].splitIdx = currSplitIdx;
        splitEvent.rpIdx = bestRpIdx;

        // printf("chr%s\t%d\t%d\t%d\t%d\n", ref.refHeader.names[splitEvent.refID], splitEvent.pos, splitEvent.svType, rpSpecials[bestRpIdx].pos5[1], rpSpecials[bestRpIdx].pos3[0]);
    }
    else
    {
        SplitEvent& competitor = splitEvents[competitorIdx];
        int compDist3 = abs(competitor.pos - rpSpecials[bestRpIdx].pos3[0]);
        int compDist5 = abs(competitor.pos - rpSpecials[bestRpIdx].pos5[1]);

        int bestCompDist = compDist3 < compDist5 ? compDist3 : compDist5;
        if (bestCompDist > bestDist)
        {
            competitor.rpIdx = -1;
            rpSpecials[bestRpIdx].splitIdx = currSplitIdx;
            splitEvent.rpIdx = bestRpIdx;
        }
    }
}
