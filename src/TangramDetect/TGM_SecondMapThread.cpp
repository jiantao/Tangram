/*
 * =====================================================================================
 *
 *       Filename:  TGM_SecondMapThread.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  08/01/2012 04:53:02 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include <string>
#include "TGM_Utilities.h"
#include "TGM_SecondMapThread.h"

using namespace std;
using namespace Tangram;

#ifdef DEBUG_1

static void CigarToString(std::string& cigarStr, const uint32_t* cigar, unsigned int cigarLen)
{
    cigarStr.clear();
    char buff[50];
    for (unsigned int i = 0; i != cigarLen; ++i)
    {
        sprintf(buff, "%d", (cigar[i] >> BAM_CIGAR_SHIFT));
        cigarStr += buff;
        switch ((cigar[i] & BAM_CIGAR_MASK))
        {
            case BAM_CINS:
                cigarStr += "I";
                break;
            case BAM_CDEL:
                cigarStr += "D";
                break;
            case BAM_CMATCH:
                cigarStr += "M";
                break;
            case BAM_CSOFT_CLIP:
                cigarStr += "S";
                break;
            default:
                cigarStr += "X";
                break;
        }
    }
}

#endif

SecondMapThread::SecondMapThread(Array<SplitEvent>& events, const AlignerPars& alignerPars, const LibTable& libInfoTable, const BamPairTable& pairTable, const Reference& reference)
                               : splitEvents(events), alignerPars(alignerPars), libTable(libInfoTable), bamPairTable(pairTable), ref(reference)
{

}

SecondMapThread::~SecondMapThread()
{

}

void* SecondMapThread::StartThread(void* threadData)
{
    SecondMapData* pMapData = (SecondMapData*) threadData;

    SecondMapThread& secondMap = *(pMapData->pSecondMap);

    Array<EventTry> eventTries;
    eventTries.Init(DEFAULT_SV_TYPE_NUM);

    // a holder for the minus strand sequence
    TGM_Sequence minusSeq = {NULL, 0, 0};
    RescuePartial rescuePartial(secondMap.alignerPars);

    vector< vector<unsigned int> > poll(secondMap.ref.familyName.size() * 2);

    while (true)
    {
        int status = pthread_mutex_lock(&(pMapData->mutex));
        if (status != 0)
            TGM_ErrQuit("ERROR: Cannot lock the mutex.\n");

        unsigned int i = pMapData->currIdx;
        ++(pMapData->currIdx);

        status = pthread_mutex_unlock(&(pMapData->mutex));
        if (status != 0)
            TGM_ErrQuit("ERROR: Cannot unlock the mutex.\n");

        if (i >= pMapData->totalWorks)
            break;

        SplitEvent& event = secondMap.splitEvents[i];
        event.pSpecialData = NULL;
        secondMap.SearchEventTries(eventTries, event);
        secondMap.InitSecondPartial(event);

        if (eventTries.Size() != 0)
        {
            for (unsigned int j = 0; j != eventTries.Size(); ++j)
            {
                bool isSucess = false;
                event.svType = SV_NORMAL;
                switch (eventTries[j].svType)
                {
                    case SV_SPECIAL:
                        isSucess = secondMap.TrySpecial(event, minusSeq, rescuePartial);
                        if (isSucess)
                        {
                            bool isGood = secondMap.ProcessSpecial(event, poll);
                            if (isGood)
                            {
                                event.svType = SV_SPECIAL;
                                /*  
                                printf("chr%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\n", secondMap.ref.refHeader.names[event.refID], event.pos, event.pos + event.len, event.strand,
                                        secondMap.ref.familyName[event.pSpecialData->familyID].c_str(), event.pSpecialData->pos, event.pSpecialData->end,
                                        event.pSpecialData->polyALen, event.pSpecialData->tsdLen);
                                */
                            }
                        }
                        break;
                    default:
                        break;
                }

                if (event.svType != SV_NORMAL)
                    break;
            }
        }
        else
        {

        }

        eventTries.Clear();
    }

    // clean up the memory
    TGM_SeqClean(&minusSeq);

    pthread_exit(NULL);
}

void SecondMapThread::InitSecondPartial(SplitEvent& splitEvent)
{
    if (splitEvent.size3 > 0)
    {
        splitEvent.second5 = (PrtlAlgnmnt*) calloc(sizeof(PrtlAlgnmnt),  splitEvent.size3);
        if (splitEvent.second5 == NULL)
            TGM_ErrQuit("ERROR: Not enough memory for the second 5 prime partial alignment.\n");
    }

    if (splitEvent.size5 > 0)
    {
        splitEvent.second3 = (PrtlAlgnmnt*) calloc(sizeof(PrtlAlgnmnt),  splitEvent.size5);
        if (splitEvent.second3 == NULL)
            TGM_ErrQuit("ERROR: Not enough memory for the second 5 prime partial alignment.\n");
    }
}

void SecondMapThread::SearchEventTries(Array<EventTry>& eventTries, const SplitEvent& splitEvent)
{
    // FIXME: try different SV candidates according to detected read-pair events
    if (splitEvent.size3 > 0 && splitEvent.size5 > 0)
    {
        eventTries[0].svType = SV_SPECIAL;
        eventTries.Increment();
    }
    else
    {
        eventTries[0].svType = SV_SPECIAL;
        eventTries.Increment();
    }
}

bool SecondMapThread::TrySpecial(SplitEvent& splitEvent, TGM_Sequence& minusSeq, RescuePartial& rescuePartial)
{
    if (splitEvent.majorCount == 0)
        return false;

    // increase the failLimit to gain more sensitivity
    unsigned int failLimit = DoubleRoundToInt((double) splitEvent.majorCount * 0.6);
    if (failLimit == splitEvent.majorCount)
        failLimit = splitEvent.majorCount - 1;

    unsigned int failCount = 0;

    // align the second half read to the entire special reference
    RefRegion refRegion = {ref.spRefSeq.GetPointer(0), ref.spRefSeq.Size(), 0};
    SecondFilter secondFilter = &SecondMapThread::SecondFilterSpecial;

    // true if we want to align the read with a different orientation from the first partial
    bool doOtherFirst = false;

    // true if the second partial is rescued
    bool isRescued = false;

    uint8_t polyALen = 0;

    for (unsigned int i = 0; i != splitEvent.size3; ++i)
    {
        const PrtlAlgnmnt& firstPartial = splitEvent.first3[i];
        uint8_t isReversed = false;

        s_align* pSecAlignment = AlignSecPartial(isRescued, rescuePartial, doOtherFirst, isReversed, polyALen, firstPartial, refRegion, minusSeq, secondFilter);
        if (pSecAlignment == NULL)
        {
            // failure count only increase when this is a major partial alignment (unaligned part is long enough)
            if (firstPartial.isMajor != 0)
                ++failCount;

            // quit if more than half of the first partial alingments cannot find their second partials
            if (failCount > failLimit)
            {
                CleanUpSecond(splitEvent);
                return false;
            }
        }
        else
        {
            UpdateSecPartial(splitEvent.second5[i], isRescued, rescuePartial, pSecAlignment, firstPartial, isReversed, polyALen, SV_SPECIAL);
            if (isRescued)
            {
                isRescued = false;
                rescuePartial.Clear();
                align_destroy(pSecAlignment);
            }
            else
                free(pSecAlignment);
        }
    }

    for (unsigned int i = 0; i != splitEvent.size5; ++i)
    {
        const PrtlAlgnmnt& firstPartial = splitEvent.first5[i];
        uint8_t isReversed = false;

        s_align* pSecAlignment = AlignSecPartial(isRescued, rescuePartial, doOtherFirst, isReversed, polyALen, firstPartial, refRegion, minusSeq, secondFilter);
        if (pSecAlignment == NULL)
        {
            // failure count only increase when this is a major partial alignment (unaligned part is long enough)
            if (firstPartial.isMajor != 0)
                ++failCount;

            if (failCount > failLimit)
            {
                CleanUpSecond(splitEvent);
                return false;
            }
        }
        else
        {
            UpdateSecPartial(splitEvent.second3[i], isRescued, rescuePartial, pSecAlignment, firstPartial, isReversed, polyALen, SV_SPECIAL);
            if (isRescued)
            {
                isRescued = false;
                rescuePartial.Clear();
                align_destroy(pSecAlignment);
            }
            else
                free(pSecAlignment);
        }
    }

    return true;
}

s_align* SecondMapThread::AlignSecPartial(bool& isRescued, RescuePartial& rescuePartial, bool& doOtherFirst, uint8_t& isReversed, uint8_t& polyALen,
                                          const PrtlAlgnmnt& firstPartial, const RefRegion& refRegion, TGM_Sequence& minusSeq, SecondFilter secondFilter)
{
    unsigned int idx = firstPartial.origIdx;

    isReversed = firstPartial.isReversed;

    const TGM_Sequence* pRead = NULL;
    const int8_t* readSeq = NULL;
    unsigned int readLen = 0;

    if (!firstPartial.isSoft)
        pRead = &(bamPairTable.orphanPairs[idx].read);
    else
        pRead = &(bamPairTable.softPairs[idx].read);

    readSeq = pRead->seq;
    readLen = pRead->len;

    if (doOtherFirst)
    {
        isReversed ^= 1;
        TGM_SeqDup(&minusSeq, pRead);
        TGM_SeqRevComp(&minusSeq);

        readSeq = minusSeq.seq;
    }

    // TODO: adjust the gap open and extention penalty
    s_profile* pProfile = ssw_init(readSeq, readLen, alignerPars.secMat, 5, 2);
    s_align* pAlignment = ssw_align(pProfile, refRegion.pRef, refRegion.len, alignerPars.secGapOpen, alignerPars.secGapExt, 
                                    alignerPars.flag, alignerPars.scoreFilter, alignerPars.distFilter, readLen);

    init_destroy(pProfile);

#ifdef DEBUG_1

    std::string cigarStr;
    std::string seqStr;

#endif

    bool passFilter = (this->*secondFilter)(isRescued, rescuePartial, polyALen, pAlignment, isReversed, firstPartial, refRegion, pRead->seq, pRead->len);
    if (passFilter)
    {

#ifdef DEBUG_1

            SeqToString(seqStr, readSeq, readLen);
            const char* strand = isReversed == 1 ? "minus" : "plus";
            const char* fStrand = firstPartial.isReversed == 1 ? "minus" : "plus";
            const char* rescue = isRescued ? "rescued" : "lost";

            int refPos = pAlignment->ref_begin1;
            int refEnd = pAlignment->ref_end1;
            int readPos = pAlignment->read_begin1;
            int readEnd = pAlignment->read_end1;
            int score = pAlignment->score1;
            const uint32_t* cigar = pAlignment->cigar;
            int cigarLen = pAlignment->cigarLen;

            if (isRescued)
            {
                refPos = rescuePartial.refPos;
                refEnd = rescuePartial.refEnd;
                readPos = rescuePartial.readPos;
                readEnd = rescuePartial.readEnd;
                score = rescuePartial.bestScore;

                cigar = rescuePartial.cigar;
                cigarLen = rescuePartial.cigarLen;
            }

            CigarToString(cigarStr, cigar, cigarLen);

            printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s", firstPartial.refPos, firstPartial.refEnd, refPos + refRegion.start, 
                    refEnd + refRegion.start, readPos, readEnd, score, strand, rescue, cigarStr.c_str(), seqStr.c_str());

            CigarToString(cigarStr, firstPartial.cigar, firstPartial.cigarLen);
            printf("\t%s\t%s\n", cigarStr.c_str(), fStrand);
#endif

        return pAlignment;
    }
    else
    {
        // second chance: the other orientation
        align_destroy(pAlignment);
        if (isRescued)
        {
            rescuePartial.Clean();
            isRescued = false;
        }

        pProfile = NULL;
        pAlignment = NULL;

        if (isReversed == firstPartial.isReversed)
        {
            TGM_SeqDup(&minusSeq, pRead);
            TGM_SeqRevComp(&minusSeq);

            readSeq = minusSeq.seq;
        }
        else
        {
            readSeq = pRead->seq;
        }

        // TODO: adjust the gap open and extention penalty
        pProfile = ssw_init(readSeq, readLen, alignerPars.secMat, 5, 2);
        pAlignment = ssw_align(pProfile, refRegion.pRef, refRegion.len, alignerPars.secGapOpen, alignerPars.secGapExt, 
                               alignerPars.flag, alignerPars.scoreFilter, alignerPars.distFilter, readLen);

        isReversed ^= 1;
        init_destroy(pProfile);

        passFilter = (this->*secondFilter)(isRescued, rescuePartial, polyALen, pAlignment, isReversed, firstPartial, refRegion, pRead->seq, pRead->len);
        if (passFilter)
        {
            if (isReversed != firstPartial.isReversed)
                doOtherFirst = true;
            else
                doOtherFirst = false;

            return pAlignment;
        }
        else
        {
            if (isRescued)
            {
                rescuePartial.Clean();
                isRescued = false;
            }

            align_destroy(pAlignment);
            return NULL;
        }
    }
}

bool SecondMapThread::SecondFilterSpecial(bool& isRescued, RescuePartial& rescuePartial, uint8_t& polyALen, const s_align* pAlignment, uint8_t isReversed, 
                                          const PrtlAlgnmnt& firstPartial, const RefRegion& refRegion, const int8_t* readSeq, int readLen)
{
    isRescued = false;

    // best alignment is not available
    if (pAlignment->ref_begin1 < 0 || pAlignment->read_begin1 < 0)
        return false;

    // length filter
    int alignedReadLen = pAlignment->read_end1 - pAlignment->read_begin1 + 1;
    if (alignedReadLen < alignerPars.minAlignedLen)
        return false;

    double scoreRate = (double) pAlignment->score1 / (alignerPars.mat[0] * alignedReadLen);
    PartialType partialType = GetPartialType(pAlignment->read_begin1, pAlignment->read_end1, readLen);

    if (scoreRate >= alignerPars.minScoreRate)
    {
        if (readLen - alignedReadLen < alignerPars.minAlignedLen)
            return false;
    }
    else
    {
        // rescue those low score alignments
        isRescued = rescuePartial.RescueLowScore(partialType, pAlignment, readSeq, readLen, refRegion);
        if (!isRescued)
            return false;
    }

    // partial type filter
    if (partialType == PARTIAL_UNKNOWN 
        || (isReversed == firstPartial.isReversed && partialType == firstPartial.partialType) 
        || (isReversed != firstPartial.isReversed && partialType != firstPartial.partialType))
    {
        return false;
    }

    partialType = (PartialType) (firstPartial.partialType ^ 1);


    // cover rate filter
    int lowPos = 0;
    int lowEnd = 0;
    int highPos = 0;
    int highEnd = 0;
    int midPos = 0;
    int midEnd = 0;
    int coverLen = 0;
    int actualPos = 0;
    int actualEnd = 0;
    int8_t tailChar = 0;

    if (isReversed == firstPartial.isReversed)
    {
        if (!isRescued)
        {
            actualPos = pAlignment->read_begin1;
            actualEnd = pAlignment->read_end1;
        }
        else
        {
            actualPos = rescuePartial.readPos;
            actualEnd = rescuePartial.readEnd;
        }

        // "A"
        tailChar = 0;
    }
    else
    {
        if (!isRescued)
        {
            actualPos = readLen - pAlignment->read_end1 - 1;
            actualEnd = readLen - pAlignment->read_begin1 - 1;
        }
        else
        {
            actualPos = readLen - rescuePartial.readEnd - 1;
            actualEnd = readLen - rescuePartial.readPos - 1;
        }

        // "T"
        tailChar = 3;
    }

    double entropy = SeqGetEntropy(readSeq + actualPos, actualEnd - actualPos + 1);
    if (entropy < alignerPars.minEntropy)
        return false;

    if (actualPos < firstPartial.readPos)
    {
        lowPos = actualPos;
        lowEnd = actualEnd;
        highPos = firstPartial.readPos;
        highEnd = firstPartial.readEnd;
    }
    else
    {
        lowPos = firstPartial.readPos;
        lowEnd = firstPartial.readEnd;
        highPos = actualPos;
        highEnd = actualEnd;
    }

    bool checkPolyA = true;
    if (highPos <= lowEnd)
    {
        coverLen = highEnd - lowPos + 1;
        checkPolyA = false;
    }
    else
    {
        coverLen = (lowEnd - lowPos + 1) + (highEnd - highPos + 1);

        if ((firstPartial.partialType == PARTIAL_3 && isReversed == firstPartial.isReversed)
            || (firstPartial.partialType == PARTIAL_5 && isReversed != firstPartial.isReversed))
        {
            checkPolyA = false;
        }
        else
        {
            int spRefID = 0;
            const char* spRefName = ref.GetSpRefName(spRefID, pAlignment->ref_begin1);
            if (spRefName == NULL)
                TGM_ErrQuit("ERROR: Special reference position is out of range.\n");

            // check polyA tail only for non-LTR retrotransposons
            if (strstr(spRefName, "ALU") == NULL && strstr(spRefName, "L1") == NULL && strstr(spRefName, "SVA") == NULL)
                checkPolyA = false;
        }
    }

    double coverRate = (double) coverLen /  readLen;
    polyALen = 0;

    if (coverRate < alignerPars.minCoverRate && !checkPolyA)
        return false;
    else if (coverRate >= alignerPars.minCoverRate && !checkPolyA)
        return true;
    else
    {
        midPos = lowEnd + 1;
        midEnd = highPos - 1;

        while(midPos < readLen && readSeq[midPos] != tailChar)
            ++midPos;

        while(midEnd >= 0 && readSeq[midEnd] != tailChar)
            --midEnd;

        int midLen = midEnd - midPos + 1;

        if (midLen < MIN_POLYA_SIZE)
            return false;

        int tailCount = 0;
        for (int i = midPos; i <= midEnd; ++i)
        {
            if (readSeq[i] == tailChar)
                ++tailCount;
        }

        double polyaRate = (double) tailCount / midLen;
        if (polyaRate < MIN_POLYA_RATE)
            return false;

        coverLen += midLen;
        coverRate = (double) coverLen /  readLen;
        if (coverRate < alignerPars.minCoverRate)
            return false;

        polyALen = midLen;
    }

    return true;
}


void SecondMapThread::UpdateSecPartial(PrtlAlgnmnt& secPartial, bool isRescued, const RescuePartial& rescuePartial, const s_align* pSecAlignment, 
                                       const PrtlAlgnmnt& firstPartial, uint8_t isReversed, uint8_t polyALen, SV_EventType svType)
{
    if (!isRescued)
    {
        secPartial.refPos = pSecAlignment->ref_begin1;
        secPartial.refEnd = pSecAlignment->ref_end1;

        secPartial.readPos = pSecAlignment->read_begin1;
        secPartial.readEnd = pSecAlignment->read_end1;

        secPartial.cigar = pSecAlignment->cigar;
        secPartial.cigarLen = pSecAlignment->cigarLen;
    }
    else
    {
        secPartial.refPos = rescuePartial.refPos;
        secPartial.refEnd = rescuePartial.refEnd;

        secPartial.readPos = rescuePartial.readPos;
        secPartial.readEnd = rescuePartial.readEnd;

        secPartial.cigar = rescuePartial.cigar;
        secPartial.cigarLen = rescuePartial.cigarLen;
    }

    secPartial.isSoft = firstPartial.isSoft;
    secPartial.isReversed = isReversed;

    secPartial.origIdx = firstPartial.origIdx;

    if (svType == SV_SPECIAL)
    {
        int spRefID = 0;
        ref.GetSpRefName(spRefID, secPartial.refPos);
        if (spRefID < 0)
            TGM_ErrQuit("ERROR: Special reference position is out of range.\n");

        secPartial.spRefID = spRefID;
    }
    else
        secPartial.spRefID = -1;

    secPartial.partialType = firstPartial.partialType ^ 1;
    secPartial.polyALen = polyALen;
}

void SecondMapThread::CleanUpSecond(SplitEvent& splitEvent)
{
    for (unsigned int i = 0; i != splitEvent.size3; ++i)
    {
        free(splitEvent.second5[i].cigar);
        splitEvent.second5[i].cigar = NULL;
    }

    for (unsigned int i = 0; i != splitEvent.size5; ++i)
    {
        free(splitEvent.second3[i].cigar);
        splitEvent.second3[i].cigar = NULL;
    }
}

bool SecondMapThread::ProcessSpecial(SplitEvent& splitEvent, vector< vector<unsigned int> >& poll)
{
    ClearPoll(poll);
    int validCount = 0;

    for (unsigned int i = 0; i != splitEvent.size3; ++i)
    {
        splitEvent.first3[i].isMajor = 0;
        splitEvent.second5[i].isMajor = 0;

        if (splitEvent.second5[i].cigar != NULL)
        {
            ++validCount;
            int familyID = ref.GetFamilyID(splitEvent.second5[i].spRefID);
            int strand = splitEvent.second5[i].isReversed == splitEvent.first3[i].isReversed ? 0 : 1;
            int vectorID = familyID * 2 + strand;
            poll[vectorID].push_back(i);
        }
    }

    for (unsigned int i = 0; i != splitEvent.size5; ++i)
    {
        splitEvent.first5[i].isMajor = 0;
        splitEvent.second3[i].isMajor = 0;

        if (splitEvent.second3[i].cigar != NULL)
        {
            ++validCount;
            int familyID = ref.GetFamilyID(splitEvent.second3[i].spRefID);
            int strand = splitEvent.second3[i].isReversed == splitEvent.first5[i].isReversed ? 0 : 1;
            int vectorID = familyID * 2 + strand;
            poll[vectorID].push_back(i + splitEvent.size3);
        }
    }

    unsigned int pollSize = poll.size();
    unsigned int maxCount = 0;
    unsigned int maxID = 0;
    for (unsigned int i = 0; i != pollSize; ++i)
    {
        if (poll[i].size() > maxCount)
        {
            maxCount = poll[i].size();
            maxID = i;
        }
    }

    double agreeRate = (double) maxCount / validCount;
    if (agreeRate < alignerPars.minAgreeRate)
        return false;

    bool hasBoth = false;
    if (poll[maxID][0] < splitEvent.size3 && poll[maxID][maxCount - 1] >= splitEvent.size3)
        hasBoth = true;

    splitEvent.pSpecialData = (SplitSpecial*) calloc(1, sizeof(SplitSpecial));
    if (splitEvent.pSpecialData == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for the split special data.\n");

    splitEvent.pSpecialData->familyID = maxID / 2;
    splitEvent.strand = maxID % 2;
    splitEvent.refID = ref.refID;
    
    PrtlAlgnmnt* firstPartials = splitEvent.first3;
    PrtlAlgnmnt* secondPartials = splitEvent.second5;

    splitEvent.pos5[0] = -1;
    splitEvent.pos5[1] = -1;
    splitEvent.pos3[0] = -1;
    splitEvent.pos3[1] = -1;

    int32_t* pStart = &(splitEvent.pos5[0]);
    int32_t* pEnd = &(splitEvent.pos5[1]);

    splitEvent.pSpecialData->pos = -1;
    splitEvent.pSpecialData->end = -1;
    splitEvent.pSpecialData->polyALen = 0;

    unsigned int selectedSize = poll[maxID].size();
    for (unsigned int i = 0; i != selectedSize; ++i)
    {
        unsigned int idx = poll[maxID][i];
        if (idx >= splitEvent.size3)
        {
            idx -= splitEvent.size3;
            firstPartials = splitEvent.first5;
            secondPartials = splitEvent.second3;
            pStart = &(splitEvent.pos3[0]);
            pEnd = &(splitEvent.pos3[1]);
        }

        firstPartials[idx].isMajor = 1;
        secondPartials[idx].isMajor = 1;

        if (secondPartials[idx].polyALen > splitEvent.pSpecialData->polyALen)
            splitEvent.pSpecialData->polyALen = secondPartials[idx].polyALen;

        if (hasBoth)
        {
            int32_t spPos = ref.GetSpRefPos(secondPartials[idx].spRefID, secondPartials[idx].refPos);
            int32_t spEnd = ref.GetSpRefPos(secondPartials[idx].spRefID, secondPartials[idx].refEnd);

            if (splitEvent.pSpecialData->pos < 0 || splitEvent.pSpecialData->pos > spPos)
                splitEvent.pSpecialData->pos = spPos;

            if (splitEvent.pSpecialData->end < 0 || splitEvent.pSpecialData->end < spEnd)
                splitEvent.pSpecialData->end = spEnd;
        }

        if (*pStart < 0)
            *pStart = firstPartials[idx].refPos;
        else if (*pStart > firstPartials[idx].refPos)
            *pStart = firstPartials[idx].refPos;

        if (*pEnd < 0)
            *pEnd = firstPartials[idx].refEnd;
        else if (*pEnd < firstPartials[idx].refEnd)
            *pEnd = firstPartials[idx].refEnd;
    }

    if (splitEvent.pos5[0] < 0)
        splitEvent.pos = splitEvent.pos3[0];
    else if (splitEvent.pos3[0] < 0)
        splitEvent.pos = splitEvent.pos5[1];
    else
    {
        splitEvent.pos = splitEvent.pos5[1] < splitEvent.pos3[0] ? splitEvent.pos5[1] : splitEvent.pos3[0];
        if (splitEvent.pos5[1] > splitEvent.pos3[0])
        {
            splitEvent.pSpecialData->tsdPos = splitEvent.pos3[0];
            splitEvent.pSpecialData->tsdLen = splitEvent.pos5[1] - splitEvent.pos3[0] + 1;
        }
    }

    splitEvent.len = 0;
    splitEvent.rpIdx = -1;

    return true;
}
