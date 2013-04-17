/*
 * =====================================================================================
 *
 *       Filename:  TGM_FirstMapThread.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  08/01/2012 03:40:36 PM
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
#include "TGM_FirstMapThread.h"

using namespace Tangram;

#ifdef DEBUG

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
            default:
                break;
        }
    }
}

#endif

FirstMapThread::FirstMapThread(const AlignerPars& pars, const LibTable& libInfoTable, const BamPairTable& pairTable, const Reference& reference)
                              : alignerPars(pars), libTable(libInfoTable), bamPairTable(pairTable), ref(reference)
{

}

FirstMapThread::~FirstMapThread()
{

}

void* FirstMapThread::StartThread(void* threadData)
{
    FirstMapData* mapData = (FirstMapData*) threadData;
    FirstMapThread& firstMap = *(mapData->pFirstMapThread);

    bool isUpStream = false;
    bool isOK = true;

#ifdef DEBUG

    std::string cigarStr;

#endif

    RescuePartial rescuePartial(firstMap.alignerPars);

    for (unsigned int i = 0; i != mapData->orphanSize; ++i)
    {
        OrphanPair& orphanPair = mapData->orphanPairs[i];
        RefRegion& refRegion = mapData->refRegions[i];

        switch (orphanPair.aOrient)
        {
            case TGM_1F:
            case TGM_2F:
                isUpStream = false;
                TGM_SeqRevComp(&(orphanPair.read));
                isOK = firstMap.SetFirstRefRegion(refRegion, orphanPair.anchorPos, orphanPair.anchorEnd, orphanPair.readGrpID, orphanPair.read.len, isUpStream);
                break;
            case TGM_1R:
            case TGM_2R:
                isUpStream = true;
                isOK = firstMap.SetFirstRefRegion(refRegion, orphanPair.anchorPos, orphanPair.anchorEnd, orphanPair.readGrpID, orphanPair.read.len, isUpStream);
                break;
            default:
                isOK = false;
                break;
        }

        const AlignerPars& alignerPars = firstMap.alignerPars;
        PrtlAlgnmnt& partial = mapData->firstPartials[i];

        if (isOK)
        {
            s_profile* pProfile = ssw_init(orphanPair.read.seq, orphanPair.read.len, alignerPars.mat, 5, 2);
            s_align* pAlignment = ssw_align(pProfile, refRegion.pRef, refRegion.len, alignerPars.gapOpen, 
                                            alignerPars.gapExt, alignerPars.flag, alignerPars.scoreFilter, alignerPars.distFilter, orphanPair.read.len);

            PartialType partialType;
            bool isRescued = false;

            bool passFilter = firstMap.FirstFilter(partialType, isRescued, rescuePartial, pAlignment, orphanPair, refRegion);
            
            if (passFilter)
            {
                if (!isRescued)
                {
                    partial.refPos = pAlignment->ref_begin1 + refRegion.start;
                    partial.refEnd = pAlignment->ref_end1 + refRegion.start;
                    partial.readPos = pAlignment->read_begin1;
                    partial.readEnd = pAlignment->read_end1;
                    partial.cigar = pAlignment->cigar;
                    partial.cigarLen = pAlignment->cigarLen;

#ifdef DEBUG
                CigarToString(cigarStr, pAlignment->cigar, pAlignment->cigarLen);
                printf("chr%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n", orphanPair.refID + 1, orphanPair.anchorPos, orphanPair.anchorEnd + 1, refRegion.start, 
                        refRegion.start + refRegion.len, pAlignment->ref_begin1 + refRegion.start, pAlignment->ref_end1 + refRegion.start, pAlignment->read_begin1, 
                        pAlignment->read_end1, pAlignment->score1, cigarStr.c_str());
#endif

                    free(pAlignment);
                }
                else
                {

#ifdef DEBUG
                    CigarToString(cigarStr, pAlignment->cigar, pAlignment->cigarLen);
                    printf("chr%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t", orphanPair.refID + 1, orphanPair.anchorPos, orphanPair.anchorEnd + 1, refRegion.start, 
                            refRegion.start + refRegion.len, pAlignment->ref_begin1 + refRegion.start, pAlignment->ref_end1 + refRegion.start, pAlignment->read_begin1, 
                            pAlignment->read_end1, pAlignment->score1, cigarStr.c_str());

                    CigarToString(cigarStr, rescuePartial.cigar, rescuePartial.cigarLen);
                    printf("%s\n", cigarStr.c_str());
#endif

                    align_destroy(pAlignment);

                    partial.refPos = rescuePartial.refPos + refRegion.start;
                    partial.refEnd = rescuePartial.refEnd + refRegion.start;
                    partial.readPos = rescuePartial.readPos;
                    partial.readEnd = rescuePartial.readEnd;
                    partial.cigar = rescuePartial.cigar;
                    partial.cigarLen = rescuePartial.cigarLen;

                    isRescued = false;
                    rescuePartial.Clear();
                }

                partial.origIdx = i + mapData->orphanStart;
                partial.isReversed = isUpStream ? 0 : 1;
                partial.isSoft = 0;
                partial.partialType = partialType;

            }
            else
            {

#ifdef DEBUG
                CigarToString(cigarStr, pAlignment->cigar, pAlignment->cigarLen);
                printf("chr%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n", orphanPair.refID + 1, orphanPair.anchorPos, orphanPair.anchorEnd + 1, refRegion.start, 
                        refRegion.start + refRegion.len, pAlignment->ref_begin1 + refRegion.start, pAlignment->ref_end1 + refRegion.start, pAlignment->read_begin1, 
                        pAlignment->read_end1, pAlignment->score1, cigarStr.c_str());
#endif

                partial.cigar = NULL;
                partial.refPos = INT32_MAX;
                align_destroy(pAlignment);
            }

            init_destroy(pProfile);
        }
        else
        {
            partial.cigar = NULL;
            partial.refPos = INT32_MAX;
        }
    }

    pthread_exit(NULL);
}

bool FirstMapThread::SetFirstRefRegion(RefRegion& refRegion, int32_t anchorPos, int32_t anchorEnd, int32_t readGrpID, uint32_t readLen, bool isUpStream) const
{
    uint32_t fragLenHigh = libTable.GetFragLenHigh(readGrpID);
    // uint32_t fragLenLow = libTable.GetFragLenLow(readGrpID);
    int32_t refEnd = ref.pos + ref.refSeq.Size() - 1;

    if (!isUpStream)
    {
        int32_t start = anchorEnd + 1;
        int32_t end = anchorPos + fragLenHigh + 2 * readLen;

        if (end >= refEnd)
            end = refEnd;

        int32_t regionLen = end - start + 1;

        if (regionLen <= 0 || regionLen < (int) readLen)
            return false;

        refRegion.pRef = ref.refSeq.GetPointer(start - ref.pos);
        refRegion.len = regionLen;
        refRegion.start = start;
    }
    else
    {
        int32_t start = anchorEnd - fragLenHigh - 2 * readLen;
        int32_t end = anchorPos - 1;

        if (start < 0)
            start = 0;

        if (end < 0)
            end = 0;

        int32_t regionLen = end - start + 1;
        if (regionLen <= 0 || regionLen < (int) readLen)
            return false;

        refRegion.pRef = ref.refSeq.GetPointer(start - ref.pos);
        refRegion.len = regionLen;
        refRegion.start = start;
    }

    return true;
}

bool FirstMapThread::FirstFilter(PartialType& partialType, bool& isRescued, RescuePartial& rescuePartial,
                                 const s_align* pAlignment, const OrphanPair& orphanPair, const RefRegion& refRegion) const
{
    isRescued = false;

    // best alignment is not available
    if (pAlignment->ref_begin1 < 0 || pAlignment->read_begin1 < 0)
    {
        return false;
    }

    // length filter
    int alignedReadLen = pAlignment->read_end1 - pAlignment->read_begin1 + 1;
    int readLen = orphanPair.read.len;

    if (alignedReadLen < alignerPars.minAlignedLen)
        return false;

    double scoreRate = (double) pAlignment->score1 / (alignerPars.mat[0] * alignedReadLen);
    partialType = GetPartialType(pAlignment->read_begin1, pAlignment->read_end1, orphanPair.read.len);

    if (scoreRate >= alignerPars.minScoreRate && pAlignment->cigarLen < 3)
    {
        if (readLen - alignedReadLen < alignerPars.minAlignedLen)
                return false;

        if (partialType == PARTIAL_UNKNOWN)
            return false;
    }
    else
    {
        isRescued = rescuePartial.RescueLowScore(partialType, pAlignment, orphanPair.read.seq, orphanPair.read.len, refRegion);
        if (!isRescued)
            return false;
    }

    return true;
}
