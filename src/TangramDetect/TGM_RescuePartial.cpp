/*
 * =====================================================================================
 *
 *       Filename:  TGM_RescuePartial.cpp
 *
 *    Description:  Rescue the partial alignments with low sw score
 *
 *        Created:  08/01/2012 11:02:20 PM
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Jiantao Wu (), 
 *   Inistitution:  Boston College
 *
 * =====================================================================================
 */

#include "TGM_BamPair.h"
#include "TGM_RescuePartial.h"

using namespace Tangram;


RescuePartial::RescuePartial(const AlignerPars& pars) : alignerPars(pars)
{

}
           
RescuePartial::~RescuePartial()
{

}

bool RescuePartial::RescueLowScore(PartialType& partialType, const s_align* pAlignment, const int8_t* readSeq, int readLen, const RefRegion& refRegion)
{
    Init(readLen);

    DistributeScore(pAlignment, readSeq, readLen, refRegion);

    int start = 0;
    int end = 0;

    if (FindMaxSubArray(start, end))
    {
        SetRescueData(pAlignment, start, end);
        if (!RescueFilter(partialType, readLen))
            return false;
    }
    else
        return false;

    return true;
}

void RescuePartial::Init(unsigned int readLen)
{
    score.Clear();
    len.Clear();
    type.Clear();

    cigar = NULL;
    cigarLen = 0;

    if (score.Capacity() < readLen)
    {
        score.ResizeNoCopy(readLen);
        len.ResizeNoCopy(readLen);
        type.ResizeNoCopy(readLen);
    }
}

void RescuePartial::DistributeScore(const s_align* pAlignment, const int8_t* readSeq, int readLen, const RefRegion& refRegion)
{
    const int matchScore = alignerPars.mat[0];

    // heavy mismatch and gap penalty for a cleaner alignment edge
    const int mismatchScore = -2 * alignerPars.mat[0];
    const int gapOpenScore = -3 * alignerPars.mat[0];

    int scoreIdx = 0;
    int refPos = pAlignment->ref_begin1;
    int readPos = pAlignment->read_begin1;
    for (int i = 0; i != pAlignment->cigarLen; ++i)
    {
        uint8_t cigarType = pAlignment->cigar[i] & BAM_CIGAR_MASK ;
        int cigarTypeLen = (pAlignment->cigar[i] >> BAM_CIGAR_SHIFT);
        int8_t lastScore = 0;
        int start = 0;
        switch(cigarType)
        {
            case BAM_CINS:
                readPos += cigarTypeLen;
                type[scoreIdx] = cigarType;
                len[scoreIdx] = cigarTypeLen;
                score[scoreIdx] = gapOpenScore + mismatchScore * (cigarTypeLen - 1);

                ++scoreIdx;
                type.Increment();
                len.Increment();
                score.Increment();
                break;
            case BAM_CDEL:
                refPos += cigarTypeLen;
                type[scoreIdx] = cigarType;
                len[scoreIdx] = cigarTypeLen;
                score[scoreIdx] = gapOpenScore + mismatchScore * (cigarTypeLen - 1);

                ++scoreIdx;
                type.Increment();
                len.Increment();
                score.Increment();
                break;
            case BAM_CMATCH:
                for (int j = 0; j != cigarTypeLen; ++j)
                {
                    int refIdx = refPos + j;
                    int readIdx = readPos + j;
                    int currScore = 0;
                    if (refRegion.pRef[refIdx] == readSeq[readIdx] && refRegion.pRef[refIdx] != 4)
                        currScore = matchScore;
                    else
                        currScore = mismatchScore;

                    if (j != 0)
                    {
                        if (lastScore != currScore)
                        {
                            type[scoreIdx] = lastScore > 0 ? BAM_CSEQ_MATCH : BAM_CMISMATCH;
                            len[scoreIdx] = j - start;
                            score[scoreIdx] = len[scoreIdx] * lastScore;

                            start = j;
                            ++scoreIdx;
                            type.Increment();
                            len.Increment();
                            score.Increment();
                        }
                    }

                    lastScore = currScore;
                }

                type[scoreIdx] = lastScore > 0 ? BAM_CSEQ_MATCH : BAM_CMISMATCH;
                len[scoreIdx] = cigarTypeLen - start;
                score[scoreIdx] = len[scoreIdx] * lastScore;

                ++scoreIdx;
                type.Increment();
                len.Increment();
                score.Increment();

                refPos += cigarTypeLen;
                readPos += cigarTypeLen;
                break;
            default:
                break;
        }
    }
}

bool RescuePartial::FindMaxSubArray(int& start, int& end)
{
    int maxSum = 0;
    int currSum = 0;
    int currStart = 0;

    start = 0;
    end = 0;

    unsigned int size = score.Size();
    for (unsigned int i = 0; i != size; ++i)
    {
        currSum += score[i];

        if (currSum > maxSum)
        {
            maxSum = currSum;
            start = currStart;
            end = i;
        }
        else if (currSum <= 0)
        {
            currStart = i + 1;
            currSum = 0;
        }
    }

    if (maxSum > 0)
        return true;
    else
        return false;
}

void RescuePartial::SetRescueData(const s_align* pAlignment, int start, int end)
{
    bestScore = 0;
    refPos = pAlignment->ref_begin1;
    refEnd = pAlignment->ref_end1;
    readPos = pAlignment->read_begin1;
    readEnd = pAlignment->read_end1;

    for (int i = 0; i != start; ++i)
    {
        switch(type[i])
        {
            case BAM_CINS:
                readPos += len[i];
                break;
            case BAM_CDEL:
                refPos += len[i];
                break;
            case BAM_CSEQ_MATCH:
            case BAM_CMISMATCH:
                readPos += len[i];
                refPos += len[i];
                break;
            default:
                break;
        }
    }

    for (int i = len.Size() - 1; i != end; --i)
    {
        switch(type[i])
        {
            case BAM_CINS:
                readEnd -= len[i];
                break;
            case BAM_CDEL:
                refEnd -= len[i];
                break;
            case BAM_CSEQ_MATCH:
            case BAM_CMISMATCH:
                readEnd -= len[i];
                refEnd -= len[i];
                break;
            default:
                break;
        }
    }

    bool isLastMatch = false;
    int moveIdx = 0;
    for (int i = start; i <= end; ++i)
    {
        switch(type[i])
        {
            case BAM_CINS:
            case BAM_CDEL:
                if (isLastMatch)
                {
                    ++moveIdx;
                    isLastMatch = false;
                }

                bestScore += -alignerPars.gapOpen - alignerPars.gapExt * (len[i] - 1);
                len[moveIdx] = len[i];
                type[moveIdx] = type[i];
                ++moveIdx;
                break;
            case BAM_CSEQ_MATCH:
                bestScore += alignerPars.mat[0] * len[i];
                if (isLastMatch)
                    len[moveIdx] += len[i];
                else
                {
                    isLastMatch = true;
                    len[moveIdx] = len[i];
                    type[moveIdx] = BAM_CMATCH;
                }
                break;
            case BAM_CMISMATCH:
                bestScore += -alignerPars.mat[0] * len[i];
                if (isLastMatch)
                    len[moveIdx] += len[i];
                else
                {
                    isLastMatch = true;
                    len[moveIdx] = len[i];
                    type[moveIdx] = BAM_CMATCH;
                }
                break;
            default:
                break;
        }
    }

    if (isLastMatch)
    {
        len.SetSize(moveIdx + 1);
        type.SetSize(moveIdx + 1);
    }
    else
    {
        len.SetSize(moveIdx);
        type.SetSize(moveIdx);
    }
}

bool RescuePartial::RescueFilter(PartialType& partialType, int readLen)
{
    // length filter
    int alignedLen = readEnd - readPos + 1;
    if (alignedLen < alignerPars.minAlignedLen || (int) readLen - alignedLen < alignerPars.minAlignedLen)
        return false;

    // score filter
    double scoreRate = (double) bestScore / (alignedLen * alignerPars.mat[0]);
    if (scoreRate < alignerPars.minScoreRate)
        return false;

    // type filter
    partialType = GetPartialType(readPos, readEnd, readLen);
    if (partialType == PARTIAL_UNKNOWN)
        return false;

    // create the cigar info
    cigarLen = len.Size();
    cigar = (uint32_t*) malloc(cigarLen * sizeof(uint32_t));
    if (cigar == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for the cigar data.\n");

    for (unsigned int i = 0; i != cigarLen; ++i)
    {
        cigar[i] = len[i];
        cigar[i] = cigar[i] << BAM_CIGAR_SHIFT;
        cigar[i] |= type[i];
    }

    return true;
}
