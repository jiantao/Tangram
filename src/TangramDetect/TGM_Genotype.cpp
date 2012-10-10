/*
 * =====================================================================================
 *
 *       Filename:  TGM_Genotype.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  09/21/2012 11:55:54 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include "TGM_Genotype.h"

using namespace BamTools;
using namespace Tangram;

Genotype::Genotype(BamMultiReader& reader, const GenotypePars& genotypePars, const LibTable& libTable, BamPairTable& bamPairTable)
          : reader(reader), genotypePars(genotypePars), libTable(libTable), bamPairTable(bamPairTable)
{

    lastChr = -1;
    lastPos = -1;
    lastEnd = -1;

    specialPrior[0] = 1.0/3.0;
    specialPrior[1] = 1.0/3.0;
    specialPrior[2] = 1.0/3.0;
}

Genotype::~Genotype()
{

}

void Genotype::Init(void)
{
    unsigned int numSamples = libTable.GetNumSamples();

    sampleCount.Init(numSamples);
    sampleCount.SetSize(numSamples);
}

void Genotype::SetSpecialPrior(const double* prior)
{
    specialPrior[0] = prior[0];
    specialPrior[1] = prior[1];
    specialPrior[2] = prior[2];
}

bool Genotype::Special(const SpecialEvent* pRpSpecial, const SplitEvent* pSplitEvent)
{
    // clear the fragment count
    sampleCount.MemSet(0);

    int32_t chr = -1;
    int32_t pos = -1;

    bool isPresice = false;
    if (pRpSpecial != NULL)
    {
        chr = pRpSpecial->refID;
        pos = pRpSpecial->pos;
        SetSampleCountSpecial(*pRpSpecial);
    }

    if (pSplitEvent != NULL)
    {
        isPresice = true;
        chr = pSplitEvent->refID;
        pos = pSplitEvent->pos;
        SetSampleCountSplit(*pSplitEvent);
    }

    if (genotypePars.doGenotype)
    {
        if (!SpecialFilter(pRpSpecial, pSplitEvent))
            return false;

        int32_t fragLenMax = libTable.GetFragLenMax();

        int32_t jumpPos = pos - fragLenMax;
        if (jumpPos < 0)
            jumpPos = 0;

        if (!Jump(chr, jumpPos))
            return false;

        int32_t posUpper = pos;
        int32_t posLower = pos;

        if (!isPresice)
        {
            posUpper -= pRpSpecial->posUncertainty / 2;
            posLower += pRpSpecial->posUncertainty / 2;
        }

        BamAlignment alignment;
        while (reader.GetNextAlignment(alignment))
        {
            if (alignment.Position > pos + 100)
                break;

            int32_t fragEnd = 0;
            int32_t alignEnd = alignment.GetEndPosition(false, true);

            bool isUpperMate = false;
            if (alignment.RefID != alignment.MateRefID)
                continue;
            else if (alignment.Position < alignment.MatePosition)
            {
                isUpperMate = true;
                fragEnd = alignment.Position + alignment.InsertSize - 1;
            }
            else
            {
                fragEnd = alignEnd;
                if (!isPresice)
                    continue;
            }

            if (fragEnd < posUpper)
                continue;

            int32_t readGrpID = -1;
            PairType pairType = bamPairTable.CheckPairType(readGrpID, alignment);

            if (isPresice)
            {
                if (alignment.Position < pos && alignEnd > pos && alignment.MapQuality >= genotypePars.minMQ && pairType != PT_SOFT3 && pairType != PT_SOFT5)
                {
                    int upLen = pos - alignment.Position + 1;
                    int downLen = alignEnd - pos + 1;
                    if (upLen > genotypePars.minCrossLen && downLen > genotypePars.minCrossLen)
                        UpdateNonSupport(readGrpID);
                }
            }

            if (isUpperMate)
            {
                if (pairType == PT_NORMAL && alignEnd <= posUpper && alignment.MatePosition >= posLower)
                    UpdateNonSupport(readGrpID);
                else if (pairType == PT_SHORT && alignEnd <= posUpper && alignment.MatePosition >= posLower)
                    UpdateSupport(readGrpID);
            }
        }

        SetLikelihood();

        // update the bam file stream position
        lastChr = chr;
        lastPos = jumpPos;
        lastEnd = alignment.Position;
    }

    return true;
}

bool Genotype::SpecialFilter(const SpecialEvent* pRpSpecial, const SplitEvent* pSplitEvent) const
{
    bool rpPass = false;
    if (pRpSpecial != NULL)
    {
        if (pRpSpecial->numFrag[0] >= genotypePars.minRpFrag && pRpSpecial->numFrag[1] >= genotypePars.minRpFrag)
            rpPass = true;
    }

    bool srPass = false;
    if (pSplitEvent != NULL)
    {
        if (pSplitEvent->size3 + pSplitEvent->size5 >= genotypePars.minSrFrag)
            srPass = true;
    }

    if (!srPass && !rpPass)
        return false;

    return true;
}

bool Genotype::Jump(int32_t refID, int32_t pos)
{
    if (refID != lastChr || pos >= lastEnd + genotypePars.minJumpLen || pos < lastEnd)
    {
        if (reader.Jump(refID, pos))
            return true;
        else
            return false;
    }
    else
    {
        BamAlignment alignment;
        while (reader.GetNextAlignment(alignment))
        {
            if (alignment.Position >= pos)
                break;
        }
    }

    return true;
}

void Genotype::SetSampleCountSpecial(const SpecialEvent& rpSpecial)
{
    for (unsigned int k = 0; k != rpSpecial.numFrag[0]; ++k)
    {
        const SpecialPair& specialPair = bamPairTable.specialPairs[rpSpecial.origIndex[0][k]];
        unsigned int readGrpID = specialPair.readGrpID;
        unsigned int sampleID = 0;

        if (libTable.GetSampleID(sampleID, readGrpID))
            ++(sampleCount[sampleID].support);
    }

    for (unsigned int k = 0; k != rpSpecial.numFrag[1]; ++k)
    {
        const SpecialPair& specialPair = bamPairTable.specialPairs[rpSpecial.origIndex[1][k]];
        unsigned int readGrpID = specialPair.readGrpID;
        unsigned int sampleID = 0;

        if (libTable.GetSampleID(sampleID, readGrpID))
            ++(sampleCount[sampleID].support);
    }
}

void Genotype::SetSampleCountSplit(const SplitEvent& splitEvent)
{
    for (unsigned int i = 0; i != splitEvent.size3; ++i)
    {
        unsigned int readGrpID = 0;
        unsigned int sampleID = 0;
        unsigned int origIdx = splitEvent.first3[i].origIdx;

        if (splitEvent.first3[i].isSoft)
            readGrpID = bamPairTable.softPairs[origIdx].readGrpID;
        else
            readGrpID = bamPairTable.orphanPairs[origIdx].readGrpID;

        if (libTable.GetSampleID(sampleID, readGrpID))
            ++(sampleCount[sampleID].support);
    }

    for (unsigned int i = 0; i != splitEvent.size5; ++i)
    {
        unsigned int readGrpID = 0;
        unsigned int sampleID = 0;
        unsigned int origIdx = splitEvent.first5[i].origIdx;

        if (splitEvent.first5[i].isSoft)
            readGrpID = bamPairTable.softPairs[origIdx].readGrpID;
        else
            readGrpID = bamPairTable.orphanPairs[origIdx].readGrpID;

        if (libTable.GetSampleID(sampleID, readGrpID))
            ++(sampleCount[sampleID].support);
    }
}

void Genotype::SetLikelihood(void)
{
    genotypes.clear();
    likelihoods.clear();

    unsigned int size = sampleCount.Size();
    for (unsigned int i = 0; i != size; ++i)
    {
        unsigned int refCount = sampleCount[i].nonSupport;
        unsigned int altCount = sampleCount[i].support;

        if (refCount == 0 && altCount == 0)
        {
            // no evidence at all
            genotypes.push_back(-1);
            likelihoods.push_back(-0.48);
            likelihoods.push_back(-0.48);
            likelihoods.push_back(-0.48);
        }
        else
        {
            int genotype = -1;
            long double maxLikelihood = 0;
            for (unsigned j = 0; j != 3; ++j)
            {
                long double likelihood = CalculateLikelihood(refCount, altCount, genotypePars.p[j]);
                likelihoods.push_back(log10(likelihood));

                if (likelihood > maxLikelihood)
                {
                    maxLikelihood = likelihood;
                    genotype = j;
                }
            }

            genotypes.push_back(genotype);
        }
    }
}

long double Genotype::CalculateLikelihood(unsigned int refCount, unsigned int altCount, double p) const
{
    unsigned int n = refCount + altCount;
    long double cnm = 1.0;

    unsigned int m = refCount;
    if (m > altCount)
        m = altCount;

    for (unsigned int i = 1, j = n; i <= m; --j, ++i)
    {
        cnm *= j;
        cnm /= i;
    }

    return (cnm * pow((long double) p, refCount) * pow((long double) (1.0 - p), altCount));
}
