/*
 * =====================================================================================
 *
 *       Filename:  TGM_Cluster.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  05/11/2012 10:31:50 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include "khash.h"
#include <cmath>
#include "TGM_Cluster.h"

using namespace Tangram;

#define DEFAULT_CLS_ELEMENTS_CAP 20

KHASH_MAP_INIT_INT(clusterMap, int);

Cluster::Cluster()
{

}

Cluster::~Cluster()
{

}

void Cluster::Init(const Array<PairAttrbt>* pPairAttrbts, unsigned int minClusterSize, const double* minStd)
{
    this->minClusterSize = minClusterSize;

    this->minStd[0] = minStd[0];
    this->minStd[1] = minStd[1];

    elements.Init(DEFAULT_CLS_ELEMENTS_CAP);

    this->pPairAttrbts = pPairAttrbts;

    unsigned int newSize = pPairAttrbts->Size();
    if (next.Size() < newSize)
    {
        next.ResizeNoCopy(newSize * 2);
        map.ResizeNoCopy(newSize * 2);
        count.ResizeNoCopy(newSize * 2);
    }
    else
    {
        next.Clear();
        map.Clear();
        count.Clear();
    }

    next.SetSize(newSize);
    map.SetSize(newSize);
    count.SetSize(newSize);

    for (unsigned int i = 0; i != newSize; ++i)
    {
        next[i] = i;
        map[i] = i;
        count[i] = 0;
    }
}

void Cluster::Make(void)
{
    Build();
    Finalize();
    Clean();
}

void Cluster::Build(void)
{
    unsigned int lastLowIndex = 0;
    const unsigned int attrbtSize = pPairAttrbts->Size();
    for (unsigned int i = 0; i != attrbtSize; ++i)
    {
        unsigned int j = lastLowIndex;

        double firstBound = (*pPairAttrbts)[i].firstBound;
        double secondBound = (*pPairAttrbts)[i].secondBound;

        while ((j < attrbtSize) && ((*pPairAttrbts)[i].firstAttrbt - (*pPairAttrbts)[j].firstAttrbt) > firstBound)
            ++j;

        lastLowIndex = j;

        while (fabs((*pPairAttrbts)[i].firstAttrbt - (*pPairAttrbts)[j].firstAttrbt) <= firstBound)
        {
            if (fabs((*pPairAttrbts)[i].secondAttrbt - (*pPairAttrbts)[j].secondAttrbt) <= secondBound)
                ++(count[i]);

            ++j;

            if (j == attrbtSize)
                break;
        }
    }

    lastLowIndex = 0;
    for (unsigned int i = 0; i != attrbtSize; ++i)
    {
        unsigned int j = lastLowIndex;

        double firstBound = (*pPairAttrbts)[i].firstBound;
        double secondBound = (*pPairAttrbts)[i].secondBound;

        while ((j < attrbtSize) && ((*pPairAttrbts)[i].firstAttrbt - (*pPairAttrbts)[j].firstAttrbt) > firstBound)
            ++j;

        lastLowIndex = j;

        unsigned int maxCount = 0;
        unsigned int centerIndex = 0;

        while (fabs((*pPairAttrbts)[i].firstAttrbt - (*pPairAttrbts)[j].firstAttrbt) <= firstBound)
        {
            if (fabs((*pPairAttrbts)[i].secondAttrbt - (*pPairAttrbts)[j].secondAttrbt) <= secondBound)
            {
                if (count[j] > maxCount)
                {
                    centerIndex = j;
                    maxCount = count[j];
                }
            }

            ++j;

            if (j == attrbtSize)
                break;
        }

        Connect(centerIndex, i);
    }
}

void Cluster::Connect(unsigned int centerIndex, unsigned int memberIndex)
{
    if (centerIndex == memberIndex)
        return;

    int centerID = map[centerIndex];
    int centerNext = next[centerIndex];

    next[centerIndex] = memberIndex;

    unsigned int i = memberIndex;
    for (; next[i] != memberIndex; i = next[i])
    {
        map[i] = centerID;
    }

    next[i] = centerNext;
    map[i] = centerID;
}

void Cluster::Finalize(void)
{
    // create a hash to transfer the cluster center ID to cluster element ID
    khash_t(clusterMap)* clusterHash = kh_init(clusterMap);
    kh_resize(clusterMap, clusterHash, 50);

    int khRet = 0;
    khiter_t khIter = 0;
    int clusterIndex = 0;
    unsigned int attrbtSize = pPairAttrbts->Size();

    for (unsigned int i = 0; i != attrbtSize; ++i)
    {
        // skip those very small clusters 
        unsigned int clusterSize = count[map[i]];
        if (clusterSize < minClusterSize)
            continue;

        khIter = kh_put(clusterMap, clusterHash, map[i], &khRet);
        if (khRet == 0)
        {
            clusterIndex = kh_value(clusterHash, khIter);
            UpdateElmnt(clusterIndex, i);
        }
        else
        {
            clusterIndex = AddElmnt(i);
            kh_value(clusterHash, khIter) = clusterIndex;
        }

        map[i] = clusterIndex;
    }

    kh_destroy(clusterMap, clusterHash);

    unsigned int numElmnts = elements.Size();
    for (unsigned int i = 0; i != numElmnts; ++i)
    {
        ClusterElmnt* pElmnt = elements.GetPointer(i);

        for (unsigned int j = 0; j != 2; ++j)
        {
            pElmnt->mean[j] /= pElmnt->numReadPair;
            pElmnt->std[j] = sqrt(pElmnt->std[j] / pElmnt->numReadPair - (pElmnt->mean[j] * pElmnt->mean[j]));
        }
    }
}

unsigned int Cluster::AddElmnt(unsigned int attrbtIdx)
{
    int clusterID = elements.Size();

    if (elements.IsFull())
        elements.Resize(clusterID * 2);

    elements.Increment();

    ClusterElmnt* pElmnt = elements.GetPointer(clusterID);
    const PairAttrbt* pAttrbt = pPairAttrbts->GetPointer(attrbtIdx);

    pElmnt->numReadPair = 1;
    pElmnt->startIndex = attrbtIdx;

    pElmnt->mean[0] = pAttrbt->firstAttrbt;
    pElmnt->mean[1] = pAttrbt->secondAttrbt;

    pElmnt->std[0] = pAttrbt->firstAttrbt * pAttrbt->firstAttrbt;
    pElmnt->std[1] = pAttrbt->secondAttrbt * pAttrbt->secondAttrbt;

    pElmnt->min[0] = pAttrbt->firstAttrbt;
    pElmnt->min[1] = pAttrbt->secondAttrbt;

    pElmnt->max[0] = pAttrbt->firstAttrbt;
    pElmnt->max[1] = pAttrbt->secondAttrbt;

    return clusterID;
}

void Cluster::UpdateElmnt(unsigned int clusterID, unsigned int attrbtIdx)
{
    ClusterElmnt* pElmnt = elements.GetPointer(clusterID);
    const PairAttrbt* pAttrbt = pPairAttrbts->GetPointer(attrbtIdx);

    ++(pElmnt->numReadPair);

    pElmnt->mean[0] += pAttrbt->firstAttrbt;
    pElmnt->mean[1] += pAttrbt->secondAttrbt;

    pElmnt->std[0] += pAttrbt->firstAttrbt * pAttrbt->firstAttrbt;
    pElmnt->std[1] += pAttrbt->secondAttrbt * pAttrbt->secondAttrbt;

    pElmnt->min[0] = (pAttrbt->firstAttrbt < pElmnt->min[0] ? pAttrbt->firstAttrbt : pElmnt->min[0]);
    pElmnt->min[1] = (pAttrbt->secondAttrbt < pElmnt->min[1] ? pAttrbt->secondAttrbt : pElmnt->min[1]);

    pElmnt->max[0] = (pAttrbt->firstAttrbt > pElmnt->max[0] ? pAttrbt->firstAttrbt : pElmnt->max[0]);
    pElmnt->max[1] = (pAttrbt->secondAttrbt > pElmnt->max[1] ? pAttrbt->secondAttrbt : pElmnt->max[1]);
}

void Cluster::Clean(void)
{
    actualNumElmnts = elements.Size();
    unsigned int numElmnts = elements.Size();

    // lazy deletion
    // the deleted elements' number of read pairs will be set to zero
    for (unsigned int i = 0; i != numElmnts; ++i)
    {
        if (elements[i].numReadPair < minClusterSize
            || elements[i].std[0] <= minStd[0]
            || elements[i].std[1] <= minStd[1]
            || (elements[i].std[0] != elements[i].std[0] && elements[i].std[1] != elements[i].std[1]))
        {
            elements[i].numReadPair = 0;
            --actualNumElmnts;
        }
    }
}

void Cluster::Print(int32_t refID, const char* specialID) const
{
    for (unsigned int i = 0; i != elements.Size(); ++i)
    {
        int startIdx = elements[i].startIndex;
        int firstAttrbt = (*pPairAttrbts)[startIdx].firstAttrbt;
        printf("%d\t%d\t%d\t%d\t%f\n", refID + 1, firstAttrbt, firstAttrbt + 1, elements[i].numReadPair, elements[i].std[0]);
    }

    printf("\n");
}
