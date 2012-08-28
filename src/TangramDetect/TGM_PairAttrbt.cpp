/*
 * =====================================================================================
 *
 *       Filename:  TGM_PairAttrbt.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  05/11/2012 11:19:10 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include "TGM_Error.h"
#include "TGM_PairAttrbt.h"

using namespace std;
using namespace Tangram;

const static unsigned int DEFALUT_NUM_RP = 50;

static int CompareAttrbt(const void* a, const void* b)
{
    const PairAttrbt* first = (const PairAttrbt*) a;
    const PairAttrbt* second = (const PairAttrbt*) b;

    if (first->firstAttrbt < second->firstAttrbt)
    {
        return -1;
    }
    else if (first->firstAttrbt > second->firstAttrbt)
    {
        return 1;
    }
    else
    {
        if (first->secondAttrbt < second->secondAttrbt)
            return -1;
        else if (first->secondAttrbt > second->secondAttrbt)
            return 1;
        else
        {
            return 0;
        }
    }
}

PairAttrbtTable::PairAttrbtTable()
{
    pSpecialAttrbts = NULL;
}

PairAttrbtTable::~PairAttrbtTable()
{
    delete [] pSpecialAttrbts;
}

void PairAttrbtTable::Init(const LibTable* pLibTable)
{
    this->pLibTable = pLibTable;

    const Array<char*>* pSpecialRefNames = pLibTable->GetSpecialRefNames();
    specialSize = pSpecialRefNames->Size() * 2;

    pSpecialAttrbts = new (nothrow) Array<PairAttrbt>[specialSize];
    if (pSpecialAttrbts == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for the pair attribute table.\n");
}

void PairAttrbtTable::MakeSpecial(const Array<SpecialPair>& specialPairs)
{
    for (unsigned int i = 0; i != specialSize; i += 2)
    {
        pSpecialAttrbts[i].Init(DEFALUT_NUM_RP);
        pSpecialAttrbts[i + 1].Init(DEFALUT_NUM_RP);
    }

    double boundScale = 1.25;
    unsigned int numSpecialPairs = specialPairs.Size();

    for (unsigned int i = 0; i != numSpecialPairs; ++i)
    {
        const SpecialPair* pSpecialPair = specialPairs.GetPointer(i);

        int arrayIndex = pSpecialPair->specialID * 2 + pSpecialPair->pairType - PT_SPECIAL3;

        unsigned int attrbtSize = pSpecialAttrbts[arrayIndex].Size();
        if (pSpecialAttrbts[arrayIndex].IsFull())
            pSpecialAttrbts[arrayIndex].Resize(attrbtSize * 2);

        PairAttrbt* pAttrbt = pSpecialAttrbts[arrayIndex].GetPointer(attrbtSize);
        pSpecialAttrbts[arrayIndex].Increment();

        pAttrbt->origIndex = i;

        double halfMedian = (double) pLibTable->GetFragLenMedian(pSpecialPair->readGrpID) / 2.0;
        double halfMedians[2] = {halfMedian, -halfMedian};
        uint32_t pos[2] = {pSpecialPair->pos[0], pSpecialPair->end[0]};

        int posIndex = pSpecialPair->pairType - PT_SPECIAL3;

        pAttrbt->firstAttrbt = pos[posIndex] + halfMedians[posIndex];
        pAttrbt->secondAttrbt = 0;

        pAttrbt->firstBound = (double) (pLibTable->GetFragLenHigh(pSpecialPair->readGrpID) - pLibTable->GetFragLenLow(pSpecialPair->readGrpID)) * boundScale;
        pAttrbt->secondBound = 1e-3;
    }

    // sort the attribute arrays according to the first attribute
    for (unsigned int i = 0; i != specialSize; i += 2)
    {
        pSpecialAttrbts[i].Sort(CompareAttrbt);
        pSpecialAttrbts[i + 1].Sort(CompareAttrbt);
    }
}
