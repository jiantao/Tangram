/*
 * =====================================================================================
 *
 *       Filename:  TGM_PairAttrbt.cpp
 *
 *    Description:  
 *
 *        Created:  05/11/2012 11:19:10 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *   Inistitution:  Boston College
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

PairAttrbtTable::PairAttrbtTable(const LibTable& libTable, const BamPairTable& bamPairTable)
                                : libTable(libTable), bamPairTable(bamPairTable)
{

}

PairAttrbtTable::~PairAttrbtTable()
{
    delete [] pSpecialAttrbts;
}

void PairAttrbtTable::Init(void)
{
    const Array<char*>* pSpecialRefNames = libTable.GetSpecialRefNames();
    specialSize = pSpecialRefNames->Size() * 2;

    pSpecialAttrbts = new (nothrow) Array<PairAttrbt>[specialSize];
    if (pSpecialAttrbts == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for the pair attribute table.\n");
}

void PairAttrbtTable::MakeInversion(void)
{
    const Array<LocalPair>& invertedPairs = bamPairTable.invertedPairs;

    unsigned int numInv = invertedPairs.Size();
    unsigned int numInv3 = bamPairTable.numInverted3;
    unsigned int numInv5 = numInv - numInv3;

    invertedAttrbts[0].Init(numInv3);
    invertedAttrbts[0].SetSize(numInv3);

    invertedAttrbts[1].Init(numInv5);
    invertedAttrbts[1].SetSize(numInv5);

    unsigned int idx3 = 0;
    unsigned int idx5 = 0;
    for (unsigned int i = 0; i != numInv; ++i)
    {
        const LocalPair* pInvPair = invertedPairs.GetPointer(i);
        unsigned int* pIdx = NULL;
        PairAttrbt* pAttrbt = NULL;

        double fragLen = 0.0;
        double first = 0.0;

        if (pInvPair->readPairType == PT_INVERTED3)
        {
            pIdx = &idx3;
            pAttrbt = invertedAttrbts[0].GetPointer(*pIdx);

            fragLen = pInvPair->end[1] - pInvPair->pos[0] + 1;
            first = pInvPair->pos[0] + fragLen / 2.0;
        }
        else
        {
            pIdx = &idx5;
            pAttrbt = invertedAttrbts[1].GetPointer(*pIdx);

            fragLen = pInvPair->end[0] - pInvPair->pos[1] + 1;
            first = pInvPair->pos[1] + fragLen / 2.0;
        }

        double fragLenMedian = libTable.GetFragLenMedian(pInvPair->readGrpID);

        pAttrbt->origIndex = i;
        pAttrbt->firstAttrbt = first;
        pAttrbt->secondAttrbt = fragLen - fragLenMedian;
        // bound scale is borrowed from Spanner
        pAttrbt->firstBound = fragLenMedian * 1.25;
        pAttrbt->secondBound = fragLenMedian;

        ++(*pIdx);
    }

    invertedAttrbts[0].Sort(CompareAttrbt);
    invertedAttrbts[1].Sort(CompareAttrbt);
}

void PairAttrbtTable::MakeSpecial(void)
{
    const Array<SpecialPair>& specialPairs = bamPairTable.specialPairs;

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

        double halfMedian = (double) libTable.GetFragLenMedian(pSpecialPair->readGrpID) / 2.0;
        double halfMedians[2] = {halfMedian, -halfMedian};
        uint32_t pos[2] = {pSpecialPair->pos[0], pSpecialPair->end[0]};

        int posIndex = pSpecialPair->pairType - PT_SPECIAL3;

        pAttrbt->firstAttrbt = pos[posIndex] + halfMedians[posIndex];
        pAttrbt->secondAttrbt = 0;

        pAttrbt->firstBound = (double) (libTable.GetFragLenHigh(pSpecialPair->readGrpID) - libTable.GetFragLenLow(pSpecialPair->readGrpID)) * boundScale;
        pAttrbt->secondBound = 1e-3;
    }

    // sort the attribute arrays according to the first attribute
    for (unsigned int i = 0; i != specialSize; i += 2)
    {
        pSpecialAttrbts[i].Sort(CompareAttrbt);
        pSpecialAttrbts[i + 1].Sort(CompareAttrbt);
    }
}
