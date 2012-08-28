/*
 * =====================================================================================
 *
 *       Filename:  TGM_PairAttrbt.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  05/11/2012 11:19:12 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef  TGM_PAIRATTRBT_H
#define  TGM_PAIRATTRBT_H

#include <vector>
#include "TGM_BamPair.h"

namespace Tangram
{
    struct PairAttrbt
    {
        double firstAttrbt;

        double secondAttrbt;

        double firstBound;

        double secondBound;

        uint64_t origIndex;
    };

    class PairAttrbtTable
    {
        public:

            PairAttrbtTable();

            ~PairAttrbtTable();

            void Init(const LibTable* pLibTable);

            void MakeSpecial(const Array<SpecialPair>& specialPairs);

        public:

            Array<PairAttrbt>* pSpecialAttrbts;

            unsigned int specialSize;

            Array<PairAttrbt> longAttrbts;

            Array<PairAttrbt> shortAttrbts;

            Array<PairAttrbt> invertedAttrbts[2];

        private:

            const LibTable* pLibTable;
    };
};

#endif  /*TGM_PAIRATTRBT_H*/
