/*
 * =====================================================================================
 *
 *       Filename:  TGM_PairAttrbt.h
 *
 *    Description:  
 *
 *        Created:  05/11/2012 11:19:12 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *   Inistitution:  Boston College
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

            PairAttrbtTable(const LibTable& libTable, const BamPairTable& bamPairTable);

            ~PairAttrbtTable();

            void Init(void);

            void MakeInversion(void);

            void MakeSpecial(void);

        public:

            Array<PairAttrbt>* pSpecialAttrbts;

            unsigned int specialSize;

            Array<PairAttrbt> longAttrbts;

            Array<PairAttrbt> shortAttrbts;

            Array<PairAttrbt> invertedAttrbts[2];

        private:

            const LibTable& libTable;

            const BamPairTable& bamPairTable;
    };
};

#endif  /*TGM_PAIRATTRBT_H*/
