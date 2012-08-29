/*
 * =====================================================================================
 *
 *       Filename:  TGM_Printer.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  08/16/2012 03:35:17 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef  TGM_PRINTER_H
#define  TGM_PRINTER_H

#include "TGM_LibTable.h"
#include "TGM_BamPair.h"
#include "TGM_Detector.h"
#include "TGM_Aligner.h"
#include "TGM_Reference.h"

namespace Tangram
{
    struct PrintIdx
    {
        union
        {
            const SpecialEvent* pRpSpecial;
        };

        const SplitEvent* pSplitEvent;

        unsigned int rpIdx;

        unsigned int splitIdx;

        unsigned int rpSize;

        unsigned int splitSize;
    };

    struct PrintFeatures
    {
        const char* anchorName;

        const char* spRefName;

        char strand;

        unsigned int pos;

        unsigned int len;

        unsigned int polyALen;

        unsigned int rpFrag[2];

        unsigned int splitFrag[2];

        unsigned int pos5[2];

        unsigned int pos3[2];
    };

    class Printer
    {
        public:
            Printer(const Detector* pDetector, const Aligner* pAligner, const Reference* pRef, const LibTable& libTable, const BamPairTable& bamPairTable);
            ~Printer();

            void Print();

        private:

            void PrintSpecial(void);

            void PrintSpecialHeader(void);

            void PrintSpecialBody(void);

            inline void InitFeatures(void)
            {
                features.anchorName = NULL;
                features.spRefName = NULL;

                features.strand = '/';

                features.pos = 0;
                features.len = 0;

                features.rpFrag[0] = 0;
                features.rpFrag[1] = 0;

                features.splitFrag[0] = 0;
                features.splitFrag[1] = 0;
            };

            inline void InitPrintIdx(unsigned int rpSize, unsigned int splitSize)
            {
                printIdx.pRpSpecial = NULL;
                printIdx.pSplitEvent = NULL;

                printIdx.rpIdx = 0;
                printIdx.splitIdx = 0;

                printIdx.rpSize = rpSize;
                printIdx.splitSize = splitSize;
            }

            bool GetNextSpecial(const SpecialEvent* pRpSpecials, const SplitEvent* pSplitSpecials);

            void SetSampleInfoRpSpecial(const SpecialEvent& rpSpecial);

            void SetSampleInfoSplit(const SplitEvent& splitEvent);

            // void SetSampleString(void);

            void SetSpecialFeatures(unsigned int zaSpRefID);

            void SetSpecialFeaturesFromSplit(const SplitEvent& splitEvent);

            void SetSpecialFeaturesFromRp(const SpecialEvent& rpSpecial, unsigned int zaSpRefID);

            bool SpecialPrintFilter(void);

        private:

            PrintIdx printIdx;

            PrintFeatures features;

            // std::map<unsigned int, unsigned int> sampleMap;
            
            Array<unsigned int> sampleMap;
            
            std::ostringstream formatted;

            const Detector* pDetector;

            const Aligner* pAligner;

            const Reference* pRef;

            const LibTable& libTable;

            const BamPairTable& bamPairTable;
    };
};

#endif  /*TGM_PRINTER_H*/
