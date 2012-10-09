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

#include <set>

#include "TGM_LibTable.h"
#include "TGM_BamPair.h"
#include "TGM_Detector.h"
#include "TGM_Aligner.h"
#include "TGM_Reference.h"
#include "TGM_Genotype.h"

namespace Tangram
{
    const static unsigned int PRINT_BUFF_SIZE = 1000;

    struct PrintElmnt
    {
        union
        {
            const SpecialEvent* pRpSpecial;
        };

        const SplitEvent* pSplitEvent;

        int32_t pos;

        int16_t subsetIdx;

        SV_EventType svType;

        bool operator< (const PrintElmnt& element) const
        {
            return pos < element.pos;
        }
    };

    struct PrintSubset
    {
        union
        {
            const SpecialEvent* pRpSpecials;
        };

        const SplitEvent* pSplitEvents;

        int32_t rpIdx;

        int32_t srIdx;

        int32_t rpSize;

        int32_t srSize;

        SV_EventType svType;
    };

    struct OutputGroup
    {
        FILE* fpDel;

        FILE* fpSpecial;

        FILE* fpInv;
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
            Printer(const Detector* pDetector, const DetectPars& detectPars, const Aligner* pAligner, const Reference* pRef, 
                    const LibTable& libTable, const BamPairTable& bamPairTable, Genotype& genotype);

            ~Printer();

            void Init();

            void Print();

        private:

            inline void InitOutputGrp(void)
            {
                outputGrp.fpDel = NULL;
                outputGrp.fpInv = NULL;
                outputGrp.fpSpecial = NULL;
            }

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

            void InitPrintSubset(void);

            void InitPrintElmnts(void);

            void PrintHeader(void);

            void PrintDelHeader(void);

            void PrintDupHeader(void);

            void PrintInvHeader(void);

            void PrintSpecialHeader(void);

            void PrintTransHeader(void);

            inline void CloseOutputGrp(void)
            {
                if (outputGrp.fpDel != NULL)
                    fclose(outputGrp.fpDel);

                if (outputGrp.fpInv != NULL)
                    fclose(outputGrp.fpInv);

                if (outputGrp.fpSpecial != NULL)
                    fclose(outputGrp.fpSpecial);
            }

            // void PrintSpecial(void);

            void PrintSpecial(const PrintElmnt& element);

            void PrintGenotype(const Genotype& genotype, bool hasGenotype);

            /*  
            inline void InitPrintIdx(unsigned int rpSize, unsigned int splitSize)
            {
                printIdx.pRpSpecial = NULL;
                printIdx.pSplitEvent = NULL;

                printIdx.rpIdx = 0;
                printIdx.splitIdx = 0;

                printIdx.rpSize = rpSize;
                printIdx.splitSize = splitSize;
            }
            */

            bool GetNextPrintElmnt(PrintElmnt& element, int subsetIdx);

            bool GetNextSpecial(PrintElmnt& element, int subsetIdx);

            void SetSampleInfoRpSpecial(const SpecialEvent& rpSpecial);

            void SetSampleInfoSplit(const SplitEvent& splitEvent);

            // void SetSampleString(void);

            void SetSpecialFeatures(const PrintElmnt& element);

            void SetSpecialFeaturesFromSplit(const SplitEvent& splitEvent);

            void SetSpecialFeaturesFromRp(const SpecialEvent& rpSpecial);

            bool SpecialPrintFilter(void);

        private:

            OutputGroup outputGrp;

            PrintFeatures features;

            // std::map<unsigned int, unsigned int> sampleMap;
            
            std::multiset<PrintElmnt> printElmnts;

            Array<PrintSubset> printSubset;
            
            // Array<unsigned int> sampleMap;
            
            std::ostringstream formatted;

            const Detector* pDetector;

            const DetectPars& detectPars;

            const Aligner* pAligner;

            const Reference* pRef;

            const LibTable& libTable;

            const BamPairTable& bamPairTable;

            Genotype& genotype;
    };
};

#endif  /*TGM_PRINTER_H*/
