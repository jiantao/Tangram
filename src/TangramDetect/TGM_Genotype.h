/*
 * =====================================================================================
 *
 *       Filename:  TGM_Genotype.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  09/21/2012 11:55:55 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef  TGM_GENOTYPE_H
#define  TGM_GENOTYPE_H

#include "TGM_Array.h"
#include "TGM_Detector.h"
#include "TGM_SplitData.h"
#include "TGM_LibTable.h"

#include "api/BamMultiReader.h"

namespace Tangram
{
    struct FragCount
    {
        uint32_t nonSupport;

        uint32_t support;
    };

    class Genotype
    {
        public:

            Genotype(BamTools::BamMultiReader& reader, const GenotypePars& genotypePars, const LibTable& libTable, BamPairTable& bamPairTable);

            ~Genotype();

            void Init(void);

            void SetSpecialPrior(const double* prior);

            bool Special(const SpecialEvent* pRpSepcial, const SplitEvent* pSplitEvent);

        private:

            bool SpecialFilter(const SpecialEvent* pRpSpecial, const SplitEvent* pSplitEvent) const;

            bool Jump(int32_t refID, int32_t pos);

            void SetSampleCountSpecial(const SpecialEvent& rpSpecial);

            void SetSampleCountSplit(const SplitEvent& splitEvent);

            inline void UpdateNonSupport(int32_t readGrpID)
            {
                unsigned int sampleID = 0;
                if (libTable.GetSampleID(sampleID, readGrpID))
                    ++(sampleCount[sampleID].nonSupport);
            }

            inline void UpdateSupport(int32_t readGrpID)
            {
                unsigned int sampleID = 0;
                if (libTable.GetSampleID(sampleID, readGrpID))
                    ++(sampleCount[sampleID].support);
            }

            void SetLikelihood(void);

            long double CalculateLikelihood(unsigned int refCount, unsigned int altCount, double p) const;


        public:

            std::vector<int8_t> genotypes;

            std::vector<double> likelihoods;

        private:

            BamTools::BamMultiReader& reader;

            const GenotypePars& genotypePars;

            const LibTable& libTable;

            BamPairTable& bamPairTable;

            int32_t lastChr;

            int32_t lastPos;

            int32_t lastEnd;

            double specialPrior[3];

            Array<FragCount> sampleCount;
    };
};

#endif  /*TGM_GENOTYPE_H*/
