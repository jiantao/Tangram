/*
 * =====================================================================================
 *
 *       Filename:  TGM_Detector.h
 *
 *    Description:  Detect SV using read pair method
 *
 *        Created:  05/14/2012 01:04:10 PM
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Jiantao Wu (), 
 *   Inistitution:  Boston College
 * =====================================================================================
 */

#ifndef  TGM_DETECTOR_H
#define  TGM_DETECTOR_H

#include <map>
#include <string>

#include "api/BamMultiReader.h"

#include "TGM_Parameters.h"
#include "TGM_LibTable.h"
#include "TGM_BamPair.h"
#include "TGM_PairAttrbt.h"
#include "TGM_Cluster.h"

namespace Tangram
{
    struct Inversion
    {
        int32_t refID;
        uint32_t pos;
        uint32_t posU;
        uint32_t length;
        uint32_t fragLenMax;

        int pos5[2];
        int pos3[2];

        int clusterID[2];
        uint32_t* origIndex[2];
        uint32_t numFrag[2];

        int32_t splitIdx;
    };

    struct SpecialEvent
    {
        int32_t refID;
        uint32_t pos;
        uint32_t posUncertainty;
        uint32_t length;
        uint32_t fragLenMax;

        int pos5[2];
        int pos3[2];

        int clusterID[2];
        uint32_t numFrag[2];
        uint32_t* origIndex[2];

        int32_t splitIdx;

        float sense;
    };

    static inline int CompareInversion(const void* pEvent1, const void* pEvent2)
    {
        const Inversion* pE1 = (const Inversion*) pEvent1;
        const Inversion* pE2 = (const Inversion*) pEvent2;

        if (pE1->refID < pE2->refID)
            return -1;
        else if (pE1->refID > pE2->refID)
            return 1;
        else
        {
            if (pE1->pos < pE2->pos)
                return -1;
            else if (pE1->pos > pE2->pos)
                return 1;
            else
            {
                if (pE1->length < pE2->length)
                    return -1;
                else if (pE1->length > pE2->length)
                    return 1;
                else
                    return 0;
            }
        }

        return 0;
    }

    // compare function for the special event
    static inline int CompareSpecialEvents(const void* pEvent1, const void* pEvent2)
    {
        const SpecialEvent* pE1 = (const SpecialEvent*) pEvent1;
        const SpecialEvent* pE2 = (const SpecialEvent*) pEvent2;

        if (pE1->refID < pE2->refID)
            return -1;
        else if (pE1->refID > pE2->refID)
            return 1;
        else
        {
            if (pE1->pos < pE2->pos)
                return -1;
            else if (pE1->pos > pE2->pos)
                return 1;
            else
            {
                if (pE1->length < pE2->length)
                    return -1;
                else if (pE1->length > pE2->length)
                    return 1;
                else
                    return 0;
            }
        }

        return 0;
    }

    static const double DEFAULT_MIN_STD[2] = {-1.0, -1.0};

    class Detector
    {
        public:

            Detector(const DetectPars& detectPars, const LibTable& libTable, const BamPairTable& bamPairTable);

            ~Detector();

            // initialize the read pair detector
            void Init(void);

            // call sv events
            void CallEvents(void);

        private:

            // call deletions
            void CallDeletion(void);

            // call tandem duplications
            void CallTandemDup(void);

            // call inversions
            void CallInversion(void);

            // call MEI insertions
            void CallSpecial(void);

            // call translocations
            void CallTranslocation(void);

            void MakeInversions(void);

            bool IsInvOverlapped(const Inversion& prevInv, const Inversion& newInv);

            bool IsSignificantInv(const Inversion& inv);

            void DoInversionMerge(Inversion& mergedInv, const Inversion* pHeadInv, const Inversion* pTailInv);

            void MergeInversions(void);

            // make special events from cluster information
            void MakeSpecialEvents(Array<SpecialEvent>& specialEvents);

            // merge special events
            void MergeSpecialEvents(Array<SpecialEvent>& specialEvents);

            void DoSpecialMerge(SpecialEvent* pMergedEvent, const SpecialEvent* pHeadEvent, const SpecialEvent* pTailEvent);

            void SetOrigIndices(Array<SpecialEvent>& specialEvents);

        public:

            Array<Inversion> invEvents[2];

            // special events table
            Array<SpecialEvent>* pSpecialEventsTable;

        private:

            // cluster for read pairs that the 3-prime mate is abnormal
            Cluster cluster3;

            // cluster for read pairs that the 5-prime mate is abnormal
            Cluster cluster5;

            // attribute table used for clustering
            PairAttrbtTable pairAttrbtTable;

            // parameters for read pair detector
            const DetectPars& detectPars;

            // library information table
            const LibTable& libTable;

            // bam pair table
            const BamPairTable& bamPairTable;
    };
};

#endif  /*TGM_DETECTOR_H*/
