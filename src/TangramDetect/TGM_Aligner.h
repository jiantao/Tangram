/*
 * =====================================================================================
 *
 *       Filename:  TGM_Aligner.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/10/2012 08:21:09 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef  TGM_ALIGNER_H
#define  TGM_ALIGNER_H

#include "../OutSources/stripedSW/ssw.h"
#include "TGM_Array.h"
#include "TGM_BamPair.h"
#include "TGM_Reference.h"
#include "TGM_Detector.h"
#include "TGM_Parameters.h"
#include "TGM_SplitData.h"
#include "TGM_SecondMapThread.h"

namespace Tangram
{
    class Aligner
    {
        typedef bool (Aligner::*SecondFilter)(const s_align* pAlignment, uint8_t isReversed, const PrtlAlgnmnt& firstPartial, const int8_t* readSeq, int readLen);

        public:

            Aligner(Detector& detector, BamPairTable& bamPairTable, const AlignerPars& alignerPars, const Reference& ref, const LibTable& libTable);

            ~Aligner();

            void Map(void);

            inline const SplitEvent* GetEventTypeStart(unsigned int& length, SV_EventType type) const
            {
                unsigned int start = type == SV_NORMAL ? 0 : svCount[type - 1];

                const SplitEvent* pSplitEvent = splitEvents.GetPointer(start);
                length = svCount[type] - start;

                if (length == 0)
                    pSplitEvent = NULL;

                return pSplitEvent;
            }

            inline const SplitEvent* GetSpecialStartFromZA(unsigned int& length, unsigned int zaSpRefID) const
            {
                int familyID = ref.ZAToFamily[zaSpRefID];
                if (familyID < 0)
                    return NULL;

                unsigned int start = svCount[SV_SPECIAL - 1];

                unsigned int spStart = familyID == 0 ? start : start + familyCount[familyID - 1];
                unsigned int spEnd = start + familyCount[familyID];

                length = spEnd - spStart;
                if (length == 0)
                    return NULL;

                return (splitEvents.GetPointer(spStart));
            }

            inline const SplitEvent* GetSpecialStartFromFamily(unsigned int& length, unsigned int familyID) const
            {
                if (familyID < 0)
                    return NULL;

                unsigned int start = svCount[SV_SPECIAL - 1];
                unsigned int spStart = familyID == 0 ? start : start + familyCount[familyID - 1];
                unsigned int spEnd = start + familyCount[familyID];

                length = spEnd - spStart;
                if (length == 0)
                    return NULL;

                return (splitEvents.GetPointer(spStart));
            }

        private:

            void FirstMap(void);

            void InsertFirstPartial(Array<PrtlAlgnmnt>& firstPartials, const s_align* pAlignment, const RefRegion& refRegion, 
                                    unsigned int idx, bool isUpStream, PartialType partialType);

            void InsertSoft(Array<PrtlAlgnmnt>& firstPartials, const SoftPair& softPair, unsigned int origIdx, unsigned int idx) const;

            void CleanFirstPartials(Array<PrtlAlgnmnt>& firstPartials);

            void MergeFirstPartials(const Array<PrtlAlgnmnt>& firstPartials);

            void InsertNewSplitEvent(const Array<PrtlAlgnmnt>& firstPartials, const SimpleRegion& simpleRegion, unsigned int count3, unsigned int count5, unsigned int currStart);

#ifdef DEBUG

            void Print(const PrtlAlgnmnt& partialAlignment, const TGM_Sequence& read);

            void PrintSplitEvent(void);

#endif

            void SecondMap(void);

            void InitSecondMapData(SecondMapData& mapData, SecondMapThread& secondMap);

            void DestroySecondMapData(SecondMapData& mapData);

            void Merge(void);

            void CountSplitEvents(void);

            void MergeSpecial(void);

            void DoMergeSpecial(Array<SpecialEvent>& rpSpecials, unsigned int currSplitIdx);

        public:

            Array<SplitEvent> splitEvents;

        private:

            Detector& detector;

            BamPairTable& bamPairTable;

            const AlignerPars& pars;

            const Reference& ref;

            const LibTable& libTable;

            Array<int> svCount;

            Array<int> familyCount;
    };
};

#endif  /*TGM_ALIGNER_H*/
