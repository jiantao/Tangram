/*
 * =====================================================================================
 *
 *       Filename:  TGM_SecondMapThread.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  08/01/2012 04:44:21 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef  TGM_SECONDMAPTHREAD_H
#define  TGM_SECONDMAPTHREAD_H

#include <pthread.h>
#include <vector>

#include "../OutSources/stripedSW/ssw.h"
#include "TGM_BamPair.h"
#include "TGM_Detector.h"
#include "TGM_Parameters.h"
#include "TGM_LibTable.h"
#include "TGM_Reference.h"
#include "TGM_Sequence.h"
#include "TGM_SplitData.h"
#include "TGM_RescuePartial.h"

namespace Tangram
{
    typedef class SecondMapThread SecondMapThread;

    typedef bool (SecondMapThread::*SecondFilter)(bool& isRescued, RescuePartial& rescuePartial, uint8_t& polyALen, const s_align* pAlignment, uint8_t isReversed, 
                                                  const PrtlAlgnmnt& firstPartial, const RefRegion& refRegion, const int8_t* readSeq, int readLen);

    typedef struct
    {
        unsigned int idx;

        pthread_t thread;

    }SecondMapTag;

    typedef struct
    {
        SecondMapTag* pTags;

        SecondMapThread* pSecondMap;

        unsigned int currIdx;

        unsigned int totalWorks;

        pthread_mutex_t mutex;

    }SecondMapData;

    typedef struct
    {
        SV_EventType svType;

        int32_t pos;

        unsigned int origIdx;

    }EventTry;

    class SecondMapThread
    {
        public:

            SecondMapThread(Array<SplitEvent>& splitEvents, const AlignerPars& pars, const LibTable& libInfoTable, 
                            const BamPairTable& pairTable, const Reference& reference);

            ~SecondMapThread();

            static void* StartThread(void* threadData);

        private:

            void InitSecondPartial(SplitEvent& splitEvent);

            void SearchEventTries(Array<EventTry>& eventTries, const SplitEvent& splitEvent);

            bool TrySpecial(SplitEvent& splitEvent, TGM_Sequence& minusSeq, RescuePartial& rescuePartial);

            s_align* AlignSecPartial(bool& isRescued, RescuePartial& rescuePartial, bool& doOtherFirst, uint8_t& isReversed, uint8_t& polyALen,
                                     const PrtlAlgnmnt& firstPartial, const RefRegion& refRegion, TGM_Sequence& minusSeq, SecondFilter secondFilter);

            bool SecondFilterSpecial(bool& isRescued, RescuePartial& rescuePartial, uint8_t& polyALen, const s_align* pAlignment, uint8_t isReversed, 
                                     const PrtlAlgnmnt& firstPartial, const RefRegion& refRegion, const int8_t* readSeq, int readLen);

            void CleanUpSecond(SplitEvent& splitEvent);

            void UpdateSecPartial(PrtlAlgnmnt& secPartial, bool isRescued, const RescuePartial& rescuePartial, const s_align* pAlignment, 
                                  const PrtlAlgnmnt& firstPartial, uint8_t isReversed, uint8_t polyALen, SV_EventType svType);

            bool ProcessSpecial(SplitEvent& splitEvent, std::vector< std::vector<unsigned int> >& poll);

            inline void ClearPoll(std::vector< std::vector<unsigned int> >& poll) const
            {
                unsigned int size = poll.size();
                for (unsigned int i = 0; i != size; ++i)
                    poll[i].clear();
            }

            void MergeSpecial(SplitEvent& splitEvent);

        private:

            Array<SplitEvent>& splitEvents;

            const AlignerPars& alignerPars;

            const LibTable& libTable;

            const BamPairTable& bamPairTable;

            const Reference& ref;
    };
};

#endif  /*TGM_SECONDMAPTHREAD_H*/
