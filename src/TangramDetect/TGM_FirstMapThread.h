/*
 * =====================================================================================
 *
 *       Filename:  TGM_FirstMapThread.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  08/01/2012 01:01:33 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef  TGM_FIRSTMAPTHREAD_H
#define  TGM_FIRSTMAPTHREAD_H

#include "../OutSources/stripedSW/ssw.h"
#include "TGM_Array.h"
#include "TGM_SplitData.h"
#include "TGM_BamPair.h"
#include "TGM_Parameters.h"
#include "TGM_LibTable.h"
#include "TGM_Reference.h"
#include "TGM_RescuePartial.h"

namespace Tangram
{
    typedef class FirstMapThread FirstMapThread;

    typedef struct
    {
        pthread_t thread;

        unsigned int idx;

        PrtlAlgnmnt* firstPartials;

        OrphanPair* orphanPairs;

        RefRegion* refRegions;

        FirstMapThread* pFirstMapThread;

        unsigned int orphanSize;

        unsigned int orphanStart;

    }FirstMapData;

    class FirstMapThread
    {
        public:

            FirstMapThread(const AlignerPars& pars, const LibTable& libTable, const BamPairTable& bamPairTable, const Reference& ref);

            ~FirstMapThread();

            static void* StartThread(void* threadData);

        private:

            bool SetFirstRefRegion(RefRegion& refRegion, int32_t pos, int32_t end, int32_t readGrpID, uint32_t readLen, bool isUp) const;
            
            bool FirstFilter(PartialType& partialType, bool& isRescued, RescuePartial& rescuePartial,
                             const s_align* pAlignment, const OrphanPair& orphanPair, const RefRegion& refRegion) const;

        private:

            const AlignerPars& alignerPars;

            const LibTable& libTable;

            const BamPairTable& bamPairTable;

            const Reference& ref;
    };
};

#endif  /*TGM_FIRSTMAPTHREAD_H*/
