/*
 * =====================================================================================
 *
 *       Filename:  TGM_RescuePartial.h
 *
 *    Description:  Rescue the partial alignments with low sw score
 *
 *        Created:  08/01/2012 11:02:21 PM
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Jiantao Wu (), 
 *   Inistitution:  Boston College
 *
 * =====================================================================================
 */

#ifndef  TGM_RESCUEPARTIAL_H
#define  TGM_RESCUEPARTIAL_H

#include <stdint.h>
#include "ssw.h"
#include "TGM_Array.h"
#include "TGM_SplitData.h"

namespace Tangram
{
    class RescuePartial
    {
        public:

            RescuePartial(const AlignerPars& pars);
            ~RescuePartial();

            // resuce the partial alignments with low alignment score
            bool RescueLowScore(PartialType& partialType, const s_align* pAlignment, const int8_t* readSeq, int readLen, const RefRegion& refRegion);

            inline void Clear(void)
            {
                cigar = NULL;
                cigarLen = 0;
            }

            inline void Clean(void)
            {
                free(cigar);
                Clear();
            }

        private:

            // initialize the data member
            void Init(unsigned int readLen);

            // generate the score array based on the cigar string
            void DistributeScore(const s_align* pAlignment, const int8_t* readSeq, int readLen, const RefRegion& refRegion);

            // find the sub score array that has maximum sum
            bool FindMaxSubArray(int& start, int& end);

            // set the rescue data based on the sub mamximum array
            void SetRescueData(const s_align* pAlignment, int start, int end);

            // filter the rescued alignment
            bool RescueFilter(PartialType& partialType, int readLen);

        public:

            uint32_t* cigar;     // cigar info

            uint32_t cigarLen;   // cigar length

            int32_t refPos;      // mapping position on reference

            int32_t refEnd;      // mapping end on reference

            int32_t readPos;     // mapping position on the read

            int32_t readEnd;     // mapping end on the read

            int32_t bestScore;   // best alignment score

        private:
            
            Array<int> score;                  // score array

            Array<int> len;                    // length array

            Array<uint8_t> type;               // type array

            const AlignerPars& alignerPars;    // aligner parameters
    };
}


#endif  /*TGM_RESCUEPARTIAL_H*/
