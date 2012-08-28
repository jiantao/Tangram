/*
 * =====================================================================================
 *
 *       Filename:  TGM_SplitData.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  08/01/2012 03:21:40 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef  TGM_SPLITDATA_H
#define  TGM_SPLITDATA_H

#include <stdint.h>
#include "TGM_Array.h"
#include "TGM_Types.h"
#include "TGM_Parameters.h"

namespace Tangram
{
    typedef enum
    {
        PARTIAL_UNKNOWN = -1,

        PARTIAL_3 = 0,

        PARTIAL_5 = 1

    }PartialType;

    struct RefRegion
    {
        const int8_t* pRef;

        unsigned int len;

        unsigned int start;
    };

    struct SimpleRegion
    {
        int pos;

        int end;

        unsigned int count;
    };

    struct PrtlAlgnmnt
    {
	int32_t refPos;

	int32_t refEnd;	

        uint32_t origIdx:30, isReversed:1, isSoft:1;

	uint32_t cigarLen:22, polyALen:8, partialType:1, isMajor:1;

	int64_t	readPos:24, readEnd:24, spRefID:16;

	uint32_t* cigar;	
    };

    struct SplitSpecial
    {
        int32_t pos;

        int32_t end;

        uint32_t familyID:16, polyALen:16;

        int32_t tsdPos;

        int32_t tsdLen;
    };

    struct SplitEvent
    {
        PrtlAlgnmnt* first3;

        PrtlAlgnmnt* second5;

        PrtlAlgnmnt* first5;

        PrtlAlgnmnt* second3;

        SplitSpecial* pSpecialData;

        uint32_t size3;

        uint32_t size5;

        int32_t pos5[2];

        int32_t pos3[2];

        int32_t refID;

        int32_t pos;

        int32_t rpIdx;

        uint32_t len:31, strand:1;

        SV_EventType svType;
    };

    inline PartialType GetPartialType(int readPos, int readEnd, unsigned int readLen)
    {
        if (readPos < EDGE_TOLERANCE)
            return PARTIAL_3;
        else if ((int) readLen - readEnd <= EDGE_TOLERANCE)
            return PARTIAL_5;
        else
            return PARTIAL_UNKNOWN;
    }
};

#endif  /*TGM_SPLITDATA_H*/
