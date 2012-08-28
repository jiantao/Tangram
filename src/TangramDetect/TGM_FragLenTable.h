/*
 * =====================================================================================
 *
 *       Filename:  TGM_FragLenTable.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  05/09/2012 04:08:52 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef  TGM_FRAGLENTABLE_H
#define  TGM_FRAGLENTABLE_H

#include <vector>
#include <stdint.h>
#include "TGM_Array.h"

#define INVALID_FRAG_LEN_QUAL 255

inline int CompareFragLen(const void* a, const void* b)
{
    const uint32_t* pFragLen1 = (const uint32_t*) a;
    const uint32_t* pFragLen2 = (const uint32_t*) b;
    if (*pFragLen1 < *pFragLen2)
        return -1;
    else if (*pFragLen1 > *pFragLen2)
        return 1;
    else
        return 0;
}

namespace Tangram
{
    struct FragLenHist
    {
        uint32_t* fragLen;

        uint64_t* freq;

        uint32_t size;
    };

    class FragLenTable
    {
        public:
            FragLenTable();
            ~FragLenTable();

            void Read(FILE* fpHistInput);

            void Destory(void);

            int GetQuality(uint32_t readGrpID, uint32_t targetFragLen, uint32_t median) const;

        private:

            Array<FragLenHist> fragLenTable;
    };
};

#endif  /*TGM_FRAGLENTABLE_H*/
