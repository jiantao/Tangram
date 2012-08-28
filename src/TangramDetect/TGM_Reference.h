/*
 * =====================================================================================
 *
 *       Filename:  TGM_Reference.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  05/17/2012 04:05:27 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef  TGM_REFERENCE_H
#define  TGM_REFERENCE_H

#include <cstdio>
#include <zlib.h>
#include <stdint.h>
#include <vector>
#include <string>

#include "TGM_Array.h"
#include "TGM_Types.h"
#include "TGM_Sequence.h"

// the default length of padding between two special references
#define DEFAULT_PADDING_LEN 300

namespace Tangram
{
    struct RefHeader
    {
        Array<char*> names;

        Array<int64_t> endPos;

        Array<char> md5;
    };

    static inline void SeqTransfer(char* seq, unsigned int len)
    {
        for (unsigned int i = 0; i != len; ++i)
            seq[i] = nt_table[ (int) seq[i]];
    }

    class Reference
    {
        public:

            Reference();
            ~Reference();

            void Create(gzFile fpRefFastaInput, gzFile fpSpRefFastaInput, FILE* fpRefOutput);

            void Read(FILE* fpRefInput, const int32_t& refID, const int32_t& start, const int32_t& end);

            inline int32_t GetSpRefBeginPos(int32_t specialRefID) const
            {
                return (specialRefID == 0 ? 0 : spRefHeader.endPos[specialRefID - 1] + DEFAULT_PADDING_LEN + 1);
            }

            inline int32_t GetSpRefPos(int32_t specialRefID, int32_t pos) const
            {
                int32_t refBegin = GetSpRefBeginPos(specialRefID);
                return (pos - refBegin);
            }

            inline int64_t GetRefBeginPos(int32_t refID) const
            {
                return (refID == 0 ? 0 : refHeader.endPos[refID - 1] + 1);
            }

            const char* GetSpRefName(int& spRefID, int32_t spRefPos) const;

            inline int GetFamilyID(int spRefID) const
            {
                return familyMap[spRefID];
            }

            void CreateFamilyToZA(const Array<char*>& spZANames);

        private:

            void CreateRef(RefHeader& header, gzFile fpRefFastaInput, FILE* fpRefOutput, bool hasPadding);

            void InitRefHeader(RefHeader& header);

            void SetName(RefHeader& header, const char* buff, int len);

            void SetMd5(RefHeader& refHeader, char* sequence, uint32_t seqLen);

            void WriteHeader(FILE* fpRefOutput);

            void ReadRefHeader(FILE* fpRefInput);

            void ReadSpecialRef(FILE* fpRefInput);

            void ReadRef(FILE* fpRefInput, const int32_t& refID, int32_t start, int32_t end);

            void CreatFamily(void);

        public:

            Array<int8_t> refSeq;

            Array<int8_t> spRefSeq;

            uint16_t refID;

            int32_t pos;

            RefHeader refHeader;

            RefHeader spRefHeader;

            Array<int> familyMap;

            Array<int> familyToZA;

            Array<int> ZAToFamily;

            std::vector<std::string> familyName;
    };
};


#endif  /*TGM_REFERENCE_H*/
