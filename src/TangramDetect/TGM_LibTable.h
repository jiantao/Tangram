/*
 * =====================================================================================
 *
 *       Filename:  TGM_LibTable.h
 *
 *    Description:  Library information table
 *
 *        Created:  05/08/2012 01:14:13 PM
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Jiantao Wu (), 
 *   Inistitution:  Boston College
 * 
 * =====================================================================================
 */

#ifndef  TGM_LIBTABLE_H
#define  TGM_LIBTABLE_H

#include <map>
#include <stdint.h>
#include <cstring>
#include <vector>
#include <string>

#include "TGM_Array.h"

static const short MD5_STR_LEN = 32;

namespace Tangram
{
    typedef struct LibInfo
    {
        int32_t fragLenMedian;

        int32_t fragLenHigh;

        int32_t fragLenLow;

    }LibInfo;

    typedef enum
    {
        ST_ILLUMINA = 0,

        ST_454 = 1,

        ST_SOLID = 2,

        ST_ILLUMINA_LONG = 3

    }SeqTech;

    class LibTable
    {
        public:

            LibTable();
            ~LibTable();

            bool Read(FILE* fpLibInput, int maxFragDiff);

            bool GetReadGrpID(uint32_t& readGrpID, const char* readGrpName) const;

            bool GetSampleID(uint32_t& sampleID, const char* sampleName) const;

            bool GetSampleID(uint32_t& sampleID, uint32_t readGrpID) const;

            bool GetSpecialRefID(uint32_t& specialRefID, const char* specialRefName) const;

            inline unsigned int GetNumSamples(void) const
            {
                return sampleNames.Size();
            }

            inline const Array<char*>& GetAnchorNames(void) const
            {
                return anchorNames;
            }

            inline const Array<char*>& GetSampleNames(void) const
            {
                return sampleNames;
            }

            inline const Array<char*>* GetSpecialRefNames(void) const
            {
                return &specialRefNames;
            }

            inline const unsigned int GetNumSpecialRef(void) const
            {
                return specialRefNames.Size();
            }

            inline int32_t GetFragLenMax(void) const
            {
                return fragLenMax;
            }

            inline SeqTech GetSeqTech(uint32_t readGrpID) const
            {
                return (SeqTech) seqTech[readGrpID];
            }

            inline int32_t GetFragLenHigh(uint32_t readGrpID) const
            {
                return libInfo[readGrpID].fragLenHigh;
            }

            inline int32_t GetFragLenLow(uint32_t readGrpID) const
            {
                return libInfo[readGrpID].fragLenLow;
            }

            inline int32_t GetFragLenMedian(uint32_t readGrpID) const
            {
                return libInfo[readGrpID].fragLenMedian;
            }


        private:
            
            void ReadAnchors(uint32_t sizeAC, FILE* fpLibInput);

            void ReadSamples(uint32_t sizeSM, FILE* fpLibInput);

            void ReadReadGrps(uint32_t sizeRG, FILE* fpLibInput);

            void CheckReadGrps(int maxFragDiff);

            void ReadSpecialRef(uint32_t sizeSP, FILE* fpLibInput);

            void UpdateReadGrpHash(uint32_t readGrpID, const char* readGrpName);

            void UpdateSpecialRefHash(uint32_t sampleID, const char* readGrpName);

        private:

            Array<char*> anchorNames;

            Array<char*> specialRefNames;

            Array<char*> sampleNames;

            Array<char*> readGrpNames;

            Array<int32_t> anchorLength;

            Array<char> md5;

            std::string specialPrefix;

            Array<int32_t> sampleMap;

            Array<int8_t> seqTech;

            Array<LibInfo> libInfo;

            void* specialRefHash;

            void* readGrpHash;

            uint32_t fragLenMax;

            double cutoff;

            double trimRate;
    };
};

#endif  /*TGM_LIBTABLE_H*/

