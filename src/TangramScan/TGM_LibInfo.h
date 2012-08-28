/*
 * =====================================================================================
 *
 *       Filename:  TGM_LibInfo.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  02/03/2012 02:17:44 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef  TGM_LIBINFO_H
#define  TGM_LIBINFO_H

#include <stdint.h>
#include <stdio.h>

#include "TGM_Types.h"
#include "TGM_BamHeader.h"
#include "TGM_BamInStream.h"
#include "TGM_FragLenHist.h"

#define TGM_NUM_READ_PAIR_TYPES 8

#define TGM_NUM_SV_TYPES 5

#define DEFAULT_NUM_SPECIAL_REF 30

// a map used to map the pair mode into its corresponding number
// negative value means invalid mode
static const int8_t TGM_PairModeMap[16] = 
{ 
     -1, -1, 0, 1,

     -1, -1, 2, 3,

     4, 5, -1, -1,

     6, 7, -1, -1
};  

static const int8_t TGM_SeqTechMap[4][8] = 
{
    {0, 1, 2, 3, 4, 5, 6, 7},

    {2, 3, 0, 1, 5, 4, 7, 6},

    {1, 0, 3, 2, 6, 7, 4, 5},

    {3, 2, 1, 0, 7, 6, 5, 4}
};

typedef enum
{
    PT_UNKNOWN    = -1,

    PT_NORMAL     = 0,

    PT_LONG       = 1,

    PT_SHORT      = 2,

    PT_REVERSED   = 3,

    PT_INVERTED3  = 4,

    PT_INVERTED5  = 5,

    PT_SPECIAL3   = 6,

    PT_SPECIAL5   = 7,

    PT_CROSS      = 8,

}SV_ReadPairType;

static const int8_t SV_ReadPairTypeMap[2][8] = 
{
    {PT_INVERTED3, PT_NORMAL, PT_UNKNOWN, PT_INVERTED5, PT_UNKNOWN, PT_UNKNOWN, PT_REVERSED, PT_UNKNOWN},

    {PT_UNKNOWN, PT_UNKNOWN, PT_REVERSED, PT_UNKNOWN, PT_INVERTED3, PT_NORMAL, PT_UNKNOWN, PT_INVERTED5}
};

typedef enum
{
    SV_NORMAL             = 0,

    SV_DELETION           = 1,

    SV_TANDEM_DUP         = 2,

    SV_INVERSION          = 3,

    SV_SPECIAL            = 4,

    SV_INTER_CHR_TRNSLCTN = 5

}SV_EventType;

typedef enum
{
    ST_ILLUMINA = 0,

    ST_454 = 1,

    ST_SOLID = 2,

    ST_ILLUMINA_LONG = 3

}TGM_SeqTech;

// the object used to hold the basic statistics of a pair of alignments
typedef struct TGM_PairStats
{
    int32_t readGrpID;        // name of the read group

    int fragLen;              // fragment length of the pair

    TGM_PairMode pairMode;     // orientation mode of the pair

}TGM_PairStats;

typedef struct TGM_ZAtag
{
    uint8_t bestMQ[2];

    uint8_t secMQ[2];

    int16_t numMM[2];

    uint16_t numMappings[2];

    int32_t end[2];

    char spRef[2][3];

    int8_t readPairType;

}TGM_ZAtag;

typedef struct TGM_SpecialID
{
    char (*names)[3];

    void* pHash;
    
    uint32_t size;

    uint32_t capacity;

}TGM_SpecialID;

typedef struct TGM_AnchorInfo
{
    char** pAnchors;

    void* pAnchorHash;

    char* pMd5s;

    char* pSpecialPrefix;

    uint32_t specialPrefixLen;

    int32_t* pLength;

    uint32_t size;

    uint32_t capacity;

}TGM_AnchorInfo;

typedef struct TGM_SampleInfo
{
    char** pSamples;

    void* pSampleHash;

    double* pReadFraction;

    uint32_t size;

    uint32_t capacity;

}TGM_SampleInfo;

typedef struct TGM_LibInfo
{
    int32_t fragLenMedian;

    int32_t fragLenHigh;

    int32_t fragLenLow;

}TGM_LibInfo;

typedef struct TGM_LibInfoTable
{
    TGM_AnchorInfo* pAnchorInfo;

    TGM_SampleInfo* pSampleInfo;

    char** pReadGrps;

    void* pReadGrpHash;

    TGM_LibInfo* pLibInfo;

    int32_t* pSampleMap;

    int8_t* pSeqTech;

    uint32_t fragLenMax;

    uint32_t size;

    uint32_t capacity;

    double cutoff;

    double trimRate;

}TGM_LibInfoTable;


TGM_AnchorInfo* TGM_AnchorInfoAlloc(uint32_t capacity, const char* pSpecialPrefix);

void TGM_AnchorInfoFree(TGM_AnchorInfo* pInfo);

TGM_SampleInfo* TGM_SampleInfoAlloc(uint32_t capacity);

void TGM_SampleInfoFree(TGM_SampleInfo* pInfo);

TGM_SpecialID* TGM_SpecialIDAlloc(unsigned int capacity);

void TGM_SpecialIDFree(TGM_SpecialID* pSpecialID);

TGM_LibInfoTable* TGM_LibInfoTableAlloc(uint32_t capAnchor, uint32_t capSample, uint32_t capReadGrp, const char* pSpecialPrefix);

void TGM_LibInfoTableFree(TGM_LibInfoTable* pTable);

int TGM_GetNumMismatchFromBam(const bam1_t* pAlgn);

void TGM_GetNumMismatchFromZA(int16_t* pNumMM, int32_t* pLen, const char* cigarStr, unsigned int cigarLen, const char* mdStr, unsigned int mdLen);

TGM_Status TGM_LoadZAtag(TGM_ZAtag* pZAtag, const bam1_t* pAlignment);

TGM_Status TGM_LoadPairStats(TGM_PairStats* pPairStats, const bam1_t* pAlignment, const TGM_LibInfoTable* pTable);

//=================================================================
// function:
//      write the special reference name into the end of the
//      library information file
//
// args:
//      1. pSpecialTable: a pointer to a special pair table
//      2. libOutput: a file pointer to a library information file
//=================================================================
void TGM_SpecialIDWrite(const TGM_SpecialID* pSpecialID, FILE* pLibOutput);


void TGM_SpecialIDRead(TGM_SpecialID* pSpecialID, FILE* pLibInput);

//====================================================================
// function:
//      check if a pair of read is normal
//
// args:
//      1. ppUpAlgn: a pointer of the pointer to a bam node object
//                   for the alignment with smaller coordinate
//      1. ppDownAlgn: a pointer of the pointer to a bam node object
//                   for the alignment with greater coordinate
//
// return:
//      if the pair of read is normal, return TRUE; else, return
//      FALSE
//=====================================================================
TGM_Bool TGM_IsNormalPair(TGM_PairStats* pPairStats, TGM_ZAtag* pZAtag, TGM_Status* pZAstatus, const TGM_MateInfo* pMateInfo, 
                        unsigned int* pBackHistIndex, const bam1_t* pAlgns[3], int retNum, const TGM_LibInfoTable* pTable, unsigned short minMQ);

static inline void TGM_LibInfoTableSetCutoff(TGM_LibInfoTable* pTable, double cutoff)
{
    pTable->cutoff = cutoff;
}

static inline void TGM_LibInfoTableSetTrimRate(TGM_LibInfoTable* pTable, double trimRate)
{
    pTable->trimRate = trimRate;
}

TGM_Status TGM_LibInfoTableSetRG(TGM_LibInfoTable* pTable, unsigned int* oldSize, const TGM_BamHeader* pBamHeader);

TGM_Status TGM_LibInfoTableSetRGSplit(TGM_LibInfoTable* pTable, const TGM_BamHeader* pBamHeader, const char* pSpecialPrefix, unsigned int prefixLen);

TGM_Status TGM_LibInfoTableGetRGIndex(int32_t* pReadGrpIndex, const TGM_LibInfoTable* pTable, const char* pRreadGrpName);

TGM_Status TGM_LibInfoTableCheckPair(unsigned int* pHistIndex, const TGM_LibInfoTable* pTable, TGM_PairStats* pPairStats);

void TGM_LibInfoTableUpdate(TGM_LibInfoTable* pTable, const TGM_FragLenHistArray* pHistArray, unsigned int oldSize);

void TGM_LibInfoTableWrite(const TGM_LibInfoTable* pTable, TGM_Bool writeInfo, FILE* libFile);

TGM_LibInfoTable* TGM_LibInfoTableRead(FILE* libFile);

uint32_t TGM_LibInfoTableCountNormalChr(const TGM_LibInfoTable* pLibTable);

void TGM_LibInfoTableMerge(const char* workingDir);

TGM_Status TGM_LibInfoTableDoMerge(TGM_LibInfoTable* pDstLibTable, TGM_LibInfoTable* pSrcLibTable);

#endif  /*TGM_LIBINFO_H*/
