/*
 * =====================================================================================
 *
 *       Filename:  TGM_BamInStream.h
 * *    Description:  *
 *        Version:  1.0
 *        Created:  08/18/2011 05:39:37 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef  TGM_BAMINSTREAM_H
#define  TGM_BAMINSTREAM_H


#include "bam.h"
#include "TGM_BamHeader.h"
#include "TGM_BamMemPool.h"

//===============================
// Type and constant definition
//===============================

typedef TGM_StreamCode (*TGM_BamFilter) (const bam1_t* pAlignment, void* pFilterData);

// control parameters for bam in stream
typedef enum TGM_StreamControlFlag
{
    TGM_NO_SPECIAL_CONTROL = 0,

    TGM_USE_BAM_INDEX      = 1,   // bam in stream will open the bam index file and load it into memory

    TGM_READ_PAIR_MODE     = 2    // bam in stream will find read pairs genomically

}TGM_StreamControlFlag;

typedef enum TGM_SortMode
{
    TGM_ERR_SORTED = -1,

    TGM_UNSORTED = 0,

    TGM_SORTED_COORDINATE_NO_ZA = 1,

    TGM_SORTED_COORDINATE_ZA = 2,

    TGM_SORTED_NAME = 3,

    TGM_SORTED_SPLIT = 4,

    TGM_SORTED_UNKWON = 5

}TGM_SortMode;

typedef struct TGM_StreamMode
{
    TGM_BamFilter filterFunc;             // a filter function used to skip those uninterested reads

    void* filterData;                    // parameters for the filter function

    TGM_StreamControlFlag controlFlag;    // flag used to control the bam in stream

}TGM_StreamMode;

// alignment type
typedef enum
{
    TGM_UNIQUE_ORPHAN = 0,      // one mate is uniquely aligned the other mate is unaligned

    TGM_UNIQUE_SOFT = 1,        // one mate is uniquely aligned the other mate has may soft-clipping bases

    TGM_UNIQUE_MULTIPLE = 2,    // one mate is uniquely aligned the other mate is mutiply aligned

    TGM_UNIQUE_NORMAL = 3,      // both mates are uniquely aligned

    TGM_OTHER_ALGN_TYPE = 4     // other alignment types

}TGM_AlgnType;

// iterator used to retrieve the alignments and alignment types
typedef struct TGM_BamInStreamIter
{
    TGM_BamNode* pBamNode;    // a pointer to the bam node structure

    TGM_AlgnType* pAlgnType;  // a pointer to the alignment type

}TGM_BamInStreamIter;

// private data structure that holds all bam-input-related information
typedef struct TGM_BamInStream
{
    bamFile fpBamInput;                        // file pointer to a input bam file

    bam_index_t* pBamIndex;                    // file pointer to a input bam index file

    TGM_BamFilter filterFunc;                   // customized filter function 

    void* filterData;                          // data used by the filter function

    TGM_BamMemPool* pMemPool;                   // memory pool used to allocate and recycle the bam alignments

    void* pNameHashes[2];                      // two hashes used to get a pair of alignments

    TGM_BamList* pRetLists;                     // when we find any unique-orphan pairs we push them into these lists, each thread has its own list

    TGM_AlgnType* pAlgnTypes;                   // store the alignment types of read pairs(unique-orphan, unique-multiple, unique-softclipping, unique-unique...)

    TGM_BamNode* pNewNode;                      // the just read-in bam alignment

    TGM_BamList pAlgnLists[2];                  // lists used to store those incoming alignments

    unsigned int numThreads;                   // number of threads will be used

    unsigned int reportSize;                   // number of alignments should be loaded before report

    int32_t currRefID;                         // the reference ID of the current read-in alignment

    int32_t currBinPos;                        // the start position of current bin (0-based)

    uint32_t binLen;                           // the length of the bin

    TGM_StreamControlFlag controlFlag; 

}TGM_BamInStream;

typedef struct TGM_MateInfo
{
    char* queryName;

    uint64_t upEnd:32, numMM:16, mapQ:8, nameCap:8;

}TGM_MateInfo;

typedef struct TGM_MateInfoTable
{
    TGM_MateInfo* data;

    uint32_t* pAvlStack;

    void* pNameHash;

    uint32_t top;

    uint32_t capacity;

}TGM_MateInfoTable;

typedef struct TGM_BamInStreamLite
{
    bamFile pBamInput;

    int32_t currRefID;

    TGM_BamFilter filterFunc;

    void* filterData;

    bam1_t* pBamBuff[4];

    TGM_MateInfoTable* pMateInfoTable;

    uint8_t head;

    uint8_t tail;

    uint8_t size;

    TGM_SortMode sortMode;

}TGM_BamInStreamLite;


//===============================
// Constructors and Destructors
//===============================

TGM_BamInStream* TGM_BamInStreamAlloc(uint32_t binLen,                       // search range of a pair
                                    
                                    unsigned int numThreads,               // number of threads
                                     
                                    unsigned int buffCapacity,             // the number of alignments can be stored in each chunk of the memory pool
                                    
                                    unsigned int reportSize,               // number of alignments should be cached before report
                                    
                                    TGM_StreamMode* pStreamMode);           // a structure used to control the features of the stream


void TGM_BamInStreamFree(TGM_BamInStream* pBamInStream);

TGM_BamInStreamLite* TGM_BamInStreamLiteAlloc(void);

void TGM_BamInStreamLiteFree(TGM_BamInStreamLite* pBamInStreamLite);


//======================
// Interface functions
//======================

//===============================================================
// function:
//      jump to a certain chromosome in a bam file
//
// args:
//      1. pBamInStream: a pointer to an bam instream structure
//      2. refID : the reference ID we want to jump to
// 
// return:
//      if jumping succeeds, return TGM_OK; if not, return TGM_ERR
//=============================================================== 
TGM_Status TGM_BamInStreamJump(TGM_BamInStream* pBamInStream, int32_t refID);

//===============================================================
// function:
//      open a bam file
//
// args:
//      1. pBamInStream: a pointer to an bam instream structure
//      2. bamFileName: the name of the bam file
// 
// return:
//      if open succeeds, return TGM_OK; if not, return TGM_ERR
//=============================================================== 
TGM_Status TGM_BamInStreamOpen(TGM_BamInStream* pBamInStream, const char* bamFileName);

//===============================================================
// function:
//      clear the bam instream object(return list, 
//      alignment list and name hash)
//
// args:
//      1. pBamInStream: a pointer to an bam instream structure
//=============================================================== 
void TGM_BamInStreamClear(TGM_BamInStream* pBamInStream);

//===============================================================
// function:
//      close the current bam files and clear the bam instream
//
// args:
//      1. pBamInStream: a pointer to an bam instream structure
//=============================================================== 
void TGM_BamInStreamClose(TGM_BamInStream* pBamInStream);

//===============================================================
// function:
//      tell the current position in the bam file
//
// args:
//      1. pBamInStream: a pointer to an bam instream structure
//
// return:
//      the current position in the bam file
//=============================================================== 
#define TGM_BamInStreamTell(pBamInStream) bam_tell((pBamInStream)->fpBamInput)

//===============================================================
// function:
//      seek to a specific position in the bam file
//
// args:
//      1. pBamInStream: a pointer to an bam instream structure
//      2. pos: file position
//      3. where: from where we should seek 
//                (only SEEK_SET is allowed)
//=============================================================== 
#define TGM_BamInStreamSeek(pBamInStream, pos, where) bam_seek((pBamInStream)->fpBamInput, pos, where)

//================================================================
// function:
//      read the header of a bam file and load necessary
//      information from the header text
//
// args:
//      1. pBamInStream: a pointer to an bam instream structure
// 
// return:
//      if reading succeeds, return a pointer to a header
//      structure; if error is found, return NULL
//
// discussion: 
//      The file position indicator must be placed at the 
//      beginning of the file. Upon success, the position 
//      indicator will be set at the start of the first alignment
//================================================================ 
TGM_BamHeader* TGM_BamInStreamLoadHeader(TGM_BamInStream* pBamInStream);

//===============================================================
// function:
//      set the mode of the bam in stream
//
// args:
//      1. pBamInStream: a pointer to an bam instream structure
//      2. filterFunc: filter function for the bam in stream
//      3. filterData: filter data passed to filter function
//      4. controlFlag: flags used to control the stream
// 
//=============================================================== 
static inline void TGM_SetStreamMode(TGM_StreamMode* pStreamMode, TGM_BamFilter filterFunc, void* filterData, TGM_StreamControlFlag controlFlag)
{
    pStreamMode->filterFunc = filterFunc;
    pStreamMode->filterData = filterData;
    pStreamMode->controlFlag = controlFlag;
}

//================================================================
// function:
//      read an alignment from the bam file
//
// args:
//      1. pAlignment: a pointer to an alignment
//      2. pBamInStream: a pointer to an bam instream structure
// 
// return:
//      if reading succeeds, return TGM_OK; if reach the end of
//      file, return TGM_EOF; if get error, return TGM_ERR
//
// discussion: 
//      The file position indicator must be placed at the 
//      beginning of an alignment. Upon success, the position 
//      indicator will be set at the start of the next alignment
//================================================================ 
static inline TGM_Status TGM_BamInStreamRead(bam1_t* pAlignment, TGM_BamInStream* pBamInStream)
{
    int ret = bam_read1(pBamInStream->fpBamInput, pAlignment);

    if (ret > 0)
        return TGM_OK;
    else if (ret == -1)
        return TGM_EOF;
    else
        return TGM_ERR;
}

//================================================================
// function:
//      get the current reference ID
//
// args:
//      1. pBamInStream: a pointer to an bam instream structure
// 
// return:
//      reference ID
//================================================================ 
#define TGM_BamInStreamGetRefID(pBamInStream) ((pBamInStream)->currRefID)

//==================================================================
// function:
//      load a pair of bam alignments
//
// args:
//      1. ppAlgnOne: a pointer to the pointer of an alignment
//      2. ppAlgnTwo: a pointer to the pointer of an alignment
//      3. pBamInStream : a pointer to an bam instream structure
//
// return:
//      if we get enough unique-orphan pair, return TGM_OK; 
//      if we reach the end of file, return TGM_EOF; if we finish 
//      the current chromosome, return TGM_OUT_OF_RANGE; 
//      else, return TGM_ERR
//==================================================================
TGM_Status TGM_BamInStreamLoadPair(TGM_BamNode** ppAlgnOne, TGM_BamNode** ppAlgnTwo, TGM_BamInStream* pBamInStream);

//==================================================================
// function:
//      get the alignment type from a read pair
//
// args:
//      ppAnchor: a pointer to a pointer of bam node with anchor
//                alignment
//      ppOrphan: a pointer to a pointer of bam node with orphan
//                alignment
//      scTolerance: soft clipping tolerance
//
// return:
//      the alignment type of the read pair
//==================================================================
TGM_AlgnType TGM_GetAlignmentType(TGM_BamNode** ppAlgnOne, TGM_BamNode** ppAlgnTwo, double scTolerance, double maxMismatchRate, unsigned char minMQ);

//==================================================================
// function:
//      load a certain number of unique orphan pairs into a buffer
//      associated with a thread
//
// args:
//      1. pBamInStream: a pointer to a bam in stream object
//      2. threadID: the ID of the thread
//      3. scTolerance: soft clipping tolerance
//
// return:
//      status of the bam in stream. if we reach the end of a
//      chromosome, return TGM_OUT_OF_RANGE; if we reach the end of
//      the file, return TGM_EOF; if an error happens, return
//      TGM_ERR; else return TGM_OK
//==================================================================
TGM_Status TGM_LoadAlgnPairs(TGM_BamInStream* pBamInStream, unsigned int threadID, double scTolerance, double maxMismatchRate, unsigned char minMQ);

//================================================================
// function:
//      get the size of the memory pool in the bam in stream
//      structure
//
// args:
//      1. pBamInStream: a pointer to an bam instream structure
// 
// return:
//      the size of memory pool in the bam in stream object
//================================================================ 
#define TGM_BamInStreamGetPoolSize(pBamInStream) ((pBamInStream)->pMemPool->numBuffs)

//================================================================
// function:
//      get a iterator to a certain buffer of a thread
//
// args:
//      1. pBamInStream: a pointer to an bam instream structure
//      2. threadID: ID of a thread
// 
// return:
//      iterator to the buffer of a thread
//================================================================ 
#define TGM_BamInStreamSetIter(pIter, pBamInStream, threadID)                                           \
    do                                                                                                 \
    {                                                                                                  \
        (pIter)->pBamNode = TGM_BamListGetIter((pBamInStream)->pRetLists + (threadID));                 \
        (pIter)->pAlgnType = (pBamInStream)->pAlgnTypes + (pBamInStream)->reportSize * (threadID);     \
                                                                                                       \
    }while(0)

#define TGM_BamInStreamSetAlgnType(pBamInStream, threadID, algnType)                                                            \
    do                                                                                                                         \
    {                                                                                                                          \
        unsigned int index = (pBamInStream)->reportSize * (threadID) + (pBamInStream)->pRetLists[(threadID)].numNode / 2 -1;   \
        (pBamInStream)->pAlgnTypes[index] = (algnType);                                                                        \
                                                                                                                               \
    }while(0)

//================================================================
// function:
//      push the qualified alignment into a given thread buffer
//
// args:
//      1. pBamInStream: a pointer to an bam instream structure
//      2. pAlignment: a pointer to a qualified alignment
//      2. threadID: ID of a thread
// 
// return:
//      status of the thread buffer. if the thread buffer is full
//      then return TGM_FULL else return TGM_OK
//================================================================ 
static inline TGM_Status TGM_BamInStreamPush(TGM_BamInStream* pBamInStream, TGM_BamNode* pAlignment, unsigned int threadID)
{
    TGM_BamListPushBack(pBamInStream->pRetLists + threadID, pAlignment);

    if (pBamInStream->pRetLists[threadID].numNode == pBamInStream->reportSize)
        return TGM_FULL;

    return TGM_OK;
}

//================================================================
// function:
//      recycle an unwanted bam node
//
// args:
//      1. pBamInStream: a pointer to an bam instream structure
//      2. pBamNode: a pointer to a bam node 
//================================================================ 
#define TGM_BamInStreamRecycle(pBamInStream, pBamNode) TGM_BamListPushHead(&((pBamInStream)->pMemPool->avlNodeList), (pBamNode))

//================================================================
// function:
//      clear the return list  associated with a certain thread
//
// args:
//      1. pBamInStream: a pointer to an bam instream structure
//      2. threadID: the id of the thread that needs to be clear
//================================================================ 
#define TGM_BamInStreamClearRetList(pBamInStream, threadID) TGM_BamListReset((pBamInStream)->pRetLists + (threadID), (pBamInStream)->pMemPool)

//================================================================
// function:
//      decrease the size of memory pool inside the bam in stream
//      object to save memory
//
// args:
//      1. pBamInStream : a pointer to an bam instream structure
//      2. newSize: the new size of the memory pool that user 
//                  want to set (which can be only smaller than
//                  current size otherwise nothing will be done)
//
// return:
//      the actual size of the memory pool after shrinking
//
// discussion:
//      this function should be called after the processing of 
//      a certain chromosome. The memory allocated for the 
//      return lists should be freed before you can call this
//      function. the desired size may not be achieved since
//      there may not be enough empty buffer chunks. check the
//      return value for the actual size of the memory pool
//================================================================
unsigned int TGM_BamInStreamShrinkPool(TGM_BamInStream* pBamInStream, unsigned int newSize);

TGM_SortMode TGM_BamHeaderGetSortMode(const TGM_BamHeader* pBamHeader);

static inline void TGM_BamInStreamLiteSetSortMode(TGM_BamInStreamLite* pBamInStreamLite, TGM_SortMode sortMode)
{
    pBamInStreamLite->sortMode = sortMode;
}

static inline void TGM_BamInStreamLiteSetFilter(TGM_BamInStreamLite* pBamInStreamLite, TGM_BamFilter filterFunc)
{
    pBamInStreamLite->filterFunc = filterFunc;
}

static inline void TGM_BamInStreamLiteSetFilterData(TGM_BamInStreamLite* pBamInStreamLite, void* filterData)
{
    pBamInStreamLite->filterData = filterData;
}

void TGM_BamInStreamLiteOpen(TGM_BamInStreamLite* pBamInStreamLite, const char* fileName);

void TGM_BamInStreamLiteClose(TGM_BamInStreamLite* pBamInStreamLite);

TGM_Status TGM_BamInStreamLiteTestZA(TGM_BamInStreamLite* pBamInStreamLite);

void TGM_BamInStreamLiteClear(TGM_BamInStreamLite* pBamInStreamLite);

#define TGM_BamInStreamLiteTell(pBamInStreamLite) bam_tell((pBamInStreamLite)->pBamInput)

#define TGM_BamInStreamLiteSeek(pBamInStreamLite, pos, where) bam_seek((pBamInStreamLite)->pBamInput, pos, where)

TGM_BamHeader* TGM_BamInStreamLiteLoadHeader(TGM_BamInStreamLite* pBamInStreamLite);

TGM_Status TGM_BamInStreamLiteRead(const bam1_t* pAlgns[3], int* retNum, int64_t* pMateInfoIndex, TGM_BamInStreamLite* pBamInStreamLite);

static inline const TGM_MateInfo* TGM_BamInStreamLiteGetMateInfo(const TGM_BamInStreamLite* pBamInStreamLite, int64_t index)
{
    if (pBamInStreamLite->pMateInfoTable != NULL && index >= 0)
        return (pBamInStreamLite->pMateInfoTable->data + index);
    else
        return NULL;
}

#endif  /*TGM_BAMINSTREAM_H*/
