/*
 * =====================================================================================
 *
 *       Filename:  TGM_BamInStream.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  08/18/2011 05:41:34 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include <stdlib.h>
#include <assert.h>

#include "khash.h"
#include "TGM_Error.h"
#include "TGM_Utilities.h"
#include "TGM_BamInStream.h"


//===============================
// Type and constant definition
//===============================

// value indicate that we did not load any read
// into the bam array
#define NO_QUERY_YET (-2)

// default capacity of a bam array
#define DEFAULT_BAM_ARRAY_CAP 200

#define DEFAULT_MATE_INFO_CAP 500

#define TGM_MAX_BIN_LEN 500000000

// a mask used to filter out those unwanted reads for split alignments
// it includes proper paired reads, secondar reads, qc-failed reads and duplicated reads
#define TGM_BAM_FMASK (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP)

// bin number of buffers
#define PREV_BIN 0
#define CURR_BIN 1

// default soft clipping tolerance
#define DEFAULT_SC_TOLERANCE 0.2

// default mismatch tolerance
#define DEFAULT_MAX_MISMATCH_RATE 0.1


// alignment status
enum AlignmentStatus
{
    NONE_GOOD   = -1,       // neither a good anchor nor a good orphan candidate

    GOOD_ANCHOR    = 0,     // a good anchor candidate

    GOOD_ORPHAN    = 1,     // a good orphan candidate

    GOOD_SOFT      = 2,     // a good soft clipping candidate

    GOOD_MULTIPLE  = 3,     // a good multiple aligned candidate
};

// initialize a string-to-bam hash table used to retrieve a pair of read
KHASH_MAP_INIT_STR(queryName, TGM_BamNode*);

KHASH_MAP_INIT_STR(name, uint32_t);

KHASH_SET_INIT_INT64(buffAddress);


//===================
// Static functions
//===================

static inline int TGM_BamInStreamLoadNext(TGM_BamInStream* pBamInStream)
{
    // for the bam alignment array, if we need to expand its space
    // we have to initialize those newly created bam alignment 
    // and update the query name hash since the address of those
    // bam alignments are changed after expanding
    pBamInStream->pNewNode = TGM_BamNodeAlloc(pBamInStream->pMemPool);
    if (pBamInStream->pNewNode == NULL)
        TGM_ErrQuit("ERROR: Too many unpaired reads are stored in the memory. Please use smaller bin size or disable searching pair genomically.\n");

    int ret = bam_read1(pBamInStream->fpBamInput, &(pBamInStream->pNewNode->alignment));

    return ret;
}


static void TGM_BamInStreamReset(TGM_BamInStream* pBamInStream)
{
    pBamInStream->pNewNode = NULL;

    pBamInStream->currBinPos = NO_QUERY_YET;
    pBamInStream->currRefID = NO_QUERY_YET;

    kh_clear(queryName, pBamInStream->pNameHashes[PREV_BIN]);
    kh_clear(queryName, pBamInStream->pNameHashes[CURR_BIN]);

    TGM_BamListReset(&(pBamInStream->pAlgnLists[PREV_BIN]), pBamInStream->pMemPool);
    TGM_BamListReset(&(pBamInStream->pAlgnLists[CURR_BIN]), pBamInStream->pMemPool);
}

static double TGM_GetMismatchRate(const bam1_t* pAlignment)
{
    uint32_t* cigar = bam1_cigar(pAlignment);
    unsigned int numMismatch = 0;
    for (unsigned i = 0; i != pAlignment->core.n_cigar; ++i)
    {
        int type = (cigar[i] & BAM_CIGAR_MASK);
        if (type == BAM_CINS || type == BAM_CDEL || type == BAM_CSOFT_CLIP || type == BAM_CMISMATCH)
        {
            numMismatch += (cigar[i] >> BAM_CIGAR_SHIFT);
        }
    }

    return ((double) numMismatch / pAlignment->core.l_qseq);
}

static int TGM_CheckAlignment(const bam1_t* pAlignment, double scTolerance, double maxMismatchRate, unsigned char minMQ)
{
    if ((pAlignment->core.flag & BAM_FUNMAP) != 0)
        return GOOD_ORPHAN;

    unsigned int scLimit = scTolerance * pAlignment->core.l_qseq;
    uint32_t* cigar = bam1_cigar(pAlignment);

    TGM_Bool isHeadSC = FALSE;
    TGM_Bool isTailSC = FALSE;

    if ((cigar[0] & BAM_CIGAR_MASK) == BAM_CSOFT_CLIP 
        && (cigar[0] >> BAM_CIGAR_SHIFT) >= scLimit)
    {
        isHeadSC = TRUE;
    }

    unsigned int lastIndex = pAlignment->core.n_cigar - 1;
    if ((cigar[lastIndex] & BAM_CIGAR_MASK) == BAM_CSOFT_CLIP 
        && (cigar[lastIndex] >> BAM_CIGAR_SHIFT) >= scLimit)
    {
        isTailSC = TRUE;
    }

    if (!isHeadSC && !isTailSC) 
    {
        if (pAlignment->core.qual >= minMQ)
            return GOOD_ANCHOR;
        else if (maxMismatchRate > 0.0)
        {
            double mismatchRate = TGM_GetMismatchRate(pAlignment);
            if (mismatchRate <= maxMismatchRate)
                return GOOD_MULTIPLE;
        }

        return NONE_GOOD;
    }
    else if (isHeadSC && isTailSC)
        return NONE_GOOD;
    else
        return GOOD_SOFT;
}

static TGM_MateInfoTable* TGM_MateInfoTableAlloc(void)
{
    TGM_MateInfoTable* pMateInfoTable = (TGM_MateInfoTable*) malloc(sizeof(TGM_MateInfoTable));
    if (pMateInfoTable == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for the mate information stack.\n");

    pMateInfoTable->capacity = DEFAULT_MATE_INFO_CAP;
    pMateInfoTable->data = (TGM_MateInfo*) calloc(sizeof(TGM_MateInfo), DEFAULT_MATE_INFO_CAP);
    if (pMateInfoTable->data == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for the mate information.\n");

    pMateInfoTable->pNameHash = kh_init(name);
    kh_resize(name, pMateInfoTable->pNameHash, DEFAULT_MATE_INFO_CAP);

    pMateInfoTable->pAvlStack = (uint32_t*) malloc(sizeof(uint32_t) * DEFAULT_MATE_INFO_CAP);
    if (pMateInfoTable->pAvlStack == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for index stack.\n");

    for (int i = 0, j = DEFAULT_MATE_INFO_CAP - 1; i != DEFAULT_MATE_INFO_CAP; ++i, --j)
        pMateInfoTable->pAvlStack[i] = j;

    pMateInfoTable->capacity = DEFAULT_MATE_INFO_CAP;
    pMateInfoTable->top = DEFAULT_MATE_INFO_CAP;

    return pMateInfoTable;
}

static void TGM_MateInfoTableFree(TGM_MateInfoTable* pMateInfoTable)
{
    if (pMateInfoTable != NULL)
    {
        for (unsigned int i = 0; i != pMateInfoTable->capacity; ++i)
            free(pMateInfoTable->data[i].queryName);

        free(pMateInfoTable->data);
        free(pMateInfoTable->pAvlStack);

        kh_destroy(name, pMateInfoTable->pNameHash);

        free(pMateInfoTable);
    }
}

static TGM_Status TGM_MateInfoTablePut(TGM_MateInfoTable* pMateInfoTable, int64_t* pIndex, const bam1_t* pAlgn)
{
    // if we do not have any available mate information object
    // we have to reallocate the mate information array
    if (pMateInfoTable->top == 0)
    {
        uint32_t oldCap = pMateInfoTable->capacity;
        pMateInfoTable->capacity *= 2;

        // reallocate the available element stack
        pMateInfoTable->pAvlStack = (uint32_t*) realloc(pMateInfoTable->pAvlStack, sizeof(uint32_t) * pMateInfoTable->capacity);
        if (pMateInfoTable->pAvlStack == NULL)
            TGM_ErrQuit("ERROR: Not enough memory for the index stack.\n");

        // reallocate the mate information array
        pMateInfoTable->data = (TGM_MateInfo*) realloc(pMateInfoTable->data, sizeof(TGM_MateInfo) * pMateInfoTable->capacity);
        if (pMateInfoTable->data == NULL)
            TGM_ErrQuit("ERROR: Not enough memory for the mate information.\n");

        // initialize the newly allocated space
        memset(pMateInfoTable->data + oldCap, 0, sizeof(TGM_MateInfo) * oldCap);

        // put the newly allocated elements into the stack
        for (unsigned int i = 0, j = pMateInfoTable->capacity - 1; i != oldCap; ++i, --j)
            pMateInfoTable->pAvlStack[i] = j;

        // set the position of the stack top
        pMateInfoTable->top = oldCap;
    }

    // try to put the query name of the incoming alignment into the name hash
    int ret = 0;
    khiter_t khIter = 0;
    khIter = kh_put(name, pMateInfoTable->pNameHash, bam1_qname(pAlgn), &ret);

    if (ret == 0)       // we found its mate
    {
        // get the index of the element in the mate information array
        *pIndex = kh_value((khash_t(name)*) pMateInfoTable->pNameHash, khIter);
        // free this element in the name hash
        kh_del(name, pMateInfoTable->pNameHash, khIter);

        // put this element back to the available stack
        pMateInfoTable->pAvlStack[pMateInfoTable->top] = *pIndex;
        ++(pMateInfoTable->top);

        return TGM_OK;
    }
    else               // we do not find its mate
    {
        // free this element in the name hash
        kh_del(name, pMateInfoTable->pNameHash, khIter);
        
        // if the upstream mate is not found it means it is filtered out
        if (pAlgn->core.pos > pAlgn->core.mpos)
            return TGM_NOT_FOUND;

        // get a mate information element to hold the information
        --(pMateInfoTable->top);
        *pIndex = pMateInfoTable->pAvlStack[pMateInfoTable->top];

        uint8_t nameLen = strlen(bam1_qname(pAlgn));
        if (pMateInfoTable->data[*pIndex].nameCap < nameLen)
        {
            free(pMateInfoTable->data[*pIndex].queryName);
            pMateInfoTable->data[*pIndex].queryName = (char*) malloc(sizeof(char) * (nameLen + 1));
            if (pMateInfoTable->data[*pIndex].queryName == NULL)
                TGM_ErrQuit("ERROR: Not enough memory for the query name.\n");
            
            pMateInfoTable->data[*pIndex].nameCap = nameLen;
        }

        // extract useful information from the alignment and store it into the element
        strcpy(pMateInfoTable->data[*pIndex].queryName, bam1_qname(pAlgn));
        pMateInfoTable->data[*pIndex].mapQ = pAlgn->core.qual;
        pMateInfoTable->data[*pIndex].numMM = TGM_GetNumMismatchFromBam(pAlgn);
        pMateInfoTable->data[*pIndex].upEnd = bam_calend(&(pAlgn->core), bam1_cigar(pAlgn));

        // re-put the query name into the name hash
        khIter = kh_put(name, pMateInfoTable->pNameHash, pMateInfoTable->data[*pIndex].queryName, &ret);
        kh_value((khash_t(name)*) pMateInfoTable->pNameHash, khIter) = *pIndex;

        return TGM_NOT_FOUND;
    }
}

static void TGM_MateInfoTableClear(TGM_MateInfoTable* pMateInfoTable)
{
    kh_clear(name, pMateInfoTable->pNameHash);

    if (pMateInfoTable->top != pMateInfoTable->capacity)
    {
        for (int i = 0, j = pMateInfoTable->capacity - 1; i != pMateInfoTable->capacity; ++i, --j)
            pMateInfoTable->pAvlStack[i] = j;

        pMateInfoTable->top = pMateInfoTable->capacity;
    }
}

static inline int TGM_BamInStreamLiteLoadNext(TGM_BamInStreamLite* pBamInStreamLite)
{
    uint8_t loadIndex = (pBamInStreamLite->head + pBamInStreamLite->size) % 4;

    // for the bam alignment array, if we need to expand its space
    // we have to initialize those newly created bam alignment 
    // and update the query name hash since the address of those
    // bam alignments are changed after expanding
    int ret = bam_read1(pBamInStreamLite->pBamInput, pBamInStreamLite->pBamBuff[loadIndex]);
    if (ret > 0)
    {
        pBamInStreamLite->tail = loadIndex;
        ++(pBamInStreamLite->size);
    }

    return ret;
}

//===============================
// Constructors and Destructors
//===============================

TGM_BamInStream* TGM_BamInStreamAlloc(uint32_t binLen, unsigned int numThreads, unsigned int buffCapacity, 
                                    unsigned int reportSize, TGM_StreamMode* pStreamMode)
{
    TGM_BamInStream* pBamInStream = (TGM_BamInStream*) calloc(1, sizeof(TGM_BamInStream));
    if (pBamInStream == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for a bam input stream object.");

    pBamInStream->fpBamInput = NULL;
    pBamInStream->pBamIndex = NULL;

    pBamInStream->filterFunc = pStreamMode->filterFunc;
    pBamInStream->filterData = pStreamMode->filterData;
    pBamInStream->controlFlag = pStreamMode->controlFlag;

    pBamInStream->numThreads = numThreads;
    pBamInStream->reportSize = reportSize;

    pBamInStream->currRefID = NO_QUERY_YET;
    pBamInStream->currBinPos = NO_QUERY_YET;
    pBamInStream->binLen = binLen;
    pBamInStream->pNewNode = NULL;

    if (numThreads > 0)
    {
        pBamInStream->pRetLists = (TGM_BamList*) calloc(numThreads, sizeof(TGM_BamList));
        if (pBamInStream->pRetLists == NULL)
            TGM_ErrQuit("ERROR: Not enough memory for the storage of retrun alignment lists in the bam input stream object.\n");

        pBamInStream->pAlgnTypes = (TGM_AlgnType*) malloc(numThreads * reportSize * sizeof(TGM_AlgnType));
        if (pBamInStream->pAlgnTypes == NULL)
            TGM_ErrQuit("ERROR: Not enough memory for the storage of pair alignment type in the bam input stream object.\n");
    }
    else
    {
        pBamInStream->pRetLists = NULL;
        pBamInStream->pAlgnTypes = NULL;
        pBamInStream->reportSize = 0;
    }

    pBamInStream->pNameHashes[PREV_BIN] = kh_init(queryName);
    if (pBamInStream->reportSize > 0)
        kh_resize(queryName, pBamInStream->pNameHashes[PREV_BIN], reportSize);

    pBamInStream->pNameHashes[CURR_BIN] = kh_init(queryName);
    if (pBamInStream->reportSize > 0)
        kh_resize(queryName, pBamInStream->pNameHashes[CURR_BIN], reportSize);

    pBamInStream->pMemPool = TGM_BamMemPoolAlloc(buffCapacity);

    return pBamInStream;
}

void TGM_BamInStreamFree(TGM_BamInStream* pBamInStream)
{
    if (pBamInStream != NULL)
    {
        kh_destroy(queryName, pBamInStream->pNameHashes[PREV_BIN]);
        kh_destroy(queryName, pBamInStream->pNameHashes[CURR_BIN]);

        free(pBamInStream->pRetLists);
        free(pBamInStream->pAlgnTypes);
        TGM_BamMemPoolFree(pBamInStream->pMemPool);

        free(pBamInStream);
    }
}

TGM_BamInStreamLite* TGM_BamInStreamLiteAlloc(void)
{
    TGM_BamInStreamLite* pBamInStreamLite = (TGM_BamInStreamLite*) calloc(sizeof(TGM_BamInStreamLite), 1);
    if (pBamInStreamLite == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for the split bam instream object.\n");

    for (unsigned int i = 0; i != 4; ++i)
        pBamInStreamLite->pBamBuff[i] = bam_init1();

    pBamInStreamLite->pMateInfoTable = NULL;
    pBamInStreamLite->currRefID = NO_QUERY_YET;

    return pBamInStreamLite;
}

void TGM_BamInStreamLiteFree(TGM_BamInStreamLite* pBamInStreamLite)
{
    if (pBamInStreamLite != NULL)
    {
        for(unsigned int i = 0; i != 4; ++i)
            bam_destroy1(pBamInStreamLite->pBamBuff[i]);

        TGM_MateInfoTableFree(pBamInStreamLite->pMateInfoTable);
        TGM_BamInStreamLiteClose(pBamInStreamLite);

        free(pBamInStreamLite);
    }
}


//======================
// Interface functions
//======================

// open a bam file
TGM_Status TGM_BamInStreamOpen(TGM_BamInStream* pBamInStream, const char* bamFilename)
{
    pBamInStream->fpBamInput = bam_open(bamFilename, "r");
    if (pBamInStream->fpBamInput == NULL)
    {
        TGM_ErrMsg("ERROR: Cannot open bam file: %s.\n", bamFilename);
        return TGM_ERR;
    }

    if ((pBamInStream->controlFlag & TGM_USE_BAM_INDEX) != 0)
    {
        pBamInStream->pBamIndex = bam_index_load(bamFilename);
        if (pBamInStream->pBamIndex == NULL)
        {
            TGM_ErrMsg("WARNING: Cannot open bam index file for: %s\n", bamFilename);
            return TGM_ERR;
        }
    }

    return TGM_OK;
}

// clear the bam instream object(return list, alignment list and name hash)
void TGM_BamInStreamClear(TGM_BamInStream* pBamInStream)
{
    pBamInStream->pNewNode = NULL;
    pBamInStream->currRefID = NO_QUERY_YET;
    pBamInStream->currBinPos = NO_QUERY_YET;

    for (unsigned int i = 0; i != pBamInStream->numThreads; ++i)
        TGM_BamListReset(pBamInStream->pRetLists + i, pBamInStream->pMemPool);

    for (unsigned int i = 0; i != 2; ++i)
        TGM_BamListReset(pBamInStream->pAlgnLists + i, pBamInStream->pMemPool);

    kh_clear(queryName, pBamInStream->pNameHashes[PREV_BIN]);
    kh_clear(queryName, pBamInStream->pNameHashes[CURR_BIN]);
}

// close the current bam files and clear the bam instream
void TGM_BamInStreamClose(TGM_BamInStream* pBamInStream)
{
    bam_close(pBamInStream->fpBamInput);
    pBamInStream->fpBamInput = NULL;

    if (pBamInStream->pBamIndex != NULL)
    {
        bam_index_destroy(pBamInStream->pBamIndex);
        pBamInStream->pBamIndex = NULL;
    }

    TGM_BamInStreamClear(pBamInStream);
}

// jump to a certain chromosome in a bam file
TGM_Status TGM_BamInStreamJump(TGM_BamInStream* pBamInStream, int32_t refID)
{
    // if we do not have the index file return error
    if (pBamInStream->pBamIndex == NULL)
        return TGM_ERR;

    // clear the bam array before jump
    TGM_BamInStreamReset(pBamInStream);

    // jump and read the first alignment in the given chromosome
    int ret;
    bam_iter_t pBamIter = bam_iter_query(pBamInStream->pBamIndex, refID, 0, INT_MAX);

    pBamInStream->pNewNode = TGM_BamNodeAlloc(pBamInStream->pMemPool);
    ret = bam_iter_read(pBamInStream->fpBamInput, pBamIter, &(pBamInStream->pNewNode->alignment));
    bam_iter_destroy(pBamIter);

    khash_t(queryName)* pNameHashCurr = NULL;

    // see if we jump to the desired chromosome
    if (ret > 0 && pBamInStream->pNewNode->alignment.core.tid == refID)
    {
        // exclude those reads who are non-paired-end, qc-fail, duplicate-marked, proper-paired, 
        // both aligned, secondary-alignment and no-name-specified.
        if ((pBamInStream->pNewNode->alignment.core.flag & TGM_BAM_FMASK) != 0
            || strcmp(bam1_qname(&(pBamInStream->pNewNode->alignment)), "*") == 0)
        {
            TGM_BamNodeFree(pBamInStream->pNewNode, pBamInStream->pMemPool);
            pBamInStream->pNewNode = NULL;
            pBamInStream->currBinPos = NO_QUERY_YET;
        }
        else
        {
            TGM_BamListPushHead(&(pBamInStream->pAlgnLists[CURR_BIN]), pBamInStream->pNewNode);

            int khRet = 0;
            khiter_t khIter = kh_put(queryName, pBamInStream->pNameHashes[CURR_BIN], bam1_qname(&(pBamInStream->pNewNode->alignment)), &khRet);

            if (khRet != 0)
            {
                pNameHashCurr = pBamInStream->pNameHashes[CURR_BIN];
                kh_value(pNameHashCurr, khIter) = pBamInStream->pNewNode;
            }
            else
                return TGM_ERR;

            pBamInStream->currBinPos = pBamInStream->pNewNode->alignment.core.pos; 
            pBamInStream->pNewNode = NULL;
        }

        pBamInStream->currRefID = refID;
        return TGM_OK;
    }
    else if (ret == -1)
    {
        return TGM_OUT_OF_RANGE;
    }
    else
    {
        return TGM_ERR;
    }
}

// read the header of a bam file
TGM_BamHeader* TGM_BamInStreamLoadHeader(TGM_BamInStream* pBamInStream)
{
    bam_header_t* pOrigHeader = bam_header_read(pBamInStream->fpBamInput);
    if (pOrigHeader == NULL)
        return NULL;

    TGM_BamHeader* pBamHeader = TGM_BamHeaderAlloc();

    pBamHeader->pOrigHeader = pOrigHeader;


    pBamHeader->pMD5s = (const char**) calloc(pOrigHeader->n_targets, sizeof(char*));
    if (pBamHeader->pMD5s == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for md5 string");

    unsigned int numMD5 = 0;
    for (const char* md5Pos = pOrigHeader->text; numMD5 <= pOrigHeader->n_targets && (md5Pos = strstr(md5Pos, "M5:")) != NULL; ++numMD5, ++md5Pos)
    {
        pBamHeader->pMD5s[numMD5] = md5Pos + 3;
    }

    if (numMD5 != pOrigHeader->n_targets)
    {
        free(pBamHeader->pMD5s);
        pBamHeader->pMD5s = NULL;

        if (numMD5 != 0)
            TGM_ErrMsg("WARNING: Number of MD5 string is not consistent with number of chromosomes.");
    }

    return pBamHeader;
}

// load a unique-orphan pair from a bam file
TGM_Status TGM_BamInStreamLoadPair(TGM_BamNode** ppUpAlgn, TGM_BamNode** ppDownAlgn, TGM_BamInStream* pBamInStream) 
{
    (*ppUpAlgn) = NULL;
    (*ppDownAlgn) = NULL;

    khash_t(queryName)* pNameHashPrev = pBamInStream->pNameHashes[PREV_BIN];
    khash_t(queryName)* pNameHashCurr = pBamInStream->pNameHashes[CURR_BIN];

    int ret = 1;
    while(ret > 0 && (ret = TGM_BamInStreamLoadNext(pBamInStream)) > 0)
    {
        // exclude those reads who are non-paired-end, qc-fail, duplicate-marked, proper-paired, 
        // both aligned, secondary-alignment and no-name-specified.
        TGM_StreamCode filterCode = pBamInStream->filterFunc(&(pBamInStream->pNewNode->alignment), pBamInStream->filterData);

        if (filterCode != STREAM_KEEP)
        {
            TGM_BamNodeFree(pBamInStream->pNewNode, pBamInStream->pMemPool);
            pBamInStream->pNewNode = NULL;

            if (filterCode == STREAM_PASS)
                continue;
            else
            {
                ret = TGM_OK;
                break;
            }
        }

        // update the current ref ID or position if the incoming alignment has a 
        // different value. The name hash and the bam array will be reset
        if (pBamInStream->pNewNode->alignment.core.pos >= pBamInStream->currBinPos + 2 * pBamInStream->binLen
            || pBamInStream->pNewNode->alignment.core.tid != pBamInStream->currRefID
            || pBamInStream->currBinPos == NO_QUERY_YET)
        {
            if (pBamInStream->pNewNode->alignment.core.tid != pBamInStream->currRefID
                && pBamInStream->currRefID != NO_QUERY_YET)
            {
                ret = TGM_OUT_OF_RANGE;
            }

            pBamInStream->currRefID  = pBamInStream->pNewNode->alignment.core.tid;
            pBamInStream->currBinPos = pBamInStream->pNewNode->alignment.core.pos;

            kh_clear(queryName, pNameHashPrev);
            kh_clear(queryName, pNameHashCurr);

            TGM_BamListReset(&(pBamInStream->pAlgnLists[PREV_BIN]), pBamInStream->pMemPool);
            TGM_BamListReset(&(pBamInStream->pAlgnLists[CURR_BIN]), pBamInStream->pMemPool);

        }
        else if (pBamInStream->pNewNode->alignment.core.pos >= pBamInStream->currBinPos + pBamInStream->binLen)
        {
            pBamInStream->currBinPos += pBamInStream->binLen;

            kh_clear(queryName, pNameHashPrev);
            TGM_SWAP(pNameHashPrev, pNameHashCurr, khash_t(queryName)*);

            TGM_BamListReset(&(pBamInStream->pAlgnLists[PREV_BIN]), pBamInStream->pMemPool);

            TGM_SWAP(pBamInStream->pAlgnLists[PREV_BIN], pBamInStream->pAlgnLists[CURR_BIN], TGM_BamList);
        }

        TGM_BamListPushHead(&(pBamInStream->pAlgnLists[CURR_BIN]), pBamInStream->pNewNode);

        (*ppUpAlgn) = NULL;
        (*ppDownAlgn) = NULL;

        khiter_t khIter = 0;

        khIter = kh_get(queryName, pNameHashPrev, bam1_qname(&(pBamInStream->pNewNode->alignment)));

        if (khIter != kh_end(pNameHashPrev))
        {
            ret = TGM_OK;
            (*ppUpAlgn) = kh_value(pNameHashPrev, khIter);
            (*ppDownAlgn) = pBamInStream->pNewNode;

            kh_del(queryName, pNameHashPrev, khIter);

            TGM_BamListRemove(&(pBamInStream->pAlgnLists[PREV_BIN]), (*ppUpAlgn));
            TGM_BamListRemove(&(pBamInStream->pAlgnLists[CURR_BIN]), (*ppDownAlgn));
        }
        else
        {
            int khRet = 0;
            khIter = kh_put(queryName, pNameHashCurr, bam1_qname(&(pBamInStream->pNewNode->alignment)), &khRet);

            if (khRet == 0) // we found a pair of alignments 
            {
                ret = TGM_OK;
                (*ppUpAlgn) = kh_value(pNameHashCurr, khIter);
                (*ppDownAlgn) = pBamInStream->pNewNode;

                kh_del(queryName, pNameHashCurr, khIter);

                TGM_BamListRemove(&(pBamInStream->pAlgnLists[CURR_BIN]), (*ppUpAlgn));
                TGM_BamListRemove(&(pBamInStream->pAlgnLists[CURR_BIN]), (*ppDownAlgn));
            }
            else // not finding corresponding mate, save the current value and move on
            {
                kh_value(pNameHashCurr, khIter) = pBamInStream->pNewNode;
            }
        }
    }

    pBamInStream->pNameHashes[PREV_BIN] = pNameHashPrev;
    pBamInStream->pNameHashes[CURR_BIN] = pNameHashCurr;

    if (ret < 0)
    {
        if (ret != TGM_OUT_OF_RANGE)
            TGM_BamNodeFree(pBamInStream->pNewNode, pBamInStream->pMemPool);

        if ( ret != TGM_OUT_OF_RANGE && ret != TGM_EOF)
            ret = TGM_ERR;
    }

    pBamInStream->pNewNode = NULL;
    return ret;
}

// get the alignment type from a read pair
TGM_AlgnType TGM_GetAlignmentType(TGM_BamNode** ppAlgnOne, TGM_BamNode** ppAlgnTwo, double scTolerance, double maxMismatchRate, unsigned char minMQ)
{
    int firstType = TGM_CheckAlignment(&((*ppAlgnOne)->alignment), scTolerance, maxMismatchRate, minMQ);
    int secondType = TGM_CheckAlignment(&((*ppAlgnTwo)->alignment), scTolerance, maxMismatchRate, minMQ);

    if (firstType != GOOD_ANCHOR && secondType == GOOD_ANCHOR)
    {
        TGM_SWAP(*ppAlgnOne, *ppAlgnTwo, TGM_BamNode*);
        TGM_SWAP(firstType, secondType, int);
    }

    if (firstType == GOOD_ANCHOR)
    {
        switch(secondType)
        {
            case GOOD_ORPHAN:
                return TGM_UNIQUE_ORPHAN;
                break;
            case GOOD_SOFT:
                return TGM_UNIQUE_SOFT;
                break;
            case GOOD_MULTIPLE:
                return TGM_UNIQUE_MULTIPLE;
                break;
            case GOOD_ANCHOR:
                return TGM_UNIQUE_NORMAL;
                break;
            default:
                return TGM_OTHER_ALGN_TYPE;
                break;
        }
    }

    return TGM_OTHER_ALGN_TYPE;
}

// load a pair of bam alignments
TGM_Status TGM_LoadAlgnPairs(TGM_BamInStream* pBamInStream, unsigned int threadID, double scTolerance, double maxMismatchRate, unsigned char minMQ)
{
    TGM_BamNode* pAlgnOne = NULL;
    TGM_BamNode* pAlgnTwo = NULL;

    TGM_Status readerStatus = TGM_OK;
    TGM_Status bufferStatus = TGM_OK;
    while ((readerStatus = TGM_BamInStreamLoadPair(&pAlgnOne, &pAlgnTwo, pBamInStream)) == TGM_OK)
    {
        TGM_AlgnType algnType = TGM_GetAlignmentType(&pAlgnOne, &pAlgnTwo, scTolerance, maxMismatchRate, minMQ);
        if ((algnType == TGM_UNIQUE_ORPHAN || algnType == TGM_UNIQUE_SOFT || algnType == TGM_UNIQUE_MULTIPLE) && pBamInStream->numThreads > 0)
        {
            TGM_BamInStreamPush(pBamInStream, pAlgnOne, threadID);
            bufferStatus = TGM_BamInStreamPush(pBamInStream, pAlgnTwo, threadID);

            TGM_BamInStreamSetAlgnType(pBamInStream, threadID, algnType);
        }
        else
        {
            TGM_BamInStreamRecycle(pBamInStream, pAlgnOne);
            TGM_BamInStreamRecycle(pBamInStream, pAlgnTwo);
        }

        if (bufferStatus == TGM_FULL)
            break;
    }

    return readerStatus;
}

// decrease the size of memory pool inside the bam in stream object to save memory
unsigned int TGM_BamInStreamShrinkPool(TGM_BamInStream* pBamInStream, unsigned int newSize)
{
    unsigned int currSize = pBamInStream->pMemPool->numBuffs;
    if (currSize > newSize)
    {
        khash_t(buffAddress)* buffHash = kh_init(buffAddress);
        kh_resize(buffAddress, buffHash, currSize);

        int ret = 0;
        khiter_t khIter = 0;
        int64_t address = 0;

        for (TGM_BamNode* pUsedNode = pBamInStream->pAlgnLists[PREV_BIN].first; pUsedNode != NULL; pUsedNode = pUsedNode->next)
        {
            address = (int64_t) pUsedNode->whereFrom;
            kh_put(buffAddress, buffHash, address, &ret);
        }

        for (TGM_BamNode* pUsedNode = pBamInStream->pAlgnLists[CURR_BIN].first; pUsedNode != NULL; pUsedNode = pUsedNode->next)
        {
            address = (int64_t) pUsedNode->whereFrom;
            kh_put(buffAddress, buffHash, address, &ret);
        }

        unsigned int delNum = currSize - newSize;
        TGM_BamBuff* pPrevBuff = NULL;
        TGM_BamBuff* pCurrBuff = pBamInStream->pMemPool->pFirstBuff;

        while (pCurrBuff != NULL && delNum != 0)
        {
            int64_t address = (int64_t) pCurrBuff;
            khIter = kh_get(buffAddress, buffHash, address);

            if (khIter == kh_end(buffHash))
            {
                TGM_BamBuff* pDelBuff = pCurrBuff;
                pCurrBuff = pCurrBuff->nextBuff;

                TGM_BamBuffClear(pDelBuff, pBamInStream->pMemPool);
                TGM_BamBuffFree(pDelBuff, pBamInStream->pMemPool->buffCapacity);
                --delNum;
                --(pBamInStream->pMemPool->numBuffs);

                if (pPrevBuff != NULL)
                    pPrevBuff->nextBuff = pCurrBuff;
                else
                    pBamInStream->pMemPool->pFirstBuff = pCurrBuff;
            }
            else
            {
                pPrevBuff = pCurrBuff;
                pCurrBuff = pCurrBuff->nextBuff;
            }
        }

        kh_destroy(buffAddress, buffHash);
    }

    return pBamInStream->pMemPool->numBuffs;
}

void TGM_BamInStreamLiteOpen(TGM_BamInStreamLite* pBamInStreamLite, const char* fileName)
{
    pBamInStreamLite->pBamInput = bam_open(fileName, "r");
    if (pBamInStreamLite->pBamInput == NULL)
        TGM_ErrQuit("ERROR: Cannot open the bam file: \"%s\"", fileName);
}

void TGM_BamInStreamLiteClose(TGM_BamInStreamLite* pBamInStreamLite)
{
    TGM_BamInStreamLiteClear(pBamInStreamLite);

    if (pBamInStreamLite->pBamInput != NULL)
    {
        bam_close(pBamInStreamLite->pBamInput);
        pBamInStreamLite->pBamInput = NULL;
    }

    pBamInStreamLite->currRefID = NO_QUERY_YET;

    if (pBamInStreamLite->pMateInfoTable != NULL)
        TGM_MateInfoTableClear(pBamInStreamLite->pMateInfoTable);
}

TGM_Status TGM_BamInStreamLiteTestZA(TGM_BamInStreamLite* pBamInStreamLite)
{
    int64_t bamPos = TGM_BamInStreamLiteTell(pBamInStreamLite);

    TGM_Status status = TGM_OK;

    bam1_t* pAlgn = bam_init1();
    int ret = bam_read1(pBamInStreamLite->pBamInput, pAlgn);

    if (ret <= 0)
        status = TGM_ERR;
    else
    {
        uint8_t* pZAtag = bam_aux_get(pAlgn, "ZA");
        if (pZAtag == NULL)
            status = TGM_NOT_FOUND;

        TGM_BamInStreamLiteSeek(pBamInStreamLite, bamPos, SEEK_SET);
    }

    bam_destroy1(pAlgn);

    return status;
}

void TGM_BamInStreamLiteClear(TGM_BamInStreamLite* pBamInStreamLite)
{
    pBamInStreamLite->head = 0;
    pBamInStreamLite->tail = 0;
    pBamInStreamLite->size = 0;

    pBamInStreamLite->currRefID = NO_QUERY_YET;
    
    if (pBamInStreamLite->pMateInfoTable != NULL)
        TGM_MateInfoTableClear(pBamInStreamLite->pMateInfoTable);
}

TGM_BamHeader* TGM_BamInStreamLiteLoadHeader(TGM_BamInStreamLite* pBamInStreamLite)
{
    bam_header_t* pOrigHeader = bam_header_read(pBamInStreamLite->pBamInput);
    if (pOrigHeader == NULL)
        return NULL;

    TGM_BamHeader* pBamHeader = TGM_BamHeaderAlloc();

    pBamHeader->pOrigHeader = pOrigHeader;

    pBamHeader->pMD5s = (const char**) calloc(pOrigHeader->n_targets, sizeof(char*));
    if (pBamHeader->pMD5s == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for md5 string");

    unsigned int numMD5 = 0;
    for (const char* md5Pos = pOrigHeader->text; numMD5 <= pOrigHeader->n_targets && (md5Pos = strstr(md5Pos, "M5:")) != NULL; ++numMD5, ++md5Pos)
    {
        pBamHeader->pMD5s[numMD5] = md5Pos + 3;
    }

    if (numMD5 != pOrigHeader->n_targets)
    {
        free(pBamHeader->pMD5s);
        pBamHeader->pMD5s = NULL;

        if (numMD5 != 0)
            TGM_ErrMsg("WARNING: Number of MD5 string is not consistent with number of chromosomes.");
    }

    return pBamHeader;
}

TGM_SortMode TGM_BamHeaderGetSortMode(const TGM_BamHeader* pBamHeader)
{
    const char* orderStr = strstr(pBamHeader->pOrigHeader->text, "SO:");
    if (orderStr == NULL)
        return TGM_ERR_SORTED;

    orderStr += 3;

    if (strncmp(orderStr, "coordinate", 10) == 0)
        return TGM_SORTED_COORDINATE_NO_ZA;
    else if (strncmp(orderStr, "unknown", 8) == 0)
        return TGM_SORTED_UNKWON;
    else if (strncmp(orderStr, "unsorted", 8) == 0)
        return TGM_UNSORTED;
    else if (strncmp(orderStr, "queryname", 9) == 0)
        return TGM_SORTED_NAME;
    else
        return TGM_ERR_SORTED;
}

TGM_Status TGM_BamInStreamLiteRead(const bam1_t* pAlgns[3], int* retNum, int64_t* pMateInfoIndex, TGM_BamInStreamLite* pBamInStreamLite)
{
    int ret = 0;
    *retNum = 0;
    *pMateInfoIndex = -1;

    if (pBamInStreamLite->sortMode == TGM_SORTED_COORDINATE_NO_ZA)
    {
        // create the mate information table if we do not have it yet
        if (pBamInStreamLite->pMateInfoTable == NULL)
            pBamInStreamLite->pMateInfoTable = TGM_MateInfoTableAlloc();

        while ((ret = TGM_BamInStreamLiteLoadNext(pBamInStreamLite)) > 0)
        {
            int tail = pBamInStreamLite->tail;
            TGM_StreamCode streamCode = STREAM_KEEP;

            // update the current refrence ID
            if (pBamInStreamLite->pBamBuff[tail]->core.tid != pBamInStreamLite->currRefID)
            {
                if (pBamInStreamLite->currRefID != NO_QUERY_YET)
                    TGM_MateInfoTableClear(pBamInStreamLite->pMateInfoTable);

                pBamInStreamLite->currRefID = pBamInStreamLite->pBamBuff[tail]->core.tid;
            }

            // for each pair of read, we only filter once
            if (pBamInStreamLite->pBamBuff[tail]->core.tid != pBamInStreamLite->pBamBuff[tail]->core.mtid)
                streamCode = STREAM_PASS;
            else if (pBamInStreamLite->pBamBuff[tail]->core.pos < pBamInStreamLite->pBamBuff[tail]->core.mpos)
                streamCode = pBamInStreamLite->filterFunc(pBamInStreamLite->pBamBuff[tail], pBamInStreamLite->filterData);

            --(pBamInStreamLite->size);
            if (streamCode == STREAM_KEEP)
            {
                TGM_Status status = TGM_MateInfoTablePut(pBamInStreamLite->pMateInfoTable, pMateInfoIndex, pBamInStreamLite->pBamBuff[tail]);
                if (status == TGM_OK)
                {
                    pAlgns[0] = pBamInStreamLite->pBamBuff[tail];
                    *retNum = 1;
                    ret = TGM_OK;
                    break;
                }
            }
        }
    }
    else if (pBamInStreamLite->sortMode == TGM_SORTED_COORDINATE_ZA)
    {
        while ((ret = TGM_BamInStreamLiteLoadNext(pBamInStreamLite)) > 0)
        {
            int tail = pBamInStreamLite->tail;
            TGM_StreamCode streamCode = pBamInStreamLite->filterFunc(pBamInStreamLite->pBamBuff[tail], pBamInStreamLite->filterData);

            if (pBamInStreamLite->pBamBuff[tail]->core.tid == pBamInStreamLite->pBamBuff[tail]->core.mtid)
            {
                if (pBamInStreamLite->pBamBuff[tail]->core.pos > pBamInStreamLite->pBamBuff[tail]->core.mpos)
                    streamCode = STREAM_PASS;
            }
            else if (pBamInStreamLite->pBamBuff[tail]->core.tid > pBamInStreamLite->pBamBuff[tail]->core.mtid)
            {
                streamCode = STREAM_PASS;
            }

            --(pBamInStreamLite->size);

            if (streamCode == STREAM_KEEP)
            {
                pAlgns[0] = pBamInStreamLite->pBamBuff[tail];
                *retNum = 1;
                ret = TGM_OK;
                break;
            }
        }
    }
    else // not sorted by coordinate and mates from the same fragment will group together
    {
        if (pBamInStreamLite->size == 0)
        {
            while ((ret = TGM_BamInStreamLiteLoadNext(pBamInStreamLite)) > 0)
            {
                int tail = pBamInStreamLite->tail;
                TGM_StreamCode streamCode = pBamInStreamLite->filterFunc(pBamInStreamLite->pBamBuff[tail], pBamInStreamLite->filterData);

                if (streamCode == STREAM_KEEP)
                    break;
                else
                    --(pBamInStreamLite->size);
            }
        }

        if (pBamInStreamLite->size == 1 && ret >= 0)
        {
            *retNum = 1;
            const char* queryName = bam1_qname(pBamInStreamLite->pBamBuff[pBamInStreamLite->tail]);

            while ((ret = TGM_BamInStreamLiteLoadNext(pBamInStreamLite)) > 0)
            {
                if (strcmp(queryName, bam1_qname(pBamInStreamLite->pBamBuff[pBamInStreamLite->tail])) == 0)
                    ++(*retNum);
                else
                {
                    TGM_StreamCode streamCode = pBamInStreamLite->filterFunc(pBamInStreamLite->pBamBuff[pBamInStreamLite->tail], pBamInStreamLite->filterData);
                    if (streamCode == STREAM_KEEP)
                    {
                        if (*retNum > 1)
                            break;
                        else
                        {
                            pBamInStreamLite->size = 1;
                            pBamInStreamLite->head = pBamInStreamLite->tail;
                            queryName = bam1_qname(pBamInStreamLite->pBamBuff[pBamInStreamLite->tail]);
                        }
                    }
                    else
                        --(pBamInStreamLite->size);
                }
            }
        }

        if (*retNum > 1)
        {
            for (unsigned int i = 0; i != *retNum; ++i)
            {
                unsigned int loadIndex = (pBamInStreamLite->head + i) % 4;
                pAlgns[i] = pBamInStreamLite->pBamBuff[loadIndex];
            }

            if (*retNum == 2 && pBamInStreamLite->sortMode == TGM_SORTED_NAME)
            {
                if (pAlgns[0]->core.tid > pAlgns[1]->core.tid)
                    TGM_SWAP(pAlgns[0], pAlgns[1], const bam1_t*);
                else if (pAlgns[0]->core.tid == pAlgns[1]->core.tid && pAlgns[0]->core.pos > pAlgns[1]->core.pos)
                    TGM_SWAP(pAlgns[0], pAlgns[1], const bam1_t*);
            }

            if (ret >= 0)
                ret = TGM_OK;
        }

        pBamInStreamLite->size -= *retNum;
        pBamInStreamLite->head = pBamInStreamLite->tail;
    }

    if (ret < 0 && ret != TGM_EOF)
        return TGM_ERR;

    return ret;
}

