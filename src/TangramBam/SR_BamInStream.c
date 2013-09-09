/*
 * =====================================================================================
 *
 *       Filename:  SR_BamInStream.c
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
#include <limits.h>

#include "../OutSources/samtools/khash.h"
#include "SR_Error.h"
#include "SR_Utilities.h"
#include "SR_BamInStream.h"

//===============================
// Type and constant definition
//===============================

// value indicate that we did not load any read
// into the bam array
static const int NO_QUERY_YET = -2;

// default capacity of a bam array
//#define DEFAULT_BAM_ARRAY_CAP 200

static const int SR_MAX_BIN_LEN = 500000000;

// a mask used to filter out those unwanted reads for split alignments
// it includes proper paired reads, secondar reads, qc-failed reads and duplicated reads
#define SR_BAM_FMASK (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP)

// bin number of buffers
static const int PREV_BIN = 0;
static const int CURR_BIN = 1;

// initialize a string-to-bam hash table used to retrieve a pair of read
KHASH_MAP_INIT_STR(queryName, SR_BamNode*);

KHASH_SET_INIT_INT64(buffAddress);

/*  
// private data structure that holds all bam-input-related information
struct SR_BamInStreamPrvt
{
    bamFile fpBamInput;                        // file pointer to a input bam file

    bam_index_t* pBamIndex;                    // file pointer to a input bam index file

    SR_BamMemPool* pMemPool;                   // memory pool used to allocate and recycle the bam alignments

    khash_t(queryName)* pNameHashes[2];        // two hashes used to get a pair of alignments

    SR_BamList* pRetLists;                     // when we find any unique-orphan pairs we push them into these lists

    SR_BamNode* pNewNode;                      // the just read-in bam alignment

    SR_BamList pAlgnLists[2];                  // lists used to store those incoming alignments. each thread has their own lists.

    unsigned int numThreads;                   // number of threads will be used

    unsigned int reportSize;                   // number of pairs should be loaded before report

    int32_t currRefID;                         // the reference ID of the current read-in alignment

    int32_t currBinPos;                        // the start position of current bin (0-based)

    uint32_t binLen;                           // the length of bin

    double scTolerance;                        // soft clipping tolerance rate. 
};

*/


//===================
// Static functions
//===================

// Read the next bam record from the bam file and store it in pBamInStream->pNewNode
static inline int SR_BamInStreamLoadNext(SR_BamInStream* pBamInStream)
{
    if (pBamInStream->bam_cur_status < 0) return -1;

    // for the bam alignment array, if we need to expand its space
    // we have to initialize those newly created bam alignment 
    // and update the query name hash since the address of those
    // bam alignments are changed after expanding
    pBamInStream->pNewNode = SR_BamNodeAlloc(pBamInStream->pMemPool);
    if (pBamInStream->pNewNode == NULL)
        SR_ErrQuit("ERROR: Too many unpaired reads are stored in the memory. Please use smaller bin size or disable searching pair genomically.\n");
    
    int ret;
    if (pBamInStream->pBamIterator != NULL)
      ret = bam_iter_read(pBamInStream->fpBamInput, *(pBamInStream->pBamIterator), &(pBamInStream->pNewNode->alignment));
    else
      ret = bam_read1(pBamInStream->fpBamInput, &(pBamInStream->pNewNode->alignment));

    pBamInStream->bam_cur_status = ret;

    return ret;
}

static void SR_BamInStreamReset(SR_BamInStream* pBamInStream)
{
    pBamInStream->pNewNode = NULL;

    if (pBamInStream->pBamIterator != NULL) {
	bam_iter_destroy(*(pBamInStream->pBamIterator));
	free(pBamInStream->pBamIterator);
	pBamInStream->pBamIterator = NULL;
    }

    pBamInStream->currBinPos = NO_QUERY_YET;
    pBamInStream->currRefID = NO_QUERY_YET;

    kh_clear(queryName, pBamInStream->pNameHashes[PREV_BIN]);
    kh_clear(queryName, pBamInStream->pNameHashes[CURR_BIN]);

    SR_BamListReset(&(pBamInStream->pAlgnLists[PREV_BIN]), pBamInStream->pMemPool);
    SR_BamListReset(&(pBamInStream->pAlgnLists[CURR_BIN]), pBamInStream->pMemPool);
}


//===============================
// Constructors and Destructors
//===============================

SR_BamInStream* SR_BamInStreamAlloc(const char* bamFilename, uint32_t binLen, unsigned int numThreads, unsigned int buffCapacity, 
                                    unsigned int reportSize, const SR_StreamMode* pStreamMode)
{
    SR_BamInStream* pBamInStream = (SR_BamInStream*) calloc(1, sizeof(SR_BamInStream));
    if (pBamInStream == NULL)
        SR_ErrQuit("ERROR: Not enough memory for a bam input stream object.");

    pBamInStream->bam_cur_status = -1;

    pBamInStream->fpBamInput = bam_open(bamFilename, "r");
    if (pBamInStream->fpBamInput == NULL)
        SR_ErrQuit("ERROR: Cannot open bam file %s for reading.\n", bamFilename);

    if ((pStreamMode->controlFlag & SR_USE_BAM_INDEX) != 0)
    {
        pBamInStream->pBamIndex = bam_index_load(bamFilename);
	if (pBamInStream->pBamIndex == NULL) {
            SR_ErrMsg("WARNING: Cannot open bam index file for reading. Creating it......");
	    bam_index_build(bamFilename);
	    SR_ErrMsg("         The bam index is created.");
	    pBamInStream->pBamIndex = bam_index_load(bamFilename);
	}
    }

    pBamInStream->filterFunc = pStreamMode->filterFunc;
    pBamInStream->filterData = pStreamMode->filterData;
    pBamInStream->numThreads = numThreads;
    pBamInStream->reportSize = reportSize;
    pBamInStream->currRefID = NO_QUERY_YET;
    pBamInStream->currBinPos = NO_QUERY_YET;
    pBamInStream->binLen = binLen;
    pBamInStream->pNewNode = NULL;
    pBamInStream->pBamIterator = NULL;

    if (numThreads > 0)
    {
        pBamInStream->pRetLists = (SR_BamList*) calloc(numThreads, sizeof(SR_BamList));
        if (pBamInStream->pRetLists == NULL)
            SR_ErrQuit("ERROR: Not enough memory for the storage of retrun alignment lists in the bam input stream object.\n");

        pBamInStream->pAlgnTypes = (SR_AlgnType*) malloc(numThreads * reportSize * sizeof(SR_AlgnType));
        if (pBamInStream->pAlgnTypes == NULL)
            SR_ErrQuit("ERROR: Not enough memory for the storage of pair alignment type in the bam input stream object.\n");
    }
    else
    {
        pBamInStream->pRetLists = NULL;
        pBamInStream->pAlgnTypes = NULL;
        pBamInStream->reportSize = 0;
    }

    if ((pStreamMode->controlFlag & SR_PAIR_GENOMICALLY) == 0)
    {
        pBamInStream->pNameHashes[PREV_BIN] = kh_init(queryName);
        kh_resize(queryName, pBamInStream->pNameHashes[PREV_BIN], reportSize);
    }
    else
    {
        pBamInStream->pNameHashes[PREV_BIN] = NULL;
        pBamInStream->binLen = SR_MAX_BIN_LEN;
    }

    pBamInStream->pNameHashes[CURR_BIN] = kh_init(queryName);
    kh_resize(queryName, pBamInStream->pNameHashes[CURR_BIN], reportSize);

    pBamInStream->pMemPool = SR_BamMemPoolAlloc(buffCapacity);

    pBamInStream->bam_cur_status = 1;

    return pBamInStream;
}

void SR_BamInStreamFree(SR_BamInStream* pBamInStream)
{
    if (pBamInStream != NULL)
    {
        kh_destroy(queryName, pBamInStream->pNameHashes[PREV_BIN]);
        kh_destroy(queryName, pBamInStream->pNameHashes[CURR_BIN]);

        if (pBamInStream->pRetLists != NULL)
	  free(pBamInStream->pRetLists);
        if (pBamInStream->pAlgnTypes != NULL)
	  free(pBamInStream->pAlgnTypes);
        SR_BamMemPoolFree(pBamInStream->pMemPool);

        bam_close(pBamInStream->fpBamInput);
        bam_index_destroy(pBamInStream->pBamIndex);

	if (pBamInStream->pBamIterator != NULL) {
	  bam_iter_destroy(*(pBamInStream->pBamIterator));
	  free(pBamInStream->pBamIterator);
	  pBamInStream->pBamIterator = NULL;
	}

        free(pBamInStream);
    }
}


//======================
// Interface functions
//======================

// jump to a certain chromosome in a bam file
SR_Status SR_BamInStreamJump(SR_BamInStream* pBamInStream, int32_t refID, int32_t begin, int32_t end)
{
    // if we do not have the index file return error
    if (pBamInStream->pBamIndex == NULL)
        return SR_ERR;
    
    if (pBamInStream->pBamIterator != NULL) {
      bam_iter_destroy(*(pBamInStream->pBamIterator));
      free(pBamInStream->pBamIterator);
      pBamInStream->pBamIterator = NULL;
    }

    // clear the bam array before jump
    SR_BamInStreamReset(pBamInStream);

    pBamInStream->pBamIterator = (bam_iter_t*) malloc(sizeof(bam_iter_t));

    // jump and read the first alignment in the given chromosome
    int ret;
    //bam_iter_t pBamIter = bam_iter_query(pBamInStream->pBamIndex, refID, begin, end);
    *pBamInStream->pBamIterator = bam_iter_query(pBamInStream->pBamIndex, refID, begin, end);

    pBamInStream->pNewNode = SR_BamNodeAlloc(pBamInStream->pMemPool);
    ret = bam_iter_read(pBamInStream->fpBamInput, *(pBamInStream->pBamIterator), &(pBamInStream->pNewNode->alignment));
    //bam_iter_destroy(pBamIter);

    khash_t(queryName)* pNameHashCurr = NULL;

    // see if we jump to the desired chromosome
    if (ret > 0 && pBamInStream->pNewNode->alignment.core.tid == refID)
    {
        // exclude those reads who are non-paired-end, qc-fail, duplicate-marked, proper-paired, 
        // both aligned, secondary-alignment and no-name-specified.
        if ((pBamInStream->pNewNode->alignment.core.flag & SR_BAM_FMASK) != 0
            || strcmp(bam1_qname(&(pBamInStream->pNewNode->alignment)), "*") == 0)
        {
            SR_BamNodeFree(pBamInStream->pNewNode, pBamInStream->pMemPool);
            pBamInStream->pNewNode = NULL;
            pBamInStream->currBinPos = NO_QUERY_YET;
        }
        else
        {
            SR_BamListPushHead(&(pBamInStream->pAlgnLists[CURR_BIN]), pBamInStream->pNewNode);

            int khRet = 0;
            khiter_t khIter = kh_put(queryName, pBamInStream->pNameHashes[CURR_BIN], bam1_qname(&(pBamInStream->pNewNode->alignment)), &khRet);

            if (khRet != 0)
            {
                pNameHashCurr = pBamInStream->pNameHashes[CURR_BIN];
                kh_value(pNameHashCurr, khIter) = pBamInStream->pNewNode;
            }
            else
                return SR_ERR;

            pBamInStream->currBinPos = pBamInStream->pNewNode->alignment.core.pos; 
            pBamInStream->pNewNode = NULL;
        }

        pBamInStream->currRefID = refID;
        return SR_OK;
    }
    else if (ret == -1)
    {
        return SR_OUT_OF_RANGE;
    }
    else
    {
        return SR_ERR;
    }
}

// read the header of a bam file
SR_BamHeader* SR_BamInStreamLoadHeader(SR_BamInStream* pBamInStream)
{
    bam_header_t* pOrigHeader = bam_header_read(pBamInStream->fpBamInput);
    if (pOrigHeader == NULL)
        return NULL;

    SR_BamHeader* pBamHeader = SR_BamHeaderAlloc();

    pBamHeader->pOrigHeader = pOrigHeader;

    pBamHeader->pMD5s = (const char**) calloc(pOrigHeader->n_targets, sizeof(char*));
    if (pBamHeader->pMD5s == NULL)
        SR_ErrQuit("ERROR: Not enough memory for md5 string");

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
            SR_ErrMsg("WARNING: Number of MD5 string is not consistent with number of chromosomes.");
    }

    return pBamHeader;
}

// load a pair from a bam file
SR_Status SR_BamInStreamLoadPair(SR_BamNode** ppUpAlgn, 
                                 SR_BamNode** ppDownAlgn, 
                                 SR_BamInStream* pBamInStream, 
				 bamFile* bam_writer_complete_bam) 
{
    khash_t(queryName)* pNameHashPrev = pBamInStream->pNameHashes[PREV_BIN];
    khash_t(queryName)* pNameHashCurr = pBamInStream->pNameHashes[CURR_BIN];

    int ret = 1;
    while(ret > 0 && (ret = SR_BamInStreamLoadNext(pBamInStream)) > 0)
    {
	// exclude those reads who are non-paired-end, qc-fail, duplicate-marked, proper-paired?!, 
        // both aligned, secondary-alignment and no-name-specified.
        SR_Bool shouldBeFiltered = pBamInStream->filterFunc(pBamInStream->pNewNode, pBamInStream->filterData);
        if (shouldBeFiltered)
        {
	    if (bam_writer_complete_bam != NULL) bam_write1(*bam_writer_complete_bam, &(pBamInStream->pNewNode->alignment));
	    
	    SR_BamNodeFree(pBamInStream->pNewNode, pBamInStream->pMemPool);
            pBamInStream->pNewNode = NULL;
            continue;
        }

        // update the current ref ID or position if the incoming alignment has a 
        // different value. The name hash and the bam array will be reset
        if (pNameHashPrev != NULL 
            && (pBamInStream->pNewNode->alignment.core.tid != pBamInStream->currRefID
                || pBamInStream->pNewNode->alignment.core.pos >= pBamInStream->currBinPos + 2 * pBamInStream->binLen))
        {
            if (pBamInStream->pNewNode->alignment.core.tid != pBamInStream->currRefID)
            {
                ret = SR_OUT_OF_RANGE; // different chromosome id
            }

            pBamInStream->currRefID  = pBamInStream->pNewNode->alignment.core.tid;
            pBamInStream->currBinPos = pBamInStream->pNewNode->alignment.core.pos;

            // Clear the hash buffer
	    kh_clear(queryName, pNameHashPrev);
            kh_clear(queryName, pNameHashCurr);

            // Store alignments before releasing them
            if (bam_writer_complete_bam != NULL) {
	      SR_BamNode* cur = pBamInStream->pAlgnLists[PREV_BIN].first;
	      for (int i = 0; i < pBamInStream->pAlgnLists[PREV_BIN].numNode; ++i) {
	        // if the cur is not NULL, store the cur in the complete bam
		if (cur != NULL) bam_write1(*bam_writer_complete_bam, &(cur->alignment));
		cur = cur->next;
	      } // end for

	      cur = pBamInStream->pAlgnLists[CURR_BIN].first;
	      for (int i = 0; i < pBamInStream->pAlgnLists[CURR_BIN].numNode; ++i) {
	        // if the cur is not NULL, store the cur in the complete bam
		if (cur != NULL) bam_write1(*bam_writer_complete_bam, &(cur->alignment));
		cur = cur->next;
	      } // end for
	    } // end if
	      
	    SR_BamListReset(&(pBamInStream->pAlgnLists[PREV_BIN]), pBamInStream->pMemPool);
            SR_BamListReset(&(pBamInStream->pAlgnLists[CURR_BIN]), pBamInStream->pMemPool);

        }
        else if (pBamInStream->pNewNode->alignment.core.pos >= pBamInStream->currBinPos + pBamInStream->binLen)
        {
            pBamInStream->currBinPos += pBamInStream->binLen;

            kh_clear(queryName, pNameHashPrev);
            SR_SWAP(pNameHashPrev, pNameHashCurr, khash_t(queryName)*);

            // Store alignments before releasing them
	    if (bam_writer_complete_bam != NULL) {
	      SR_BamNode* cur = pBamInStream->pAlgnLists[PREV_BIN].first;
	      for (int i = 0; i < pBamInStream->pAlgnLists[PREV_BIN].numNode; ++i) {
	        // if the cur is not NULL, store the cur in the complete bam
		if (cur != NULL) bam_write1(*bam_writer_complete_bam, &(cur->alignment));
		cur = cur->next;
              }
	    } // end if

	    SR_BamListReset(&(pBamInStream->pAlgnLists[PREV_BIN]), pBamInStream->pMemPool);

            SR_SWAP(pBamInStream->pAlgnLists[PREV_BIN], pBamInStream->pAlgnLists[CURR_BIN], SR_BamList);
        }
	else
	{
	} // end if-elseif-else

        SR_BamListPushHead(&(pBamInStream->pAlgnLists[CURR_BIN]), pBamInStream->pNewNode);

        // Clear the pointers
	(*ppUpAlgn) = NULL;
        (*ppDownAlgn) = NULL;
        
	khiter_t khIter = 0;

        // try to get the mate from the previous bin
	if (pNameHashPrev != NULL)
            khIter = kh_get(queryName, pNameHashPrev, bam1_qname(&(pBamInStream->pNewNode->alignment)));

        // the mate is found in the previous bin
	if (pNameHashPrev != NULL && khIter != kh_end(pNameHashPrev))
        {
            ret = SR_OK;
            (*ppUpAlgn) = kh_value(pNameHashPrev, khIter);
            (*ppDownAlgn) = pBamInStream->pNewNode;

            kh_del(queryName, pNameHashPrev, khIter);

            SR_BamListRemove(&(pBamInStream->pAlgnLists[PREV_BIN]), (*ppUpAlgn));
            SR_BamListRemove(&(pBamInStream->pAlgnLists[CURR_BIN]), (*ppDownAlgn));
        }
	// the mate is not in the previous bin
        else
        {
            int khRet = 0;
            khIter = kh_put(queryName, pNameHashCurr, bam1_qname(&(pBamInStream->pNewNode->alignment)), &khRet);

            if (khRet == 0) // we found a pair of alignments in the current bin
            {
                ret = SR_OK;
                (*ppUpAlgn) = kh_value(pNameHashCurr, khIter);
                (*ppDownAlgn) = pBamInStream->pNewNode;

                kh_del(queryName, pNameHashCurr, khIter);

                SR_BamListRemove(&(pBamInStream->pAlgnLists[CURR_BIN]), (*ppUpAlgn));
                SR_BamListRemove(&(pBamInStream->pAlgnLists[CURR_BIN]), (*ppDownAlgn));
            }
            else // not finding corresponding mate, save the current value and move on
            {
                kh_value(pNameHashCurr, khIter) = pBamInStream->pNewNode;
            }
        }
    } // end while

    pBamInStream->pNameHashes[PREV_BIN] = pNameHashPrev;
    pBamInStream->pNameHashes[CURR_BIN] = pNameHashCurr;

    if (ret < 0)
    {
        if ((ret == SR_EOF) && (bam_writer_complete_bam != NULL)) {
            // Store alignments before releasing them
	    SR_BamNode* cur = pBamInStream->pAlgnLists[PREV_BIN].first;
	    for (int i = 0; i < pBamInStream->pAlgnLists[PREV_BIN].numNode; ++i) {
	      // if the cur is not NULL, store the cur in the complete bam
	        if (cur != NULL) bam_write1(*bam_writer_complete_bam, &(cur->alignment));
		cur = cur->next;
	      } // end for

	    cur = pBamInStream->pAlgnLists[CURR_BIN].first;
	    for (int i = 0; i < pBamInStream->pAlgnLists[CURR_BIN].numNode; ++i) {
	      // if the cur is not NULL, store the cur in the complete bam
              if (cur != NULL) bam_write1(*bam_writer_complete_bam, &(cur->alignment));
	      cur = cur->next;
	    } // end for

	    SR_BamListReset(&(pBamInStream->pAlgnLists[PREV_BIN]), pBamInStream->pMemPool);
            SR_BamListReset(&(pBamInStream->pAlgnLists[CURR_BIN]), pBamInStream->pMemPool);
        }

	if ( ret != SR_OUT_OF_RANGE && ret != SR_EOF)
            return SR_ERR;
    }

    return ret;
}

unsigned int SR_BamInStreamShrinkPool(SR_BamInStream* pBamInStream, unsigned int newSize)
{
    unsigned int currSize = pBamInStream->pMemPool->numBuffs;
    if (currSize > newSize)
    {
        khash_t(buffAddress)* buffHash = kh_init(buffAddress);
        kh_resize(buffAddress, buffHash, currSize);

        int ret = 0;
        khiter_t khIter = 0;
        int64_t address = 0;

        for (SR_BamNode* pUsedNode = pBamInStream->pAlgnLists[PREV_BIN].first; pUsedNode != NULL; pUsedNode = pUsedNode->next)
        {
            address = (int64_t) pUsedNode->whereFrom;
            kh_put(buffAddress, buffHash, address, &ret);
        }

        for (SR_BamNode* pUsedNode = pBamInStream->pAlgnLists[CURR_BIN].first; pUsedNode != NULL; pUsedNode = pUsedNode->next)
        {
            address = (int64_t) pUsedNode->whereFrom;
            kh_put(buffAddress, buffHash, address, &ret);
        }

        unsigned int delNum = currSize - newSize;
        SR_BamBuff* pPrevBuff = NULL;
        SR_BamBuff* pCurrBuff = pBamInStream->pMemPool->pFirstBuff;

        while (pCurrBuff != NULL && delNum != 0)
        {
            int64_t address = (int64_t) pCurrBuff;
            khIter = kh_get(buffAddress, buffHash, address);

            if (khIter == kh_end(buffHash))
            {
                SR_BamBuff* pDelBuff = pCurrBuff;
                pCurrBuff = pCurrBuff->nextBuff;

                SR_BamBuffClear(pDelBuff, pBamInStream->pMemPool);
                SR_BamBuffFree(pDelBuff, pBamInStream->pMemPool->buffCapacity);
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

