/*
 * =====================================================================================
 *
 *       Filename:  HashRegionTable.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/03/2011 15:23:54
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "SR_Error.h"
#include "SR_Utilities.h"
#include "SR_HashRegionTable.h"

// default capacity of a hash region array
const int DEFAULT_HASH_ARR_CAPACITY = 50;

//=========================
// Static methods
//=========================

// intialize the best hash region array
static void ResetBestRegions(HashRegionTable* pRegionTable, unsigned short queryLen)
{

    if (pRegionTable->pBestCloseRegions == NULL)
    {
        SR_ARRAY_ALLOC(pRegionTable->pBestCloseRegions, queryLen, BestRegionArray, BestRegion);
        SR_ARRAY_ALLOC(pRegionTable->pBestFarRegions, queryLen, BestRegionArray, BestRegion);

        for (unsigned int i = 0; i != queryLen; ++i)
        {
            BestRegion* tempRegion = pRegionTable->pBestCloseRegions->data + i;
            tempRegion->queryBegin = i;

            tempRegion = pRegionTable->pBestFarRegions->data + i;
            tempRegion->queryBegin = i;
        }
    }

    SR_Bool isSmall = FALSE;
    if ((pRegionTable->pBestCloseRegions)->capacity < queryLen)
    {
        isSmall = TRUE;
        SR_ARRAY_FREE(pRegionTable->pBestCloseRegions, FALSE);
        SR_ARRAY_FREE(pRegionTable->pBestFarRegions, FALSE);
    }

    if ((pRegionTable->pBestCloseRegions)->data == NULL)
    {
        unsigned int capacity = isSmall ? 2 * queryLen : queryLen;

        SR_ARRAY_INIT(pRegionTable->pBestCloseRegions, capacity, BestRegion);
        SR_ARRAY_INIT(pRegionTable->pBestFarRegions, capacity, BestRegion);

        for (unsigned int i = 0; i != capacity; ++i)
        {
            BestRegion* tempRegion = pRegionTable->pBestCloseRegions->data + i;
            tempRegion->queryBegin = i;

            tempRegion = pRegionTable->pBestFarRegions->data + i;
            tempRegion->queryBegin = i;
        }
    }

    (pRegionTable->pBestCloseRegions)->size = queryLen;
    (pRegionTable->pBestFarRegions)->size = queryLen;

    for (unsigned int i = 0; i != queryLen; ++i)
    {
        BestRegion* tempRegion = pRegionTable->pBestCloseRegions->data + i;
        tempRegion->length     = 0;
        tempRegion->numPos     = 0;
	//tempRegion->queryBegin = i;
	//for (unsigned int j = 0; j < MAX_BEST_REF_BEGINS; ++j)
	//  tempRegion->refBegins[j] = 0;

        tempRegion = pRegionTable->pBestFarRegions->data + i;
        tempRegion->length     = 0;
        tempRegion->numPos     = 0;
	//tempRegion->queryBegin = i;
	//for (unsigned int j = 0; j < MAX_BEST_REF_BEGINS; ++j)
	//  tempRegion->refBegins[j] = 0;
    }
}


// calculate the hash key for a hash in a read that starts at a certain position
static SR_Bool GetNextHashKey(const char* query, uint32_t queryLen, unsigned int* pPos, uint32_t* pHashKey, uint32_t mask, unsigned int hashSize)
{
    // table use to translate a nucleotide into its corresponding 2-bit representation
    static const char translation[26] = { 0, -1, 1, -1, -1, -1, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 3, -1, -1, -1, -1, -1, -1 };

    unsigned int endPos = *pPos + hashSize;
    unsigned int startPos = *pPos == 0 ? *pPos : endPos - 1;

    // remove the highest 2 bits in the previous hash key
    *pHashKey &= mask;

    while (endPos <= queryLen)
    {
        char tValue = 0;
        for (unsigned int i = startPos; i != endPos; ++i)
        {
            tValue = translation[query[i] - 'A'];
            if (tValue < 0)
            {
                *pHashKey = 0;
                *pPos = i + 1;
                startPos = *pPos;
                endPos = startPos + hashSize;
                break;
            }

            *pHashKey = (*pHashKey) << 2 | tValue;
        }

        if (tValue >= 0)
            return TRUE;
    }

    return FALSE;
}

// find the start index of the first hash position that is in our search region
static unsigned int GetStartHashPosIndex(const HashPosView* pHashPosView, uint32_t refStart)
{
    unsigned int min = 0;
    unsigned int max = pHashPosView->size - 1;

    unsigned int mid;
    while (min <= max)
    {
        mid = (min + max) / 2;
        
        if (refStart == pHashPosView->data[mid])
            return mid;
        else if (refStart > pHashPosView->data[mid])
            min = mid + 1;
        else
        {
            if (mid == 0)
                break;

            max = mid - 1;
        }
    }

    if (pHashPosView->data[mid] > refStart)
        return mid;
    else
        return mid + 1;
}


// merge a new hash region with existing ones
static SR_Bool MergeHashRegions(HashRegionTable* pRegionTable, HashRegion* pNewRegion)
{
    if (pNewRegion->refBegin < SR_ARRAY_GET(pRegionTable->pPrevRegions, 0).refBegin)
        return TRUE;

    uint32_t target = pNewRegion->refBegin + pNewRegion->length - 1;
    unsigned int min = pRegionTable->searchBegin;
    unsigned int max = (pRegionTable->pPrevRegions)->size - 1;

    if (min > max)
        return FALSE;

    while (min <= max)
    {
        unsigned int mid = (min + max) / 2;
        uint32_t refEnd = SR_ARRAY_GET(pRegionTable->pPrevRegions, mid).refBegin + SR_ARRAY_GET(pRegionTable->pPrevRegions, mid).length;

        if (target == refEnd)
        {
            pRegionTable->searchBegin = mid + 1;
            pNewRegion->refBegin = SR_ARRAY_GET(pRegionTable->pPrevRegions, mid).refBegin;
            pNewRegion->queryBegin = SR_ARRAY_GET(pRegionTable->pPrevRegions, mid).queryBegin;
            pNewRegion->length = SR_ARRAY_GET(pRegionTable->pPrevRegions, mid).length + 1;
            break;
        }
        else if (target > refEnd)
        {
            min = mid + 1;
            pRegionTable->searchBegin = min;
        }
        else
        {
            if (mid == 0)
                return TRUE;

            max = mid - 1;
        }
    }

    return TRUE;
}

// update the best hash regions after each merge
static void UpdateBestRegions(HashRegionTable* pRegionTable, const HashRegion* pNewRegion, const SR_QueryRegion* pQueryRegion)
{
    BestRegion* pBestClose = SR_ARRAY_GET_PT(pRegionTable->pBestCloseRegions, pNewRegion->queryBegin);
    BestRegion* pBestFar = SR_ARRAY_GET_PT(pRegionTable->pBestFarRegions, pNewRegion->queryBegin);

    if (pNewRegion->length > pBestFar->length)
    {
        pBestFar->length = pNewRegion->length;
        pBestFar->refBegins[0] = pNewRegion->refBegin;
        pBestFar->numPos = 1;
    }
    else if (pNewRegion->length == pBestFar->length)
    {
        if (pBestFar->numPos < MAX_BEST_REF_BEGINS)
            pBestFar->refBegins[pBestFar->numPos] = pNewRegion->refBegin;

        ++(pBestFar->numPos);
    }

    if (pNewRegion->refBegin >= pQueryRegion->closeRefBegin 
        && pNewRegion->refBegin <= pQueryRegion->closeRefEnd)
    {
        if (pNewRegion->length > pBestClose->length)
        {
            pBestClose->length = pNewRegion->length;
            pBestClose->refBegins[0] = pNewRegion->refBegin;
            pBestClose->numPos = 1;
        }
        else if (pNewRegion->length == pBestClose->length)
        {
            if (pBestClose->numPos < MAX_BEST_REF_BEGINS)
                pBestClose->refBegins[pBestClose->numPos] = pNewRegion->refBegin;

            ++(pBestClose->numPos);
        }
    }
}


//===============================
// Constructors and Destructors
//===============================

HashRegionTable* HashRegionTableAlloc(void)
{
    HashRegionTable* pNewTable = (HashRegionTable*) malloc(sizeof(HashRegionTable));
    if (pNewTable == NULL)
        SR_ErrSys("ERROR: not enough memory for a new hash region table object.\n");

    SR_ARRAY_ALLOC(pNewTable->pPrevRegions, DEFAULT_HASH_ARR_CAPACITY, HashRegionArray, HashRegion);
    SR_ARRAY_ALLOC(pNewTable->pCurrRegions, DEFAULT_HASH_ARR_CAPACITY, HashRegionArray, HashRegion);

    pNewTable->pBestCloseRegions = NULL;
    pNewTable->pBestFarRegions = NULL;

    pNewTable->searchBegin = 0;

    return pNewTable;
}

void HashRegionTableFree(HashRegionTable* pRegionTable)
{
    if (pRegionTable != NULL)
    {
        SR_ARRAY_FREE(pRegionTable->pPrevRegions, TRUE);
        SR_ARRAY_FREE(pRegionTable->pCurrRegions, TRUE);
        SR_ARRAY_FREE(pRegionTable->pBestCloseRegions, TRUE);
        SR_ARRAY_FREE(pRegionTable->pBestFarRegions, TRUE);

        free(pRegionTable);
    }
}


//===============================
// Interface functions
//===============================

// for each query find the best hash regions in the reference
void HashRegionTableLoad(HashRegionTable* pRegionTable, const SR_InHashTable* pHashTable, const SR_QueryRegion* pQueryRegion)
{
    unsigned int prevQueryPos = 0;
    unsigned int currQueryPos = 0;
    uint32_t hashKey = 0;

    // get the next hash key in the query
    while (GetNextHashKey(pQueryRegion->orphanSeq, SR_GetQueryLen(pQueryRegion->pOrphan), &currQueryPos, &hashKey, pHashTable->highEndMask, pHashTable->hashSize))
    {
        // an array stores the hash positions under current hash key
        HashPosView hashPosArray;
        // this struct will store the new hash region from the query
        HashRegion newRegion;

        if (SR_InHashTableSearch(&hashPosArray, pHashTable, hashKey))
        {

            // we only have to merge the hash regions when we do get some hash regions in the last round
            // and the current query position is 1bp ahead the previous query position
            SR_Bool doMerge = FALSE;
            if (currQueryPos == prevQueryPos + 1 && SR_ARRAY_GET_SIZE(pRegionTable->pPrevRegions) > 0)
                doMerge = TRUE;

            // find the start index of the first hash position that is in our search region
            unsigned int startIndex = GetStartHashPosIndex(&hashPosArray, pQueryRegion->farRefBegin);
            for (unsigned int i = startIndex; i != hashPosArray.size; ++i)
            {
                // initialize the new hash region
                newRegion.queryBegin = currQueryPos;
                newRegion.refBegin = hashPosArray.data[i];
                newRegion.length = pHashTable->hashSize;

                if (doMerge)
                    doMerge = MergeHashRegions(pRegionTable, &newRegion);

                // we will get out of the loop if the hash position in reference exceeds our search region
                if (newRegion.refBegin > pQueryRegion->farRefEnd)
                    break;

                // we will update the best hash region with this new region and push it into the current hash region array for next round merge
                UpdateBestRegions(pRegionTable, &newRegion, pQueryRegion);
                SR_ARRAY_PUSH(pRegionTable->pCurrRegions, &newRegion, HashRegion);
            }
        }

        // we are done with the previos hash region array
        // we will clear it and swap it with the current hash region array
        SR_ARRAY_RESET(pRegionTable->pPrevRegions);
        SR_SWAP(pRegionTable->pPrevRegions, pRegionTable->pCurrRegions, HashRegionArray*);
        pRegionTable->searchBegin = 0;

        // move to the next hash position in query
        prevQueryPos = currQueryPos;
        ++currQueryPos;
    }
}

// index the best hash regions with their end position
void HashRegionTableReverseBest(HashRegionTable* pRegionTable)
{
    for (int i = SR_ARRAY_GET_SIZE(pRegionTable->pBestCloseRegions) - 1; i >= 0; --i)
    {
        BestRegion* pCloseLowEnd = SR_ARRAY_GET_PT(pRegionTable->pBestCloseRegions, i);
        if (pCloseLowEnd->length > 0)
        {
            unsigned short closeHighEndPos = i + pCloseLowEnd->length - 1;
            BestRegion* pCloseHighEnd = SR_ARRAY_GET_PT(pRegionTable->pBestCloseRegions, closeHighEndPos);

            if (pCloseHighEnd->length < pCloseLowEnd->length)
                *pCloseHighEnd = *pCloseLowEnd;

            pCloseLowEnd->length = 0;
        }

        BestRegion* pFarLowEnd = SR_ARRAY_GET_PT(pRegionTable->pBestFarRegions, i);
        if (pFarLowEnd->length > 0)
        {
            unsigned short closeHighEndPos = i + pFarLowEnd->length - 1;
            BestRegion* pFarHighEnd = SR_ARRAY_GET_PT(pRegionTable->pBestFarRegions, closeHighEndPos);

            if (pFarHighEnd->length < pFarLowEnd->length)
                *pFarHighEnd = *pFarLowEnd;

            pFarLowEnd->length = 0;
        }
    }
}


// initialize the hash region table for a new query
void HashRegionTableInit(HashRegionTable* pRegionTable, uint32_t queryLen)
{
    SR_ARRAY_RESET(pRegionTable->pPrevRegions);
    SR_ARRAY_RESET(pRegionTable->pCurrRegions);

    ResetBestRegions(pRegionTable, queryLen);
}
