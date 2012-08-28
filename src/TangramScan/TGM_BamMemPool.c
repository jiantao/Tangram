/*
 * =====================================================================================
 *
 *       Filename:  TGM_BamMemPool.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  09/13/2011 04:21:38 PM
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

#include "TGM_Error.h"
#include "TGM_BamMemPool.h"

#define TGM_MAX_MEM_POOL_SIZE 2000000

static inline TGM_Bool TGM_BamListIsEmpty(TGM_BamList* pList)
{
    return (pList->numNode == 0);
}

TGM_BamNode* TGM_BamNodeAlloc(TGM_BamMemPool* pMemPool)
{
    if (TGM_BamListIsEmpty(&(pMemPool->avlNodeList)))
    {
        TGM_Status memStatus = TGM_BamMemPoolExpand(pMemPool);
        if (memStatus == TGM_OVER_FLOW)
            return NULL;
    }

    TGM_BamNode* pNewNode = TGM_BamListPopHead(&(pMemPool->avlNodeList));

    pNewNode->prev = NULL;
    pNewNode->next = NULL;

    return pNewNode;
}

void TGM_BamListPushHead(TGM_BamList* pList, TGM_BamNode* pNewFirstNode)
{
    pNewFirstNode->next = pList->first;
    pNewFirstNode->prev = NULL;

    if (pList->first != NULL)
        pList->first->prev = pNewFirstNode;
    else
        pList->last = pNewFirstNode;

    pList->first = pNewFirstNode;
    ++(pList->numNode);
}

void TGM_BamListPushBack(TGM_BamList* pList, TGM_BamNode* pNewLastNode)
{
    pNewLastNode->prev = pList->last;
    pNewLastNode->next = NULL;

    if (pList->last != NULL)
        pList->last->next = pNewLastNode;
    else
        pList->first = pNewLastNode;

    pList->last = pNewLastNode;
    ++(pList->numNode);
}

TGM_BamNode* TGM_BamListPopHead(TGM_BamList* pList)
{
    TGM_BamNode* pPopNode = pList->first;

    if (pList->first != NULL)
    {
        pList->first = pList->first->next;
        if (pList->first != NULL)
        {
            pList->first->prev = NULL;
        }
        else
        {
            pList->last = NULL;
        }

        --(pList->numNode);
    }

    return pPopNode;
}

void TGM_BamListRemove(TGM_BamList* pList, TGM_BamNode* pNode)
{
    assert(pList->numNode != 0);

    if (pNode->next != NULL)
        pNode->next->prev = pNode->prev;
    else
        pList->last = pNode->prev;

    if (pNode->prev != NULL)
        pNode->prev->next = pNode->next;
    else
        pList->first = pNode->next;

    pNode->next = NULL;
    pNode->prev = NULL;

    --(pList->numNode);
}

void TGM_BamListMergeHead(TGM_BamList* dstList, TGM_BamNode* first, TGM_BamNode* last, unsigned int numNode)
{
    if (numNode == 0)
        return;

    last->next = dstList->first;
    if (dstList->first != NULL)
    {
        dstList->first->prev = last;
    }
    else
    {
        dstList->last = last;
    }

    first->prev = NULL;
    dstList->first = first;

    dstList->numNode += numNode;
}

void TGM_BamListReset(TGM_BamList* pList, TGM_BamMemPool* pMemPool)
{
    TGM_BamListMergeHead(&(pMemPool->avlNodeList), pList->first, pList->last, pList->numNode);

    pList->first = NULL;
    pList->last = NULL;
    pList->numNode = 0;
}

TGM_BamBuff* TGM_BamBuffAlloc(unsigned int buffCapacity)
{
    TGM_BamBuff* pNewBuff = (TGM_BamBuff*) malloc(sizeof(TGM_BamBuff));
    if (pNewBuff == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for the bam memory node object.\n");

    pNewBuff->numUsed = 0;
    pNewBuff->nextBuff = NULL;

    pNewBuff->pNodeArray = (TGM_BamNode*) calloc(buffCapacity, sizeof(TGM_BamNode));
    if (pNewBuff->pNodeArray == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for the storage of alignments in the bam memory node object.\n");

    
    for (unsigned int i = 0; i != buffCapacity - 1; ++i)
    {
        pNewBuff->pNodeArray[i].whereFrom = pNewBuff;
        pNewBuff->pNodeArray[i].next = &(pNewBuff->pNodeArray[i + 1]);
        pNewBuff->pNodeArray[i + 1].prev = &(pNewBuff->pNodeArray[i]);
    }

    pNewBuff->pNodeArray[0].prev = NULL;
    pNewBuff->pNodeArray[buffCapacity - 1].next = NULL;
    pNewBuff->pNodeArray[buffCapacity - 1].whereFrom = pNewBuff;

    return pNewBuff;
}

void TGM_BamBuffFree(TGM_BamBuff* pBuff, unsigned int buffCapacity)
{
    if (pBuff != NULL)
    {
        for (unsigned int i = 0; i != buffCapacity; ++i)
            free(pBuff->pNodeArray[i].alignment.data);

        free(pBuff->pNodeArray);
        free(pBuff);
    }
}

void TGM_BamBuffClear(TGM_BamBuff* pBuff, TGM_BamMemPool* pMemPool)
{
    for (unsigned int i = 0; i != pMemPool->buffCapacity; ++i)
    {
        TGM_BamNode* pNode = pBuff->pNodeArray + i;
        TGM_BamListRemove(&(pMemPool->avlNodeList), pNode);
    }
}

TGM_Status TGM_BamMemPoolExpand(TGM_BamMemPool* pMemPool)
{
    if ((pMemPool->numBuffs + 1) * pMemPool->buffCapacity > TGM_MAX_MEM_POOL_SIZE)
        return TGM_OVER_FLOW;

    TGM_BamBuff* pNewBuff = TGM_BamBuffAlloc(pMemPool->buffCapacity);

    pNewBuff->nextBuff = pMemPool->pFirstBuff;
    pMemPool->pFirstBuff = pNewBuff;

    TGM_BamListMergeHead(&(pMemPool->avlNodeList), &(pNewBuff->pNodeArray[0]), &(pNewBuff->pNodeArray[pMemPool->buffCapacity - 1]), pMemPool->buffCapacity);
    ++(pMemPool->numBuffs);

    return TGM_OK;
}

TGM_BamMemPool* TGM_BamMemPoolAlloc(unsigned int buffCapacity)
{
    TGM_BamMemPool* pNewPool = (TGM_BamMemPool*) calloc(1, sizeof(TGM_BamMemPool));
    if (pNewPool == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for the bam memory pool object.\n");

    pNewPool->numBuffs = 0;
    pNewPool->pFirstBuff = NULL;
    pNewPool->buffCapacity = buffCapacity;

    TGM_BamMemPoolExpand(pNewPool);

    return pNewPool;
}

void TGM_BamMemPoolFree(TGM_BamMemPool* pMemPool)
{
    if (pMemPool != NULL)
    {
        TGM_BamBuff* pDelBuff = pMemPool->pFirstBuff;
        TGM_BamBuff* pNextBuff = pMemPool->pFirstBuff;
        while (pDelBuff != NULL)
        {
            pNextBuff = pDelBuff->nextBuff;
            TGM_BamBuffFree(pDelBuff, pMemPool->buffCapacity);
            pDelBuff = pNextBuff;
        }

        free(pMemPool);
    }
}
