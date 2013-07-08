/*
 * =====================================================================================
 *
 *       Filename:  SR_BamMemPool.c
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

#include "SR_BamMemPool.h"

#include <stdlib.h>
#include <assert.h>

#include "SR_Error.h"

#define SR_MAX_MEM_POOL_SIZE 2000000

static inline SR_Bool SR_BamListIsEmpty(SR_BamList* pList)
{
    return (pList->numNode == 0);
}

SR_BamNode* SR_BamNodeAlloc(SR_BamMemPool* pMemPool)
{
    if (SR_BamListIsEmpty(&(pMemPool->avlNodeList)))
    {
        SR_Status memStatus = SR_BamMemPoolExpand(pMemPool);
        if (memStatus == SR_OVER_FLOW)
            return NULL;
    }

    SR_BamNode* pNewNode = SR_BamListPopHead(&(pMemPool->avlNodeList));

    pNewNode->prev = NULL;
    pNewNode->next = NULL;

    return pNewNode;
}

void SR_BamListPushHead(SR_BamList* pList, SR_BamNode* pNewFirstNode)
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

void SR_BamListPushBack(SR_BamList* pList, SR_BamNode* pNewLastNode)
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

SR_BamNode* SR_BamListPopHead(SR_BamList* pList)
{
    SR_BamNode* pPopNode = pList->first;

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

void SR_BamListRemove(SR_BamList* pList, SR_BamNode* pNode)
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

void SR_BamListMergeHead(SR_BamList* dstList, SR_BamNode* first, SR_BamNode* last, unsigned int numNode)
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

void SR_BamListReset(SR_BamList* pList, SR_BamMemPool* pMemPool)
{
    SR_BamListMergeHead(&(pMemPool->avlNodeList), pList->first, pList->last, pList->numNode);

    pList->first = NULL;
    pList->last = NULL;
    pList->numNode = 0;
}

SR_BamBuff* SR_BamBuffAlloc(unsigned int buffCapacity)
{
    SR_BamBuff* pNewBuff = (SR_BamBuff*) malloc(sizeof(SR_BamBuff));
    if (pNewBuff == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the bam memory node object.\n");

    pNewBuff->numUsed = 0;
    pNewBuff->nextBuff = NULL;

    pNewBuff->pNodeArray = (SR_BamNode*) calloc(buffCapacity, sizeof(SR_BamNode));
    if (pNewBuff->pNodeArray == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the storage of alignments in the bam memory node object.\n");

    
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

void SR_BamBuffFree(SR_BamBuff* pBuff, unsigned int buffCapacity)
{
    if (pBuff != NULL)
    {
        for (unsigned int i = 0; i != buffCapacity; ++i)
            free(pBuff->pNodeArray[i].alignment.data);

        free(pBuff->pNodeArray);
        free(pBuff);
    }
}

void SR_BamBuffClear(SR_BamBuff* pBuff, SR_BamMemPool* pMemPool)
{
    for (unsigned int i = 0; i != pMemPool->buffCapacity; ++i)
    {
        SR_BamNode* pNode = pBuff->pNodeArray + i;
        SR_BamListRemove(&(pMemPool->avlNodeList), pNode);
    }
}

SR_Status SR_BamMemPoolExpand(SR_BamMemPool* pMemPool)
{
    if ((pMemPool->numBuffs + 1) * pMemPool->buffCapacity > SR_MAX_MEM_POOL_SIZE)
        return SR_OVER_FLOW;

    SR_BamBuff* pNewBuff = SR_BamBuffAlloc(pMemPool->buffCapacity);

    pNewBuff->nextBuff = pMemPool->pFirstBuff;
    pMemPool->pFirstBuff = pNewBuff;

    SR_BamListMergeHead(&(pMemPool->avlNodeList), &(pNewBuff->pNodeArray[0]), &(pNewBuff->pNodeArray[pMemPool->buffCapacity - 1]), pMemPool->buffCapacity);
    ++(pMemPool->numBuffs);

    return SR_OK;
}

SR_BamMemPool* SR_BamMemPoolAlloc(unsigned int buffCapacity)
{
    SR_BamMemPool* pNewPool = (SR_BamMemPool*) calloc(1, sizeof(SR_BamMemPool));
    if (pNewPool == NULL)
        SR_ErrQuit("ERROR: Not enough memory for the bam memory pool object.\n");

    pNewPool->numBuffs = 0;
    pNewPool->pFirstBuff = NULL;
    pNewPool->buffCapacity = buffCapacity;

    SR_BamMemPoolExpand(pNewPool);

    return pNewPool;
}

void SR_BamMemPoolFree(SR_BamMemPool* pMemPool)
{
    if (pMemPool != NULL)
    {
        SR_BamBuff* pDelBuff = pMemPool->pFirstBuff;
        SR_BamBuff* pNextBuff = pMemPool->pFirstBuff;
        while (pDelBuff != NULL)
        {
            pNextBuff = pDelBuff->nextBuff;
            SR_BamBuffFree(pDelBuff, pMemPool->buffCapacity);
            pDelBuff = pNextBuff;
        }

        free(pMemPool);
    }
}
