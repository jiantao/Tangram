/*
 * =====================================================================================
 *
 *       Filename:  SR_BamMemPool.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  09/13/2011 04:14:34 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef  SR_BAMMEMPOOL_H
#define  SR_BAMMEMPOOL_H

#include "../OutSources/samtools/bam.h"
#include "SR_Types.h"

typedef struct SR_BamNode SR_BamNode;

typedef SR_BamNode* SR_BamListIter;

typedef struct SR_BamList SR_BamList;

typedef struct SR_BamBuff SR_BamBuff;

typedef struct SR_BamMemPool SR_BamMemPool;

struct SR_BamNode
{
    bam1_t alignment;

    SR_BamNode* prev;

    SR_BamNode* next;

    const SR_BamBuff* whereFrom;
};

struct SR_BamList
{
    unsigned int numNode;

    SR_BamNode* first;

    SR_BamNode* last;
};

struct SR_BamBuff
{
    unsigned int numUsed;

    SR_BamBuff* nextBuff;

    SR_BamNode* pNodeArray;
};

struct SR_BamMemPool
{
    unsigned short numBuffs;

    unsigned int buffCapacity;

    SR_BamBuff* pFirstBuff;

    SR_BamList avlNodeList;
};

SR_BamNode* SR_BamNodeAlloc(SR_BamMemPool* pMemPool);

#define SR_BamNodeFree(pNode, pMemPool)  SR_BamListPushHead(&((pMemPool)->avlNodeList), (pNode))

static inline SR_BamListIter SR_BamListGetIter(SR_BamList* pList)
{
    return pList->first;
}

void SR_BamListPushHead(SR_BamList* pList, SR_BamNode* pNewFirstNode);

void SR_BamListPushBack(SR_BamList* pList, SR_BamNode* pNewLastNode);

SR_BamNode* SR_BamListPopHead(SR_BamList* pList);

void SR_BamListRemove(SR_BamList* pList, SR_BamNode* pNode);

void SR_BamListMergeHead(SR_BamList* dstList, SR_BamNode* first, SR_BamNode* last, unsigned int numNode);

void SR_BamListReset(SR_BamList* pList, SR_BamMemPool* pMemPool);


SR_BamBuff* SR_BamBuffAlloc(unsigned int buffCapacity);

void SR_BamBuffFree(SR_BamBuff* pBuff, unsigned int buffCapacity);

void SR_BamBuffClear(SR_BamBuff* pBuff, SR_BamMemPool* pMemPool);


SR_BamMemPool* SR_BamMemPoolAlloc(unsigned int buffCapacity);

void SR_BamMemPoolFree(SR_BamMemPool* pMemPool);

#define SR_BamMemPoolGetSize(pMemPool) ((pMemPool)->numBuffs)

SR_Status SR_BamMemPoolExpand(SR_BamMemPool* pMemPool);


#endif  /*SR_BAMMEMPOOL_H*/
