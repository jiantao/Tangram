/*
 * =====================================================================================
 *
 *       Filename:  TGM_BamMemPool.h
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

#ifndef  TGM_BAMMEMPOOL_H
#define  TGM_BAMMEMPOOL_H

#include "bam.h"
#include "TGM_Types.h"

typedef struct TGM_BamNode TGM_BamNode;

typedef TGM_BamNode* TGM_BamListIter;

typedef struct TGM_BamList TGM_BamList;

typedef struct TGM_BamBuff TGM_BamBuff;

typedef struct TGM_BamMemPool TGM_BamMemPool;

struct TGM_BamNode
{
    bam1_t alignment;

    TGM_BamNode* prev;

    TGM_BamNode* next;

    const TGM_BamBuff* whereFrom;
};

struct TGM_BamList
{
    unsigned int numNode;

    TGM_BamNode* first;

    TGM_BamNode* last;
};

struct TGM_BamBuff
{
    unsigned int numUsed;

    TGM_BamBuff* nextBuff;

    TGM_BamNode* pNodeArray;
};

struct TGM_BamMemPool
{
    unsigned int numBuffs;

    unsigned int buffCapacity;

    TGM_BamBuff* pFirstBuff;

    TGM_BamList avlNodeList;
};

TGM_BamNode* TGM_BamNodeAlloc(TGM_BamMemPool* pMemPool);

#define TGM_BamNodeFree(pNode, pMemPool)  TGM_BamListPushHead(&((pMemPool)->avlNodeList), (pNode))

static inline TGM_BamListIter TGM_BamListGetIter(TGM_BamList* pList)
{
    return pList->first;
}

void TGM_BamListPushHead(TGM_BamList* pList, TGM_BamNode* pNewFirstNode);

void TGM_BamListPushBack(TGM_BamList* pList, TGM_BamNode* pNewLastNode);

TGM_BamNode* TGM_BamListPopHead(TGM_BamList* pList);

void TGM_BamListRemove(TGM_BamList* pList, TGM_BamNode* pNode);

void TGM_BamListMergeHead(TGM_BamList* dstList, TGM_BamNode* first, TGM_BamNode* last, unsigned int numNode);

void TGM_BamListReset(TGM_BamList* pList, TGM_BamMemPool* pMemPool);


TGM_BamBuff* TGM_BamBuffAlloc(unsigned int buffCapacity);

void TGM_BamBuffFree(TGM_BamBuff* pBuff, unsigned int buffCapacity);

void TGM_BamBuffClear(TGM_BamBuff* pBuff, TGM_BamMemPool* pMemPool);


TGM_BamMemPool* TGM_BamMemPoolAlloc(unsigned int buffCapacity);

void TGM_BamMemPoolFree(TGM_BamMemPool* pMemPool);

#define TGM_BamMemPoolGetSize(pMemPool) ((pMemPool)->numBuffs)

TGM_Status TGM_BamMemPoolExpand(TGM_BamMemPool* pMemPool);


#endif  /*TGM_BAMMEMPOOL_H*/
