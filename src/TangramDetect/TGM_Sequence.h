/*
 * =====================================================================================
 *
 *       Filename:  TGM_Sequence.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/11/2012 10:14:09 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef  TGM_SEQUENCE_H
#define  TGM_SEQUENCE_H

#include "TGM_Error.h"
#include <stdlib.h>
#include <string.h>
#include <string>
#include <math.h>

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

/* This table is used to transform nucleotide letters into numbers. */
const static char nt_table[128] = {
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
        4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
        4, 4, 4, 4,  3, 0, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4 
};

const static int8_t rc_table[5] = {3, 2, 1, 0, 4};

const static char char_table[5] = {'A', 'C', 'G', 'T', 'N'};

typedef struct TGM_Sequence
{
    int8_t* seq;

    size_t len;

    size_t cap;

}TGM_Sequence;

#define TGM_SeqInit(sequence) {(sequence).seq = NULL, (sequence).len = 0, (sequence).cap = 0}

static inline void TGM_SeqCpy(TGM_Sequence* dst, const char* src, size_t srcLen)
{
    if (srcLen + 1 >= dst->cap) 
    {
        dst->cap = srcLen + 2;
        kroundup32(dst->cap);
        dst->seq = (int8_t*) realloc(dst->seq, dst->cap * sizeof(int8_t));
        if (dst->seq == NULL)
            TGM_ErrQuit("ERROR: Not enough memory for the sequence.\n");
    }

    for (unsigned int i = 0; i != srcLen; ++i)
        dst->seq[i] = (int8_t) nt_table[(int8_t) src[i]];

    dst->len = srcLen;
}

static inline void TGM_SeqDup(TGM_Sequence* dst, const TGM_Sequence* src)
{
    if (src->len + 1 >= dst->cap) 
    {
        free(dst->seq);

        dst->cap = src->len + 2;
        kroundup32(dst->cap);
        dst->seq = (int8_t*) malloc(dst->cap * sizeof(int8_t));
        if (dst->seq == NULL)
            TGM_ErrQuit("ERROR: Not enough memory for the sequence.\n");
    }

    memcpy(dst->seq, src->seq, sizeof(int8_t) * src->len);
    dst->len = src->len;
}

static inline void TGM_SeqRevComp(TGM_Sequence* pRead)
{
    int8_t temp = 0;
    unsigned int start = 0;
    unsigned int end = pRead->len - 1;

    while (start < end)
    {
        temp = pRead->seq[start];
        pRead->seq[start] = rc_table[pRead->seq[end]];
        pRead->seq[end] = rc_table[temp];

        ++start;
        --end;
    }

    if (start == end) 
        pRead->seq[start] = rc_table[pRead->seq[start]];			
}

static inline void TGM_SeqClean(TGM_Sequence* pSeq)
{
    if (pSeq != NULL)
        free(pSeq->seq);
}

static inline void SeqToString(std::string& seqStr, const int8_t* pSeq, unsigned int seqLen)
{
    seqStr.clear();
    for (unsigned int i = 0; i != seqLen; ++i)
    {
        seqStr += char_table[pSeq[i]];
    }
}

static inline void ReverseSeqToString(std::string& seqStr, const int8_t* pSeq, unsigned int seqLen)
{
    seqStr.clear();
    int8_t *ptr = pSeq + seqLen - 1;
    for (int i = seqLen; i > 0; --i)
    {
        seqStr += char_table[rc_table[ptr[i]]];
	--ptr;
    }
}

static inline double SeqGetEntropy(const int8_t* readSeq, int readLen)
{
    int count[4] = {0, 0, 0, 0};
    double entropy = 0.0;
    int len = readLen;
    for (int i = 0; i != readLen; ++i)
    {
        if (readSeq[i] < 4)
            ++count[readSeq[i]];
        else
            --len;
    }

    if (len == 0)
        return entropy;

    for (int i = 0; i != 4; ++i)
    {
        if (count[i] > 0)
        {
            double p = (double) count[i] / len;
            entropy += p * log2(p);
        }
    }

    return (-entropy);
}

#endif  /*TGM_SEQUENCE_H*/
