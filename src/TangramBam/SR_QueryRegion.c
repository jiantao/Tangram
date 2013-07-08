/*
 * =====================================================================================
 *
 *       Filename:  SR_QueryRegion.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  08/06/2011 12:42:45 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <limits.h>

#include "SR_Error.h"
#include "SR_Utilities.h"
#include "SR_QueryRegion.h"

#include "seq_converter.h"


//===============================
// Type and constant definition
//===============================

// map used to transfer the 4-bit representation of a nucleotide into the ascii representation
static const char SR_BASE_MAP[16]   = {'N','A','C','N','G','N','N','N','T','N','N','N','N','N','N','N'};
static const char SR_BASE_R_MAP[16] = {'N','T','G','N','C','N','N','N','A','N','N','N','N','N','N','N'};
//static const uint8_t SR_R_MAP[16]   = { 15, 8 , 4 , 15, 2 , 15, 15, 15, 1 , 15, 15, 15, 15, 15, 15, 15};


//===============================
// Constructors and Destructors
//===============================

SR_QueryRegion* SR_QueryRegionAlloc(void)
{
    SR_QueryRegion* pNewRegion = (SR_QueryRegion*) malloc(sizeof(SR_QueryRegion));
    if (pNewRegion == NULL)
        SR_ErrQuit("ERROR: Not enough memory for a orphanSeq region object");

    // be caureful, the actual data are not stored in these two pointers
    // they are stored in the bam aux data structure
    pNewRegion->pAnchor = NULL;
    pNewRegion->pOrphan = NULL;

    pNewRegion->orphanSeq                  = NULL;
    pNewRegion->orphanSeqForward           = NULL;
    pNewRegion->orphanSeqReverseComplament = NULL;
    pNewRegion->orphanSeqReverse           = NULL;
    pNewRegion->orphanSeqComplement        = NULL;
    pNewRegion->isOrphanInversed = FALSE;
    pNewRegion->capacity = 0;

    pNewRegion->closeRefBegin = 0;
    pNewRegion->closeRefEnd = 0;
    pNewRegion->farRefBegin = 0;
    pNewRegion->farRefEnd = 0;

    return pNewRegion;
}

void SR_QueryRegionFree(SR_QueryRegion* pQueryRegion)
{
    if (pQueryRegion != NULL)
    {
        free(pQueryRegion->orphanSeqForward);
	free(pQueryRegion->orphanSeqReverseComplament);
	free(pQueryRegion->orphanSeqReverse);
	free(pQueryRegion->orphanSeqComplement);
        free(pQueryRegion);
	pQueryRegion = NULL;
    }
}


//======================
// Interface functions
//======================

SR_Status SR_QueryRegionLoadPair(SR_QueryRegion* pQueryRegion, SR_BamInStreamIter* pIter)
{
    if (pIter->pBamNode == NULL)
        return SR_OUT_OF_RANGE;

    pQueryRegion->pAnchor = &(pIter->pBamNode->alignment);
    pIter->pBamNode = pIter->pBamNode->next;
    if (pIter->pBamNode == NULL)
        return SR_ERR;

    pQueryRegion->pOrphan = &(pIter->pBamNode->alignment);
    pIter->pBamNode = pIter->pBamNode->next;

    pQueryRegion->pOrphan->core.qual = pQueryRegion->pAnchor->core.qual;

    pQueryRegion->algnType = *(pIter->pAlgnType);
    ++(pIter->pAlgnType);

    return SR_OK;
}

void SR_QueryRegionLoadSeq(SR_QueryRegion* pQueryRegion)
{
    // if we don't have enough space for the orphan sequence, we need to expand current storage space
    if (pQueryRegion->capacity < (pQueryRegion->pOrphan->core.l_qseq * 4))
    {
        pQueryRegion->capacity = pQueryRegion->pOrphan->core.l_qseq * 4;
        free(pQueryRegion->orphanSeqForward);
	free(pQueryRegion->orphanSeqReverseComplament);
	free(pQueryRegion->orphanSeqReverse);
	free(pQueryRegion->orphanSeqComplement);

        int length = pQueryRegion->pOrphan->core.l_qseq;
        pQueryRegion->orphanSeqForward           = (char*) malloc(sizeof(char) * length);
        pQueryRegion->orphanSeqReverseComplament = (char*) malloc(sizeof(char) * length);
        pQueryRegion->orphanSeqReverse           = (char*) malloc(sizeof(char) * length);
        pQueryRegion->orphanSeqComplement        = (char*) malloc(sizeof(char) * length);
        
	//if (pQueryRegion->orphanSeq == NULL)
        //    SR_ErrQuit("ERROR: Not enough memory for a orphanSeq string.");
    }

    uint8_t* seq = bam1_seq(pQueryRegion->pOrphan);
    for (unsigned int i = 0; i != pQueryRegion->pOrphan->core.l_qseq; ++i)
    {
        unsigned int f_pos, r_pos;
	char base = 0, c_base = 0;
	if (bam1_strand(pQueryRegion->pOrphan)) {// reverse-complement
	  f_pos = pQueryRegion->pOrphan->core.l_qseq - i - 1;
	  r_pos = i;
	  base   = SR_BASE_R_MAP[bam1_seqi(seq, i)];
	  c_base = SR_BASE_MAP[bam1_seqi(seq, i)];
	} else {
	  f_pos = i;
	  r_pos = pQueryRegion->pOrphan->core.l_qseq - i - 1;
	  base   = SR_BASE_MAP[bam1_seqi(seq, i)];
	  c_base = SR_BASE_R_MAP[bam1_seqi(seq, i)];
	}

        pQueryRegion->orphanSeqForward[f_pos]           = base;
	pQueryRegion->orphanSeqReverseComplament[r_pos] = c_base;
	pQueryRegion->orphanSeqReverse[r_pos]           = base;
	pQueryRegion->orphanSeqComplement[f_pos]        = c_base;
    }
    pQueryRegion->orphanSeq = pQueryRegion->orphanSeqForward;

    // reverse-complement seq and qual
    if (bam1_strand(pQueryRegion->pOrphan)) {// reverse-complement
      // reverse complement seq
      int length = pQueryRegion->pOrphan->core.l_qseq;
      uint8_t* ptr = bam1_seq(pQueryRegion->pOrphan);
      int seq_length = (length + 1) / 2;
      uint8_t* r_ptr = (uint8_t*) malloc(sizeof(uint8_t) * seq_length);
      GetReverseComplementSequence(ptr, length, r_ptr);
      memcpy(ptr, r_ptr, seq_length);
      free(r_ptr);

      // reverse qual
      uint8_t* qual_ptr = bam1_qual(pQueryRegion->pOrphan);
      GetInverseQual(qual_ptr, length);
    }
}

void SR_QueryRegionChangeSeq(SR_QueryRegion* pQueryRegion, SR_SeqAction action)
{
  switch(action) {
    case SR_FORWARD: 
      pQueryRegion->orphanSeq = pQueryRegion->orphanSeqForward;
      break;
    case SR_REVERSE_COMP:
      pQueryRegion->orphanSeq = pQueryRegion->orphanSeqReverseComplament;
      break;
    case SR_INVERSE:
      pQueryRegion->orphanSeq = pQueryRegion->orphanSeqReverse;
      break;
    case SR_COMP:
      pQueryRegion->orphanSeq = pQueryRegion->orphanSeqComplement;
      break;
    default:
      pQueryRegion->orphanSeq = pQueryRegion->orphanSeqForward;
  }
/*
    //SR_QueryRegionLoadSeq(pQueryRegion);
    
    if (action == SR_INVERSE || action == SR_REVERSE_COMP)
    {
        for (unsigned int i = 0, j = pQueryRegion->pOrphan->core.l_qseq - 1; i < j; ++i, --j)
        {
            SR_SWAP(pQueryRegion->orphanSeq[i], pQueryRegion->orphanSeq[j], char);
        }
    }

    if (action == SR_COMP || action == SR_REVERSE_COMP)
    {
        for (unsigned int i = 0; i != pQueryRegion->pOrphan->core.l_qseq; ++i)
        {
            switch(pQueryRegion->orphanSeq[i])
            {
                case 'A':
                    pQueryRegion->orphanSeq[i] = 'T';
                    break;
                case 'C':
                    pQueryRegion->orphanSeq[i] = 'G';
                    break;
                case 'G':
                    pQueryRegion->orphanSeq[i] = 'C';
                    break;
                case 'T':
                    pQueryRegion->orphanSeq[i] = 'A';
                    break;
                default:
                    pQueryRegion->orphanSeq[i] = 'N';
                    break;
            }
        }
    }
*/
}

const char* SR_QueryRegionGetSeq(const SR_QueryRegion* pQueryRegion, 
                                 const uint32_t* query_begin) {
  if (pQueryRegion == NULL) return NULL;
  if (pQueryRegion->pOrphan == NULL) return NULL;

  if (*query_begin > pQueryRegion->pOrphan->core.l_qseq) return NULL;

  return pQueryRegion->orphanSeq + *query_begin;

}

SR_Bool SR_QueryRegionSetRange(SR_QueryRegion* pQueryRegion, const SR_SearchArgs* pSearchArgs, uint32_t refLen, SR_Direction direction)
{
    if (direction == SR_DOWNSTREAM) // the position of the search region is greater than that of the anchor mate
    {
	// calculate the rightmost coordinate of an alignment on the reference genome
	uint32_t onRefEnd = bam_calend(&(pQueryRegion->pAnchor->core), bam1_cigar(pQueryRegion->pAnchor));

        pQueryRegion->closeRefBegin = onRefEnd + 1 + pSearchArgs->fragLen - pSearchArgs->closeRange / 2;
        // out of range (close/far begin is larger than the reference length of current chr)
        if (pQueryRegion->closeRefBegin >= refLen)
            return FALSE;

        // in this case, the far reference begins equals to the close reference begin
        pQueryRegion->farRefBegin = pQueryRegion->closeRefBegin;

        pQueryRegion->closeRefEnd = pQueryRegion->closeRefBegin + pSearchArgs->closeRange - 1;
        if (pQueryRegion->closeRefEnd >= refLen)
            pQueryRegion->closeRefEnd = refLen - 1;

        // search region is too small to hold a read
        if (pQueryRegion->closeRefEnd < pQueryRegion->closeRefBegin + pQueryRegion->pAnchor->core.l_qseq - 1)
            return FALSE;

        pQueryRegion->farRefEnd = pQueryRegion->farRefBegin + pSearchArgs->farRange - 1;
        if (pQueryRegion->farRefEnd >= refLen)
            pQueryRegion->farRefEnd = refLen - 1;
    }
    else if (direction == SR_UPSTREAM) // the position of the search region is less than that of the anchor mate
    {
        if (pQueryRegion->pAnchor->core.pos  + pSearchArgs->closeRange / 2 > 1 + pSearchArgs->fragLen)
            pQueryRegion->closeRefEnd = pQueryRegion->pAnchor->core.pos + pSearchArgs->closeRange / 2 - 1 - pSearchArgs->fragLen;
        else
            return FALSE; // out of range (close/far end is smaller than 0)


        // search region is too small to hold a read
        if (pQueryRegion->closeRefEnd - 0 + 1 < pQueryRegion->pAnchor->core.l_qseq)
            return FALSE;

        // in this case, the far reference end equals to the close reference end
        pQueryRegion->farRefEnd = pQueryRegion->closeRefEnd;

        if (pQueryRegion->closeRefEnd + 1 >= pSearchArgs->closeRange)
            pQueryRegion->closeRefBegin = pQueryRegion->closeRefEnd + 1 - pSearchArgs->closeRange;
        else
            pQueryRegion->closeRefBegin = 0;

        if (pQueryRegion->farRefEnd + 1 >= pSearchArgs->farRange)
            pQueryRegion->farRefBegin = pQueryRegion->farRefEnd + 1 - pSearchArgs->farRange;
        else
            pQueryRegion->farRefBegin = 0;
    }

    return TRUE;
}
