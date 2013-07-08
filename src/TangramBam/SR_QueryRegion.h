/*
 * =====================================================================================
 *
 *       Filename:  SR_QueryRegion.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  08/06/2011 12:32:37 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef  SR_QUERYREGION_H
#define  SR_QUERYREGION_H

#include "../OutSources/samtools/bam.h"
#include "SR_Types.h"
#include "SR_BamInStream.h"

//===============================
// Type and constant definition
//===============================

// clear the reverse strand bit in a bam1_t structure
#define BAM_FFORWD (~BAM_FREVERSE)

// arguments used to decide the range of the query ragion
typedef struct SR_SearchArgs
{
    unsigned int fragLen;     // approximate fragment length of paired-end alignment

    unsigned int closeRange;  // region length of the one closer to the anchor mate

    unsigned int farRange;    // region length of the one further from the anchor mate

}SR_SearchArgs;

// an object that hold the unique-orphan pair and the targeting region for split alignment
typedef struct SR_QueryRegion
{
    bam1_t* pAnchor;                // the anchor alignment

    bam1_t* pOrphan;                // the orphan alignment

    char* orphanSeq;                // the pointer points to the one of following four char*
    // the sequence of the orphan alignment, not '\0' terminated
    char* orphanSeqForward;         
    char* orphanSeqReverseComplament;
    char* orphanSeqReverse;
    char* orphanSeqComplement;

    SR_AlgnType algnType;           // alignment type of the pair(unique-orphan, unique-soft or unique-multiple)

    SR_Bool isOrphanInversed;       // boolean varible used to indicate if the orphan sequence is inversed

    unsigned int capacity;          // capacity of the sequence of the orphan read in the current object

    uint32_t closeRefBegin;         // the begin position of the search region for the first partial alignment (closer to the anchor mate)

    uint32_t closeRefEnd;           // the end position of the search region for the first partial alignment (closer to the anchor mate)

    uint32_t farRefBegin;           // the begin position of search region for the second partial alignment (further from the anchor mate)

    uint32_t farRefEnd;             // the end position of search region for the second partial alignment (further from the anchor mate)

}SR_QueryRegion;


//===============================
// Constructors and Destructors
//===============================

SR_QueryRegion* SR_QueryRegionAlloc(void);

void SR_QueryRegionFree(SR_QueryRegion* pQueryRegion);


//==================
// Inline functions
//==================

// get and set the length of an alignment
#define SR_GetQueryLen(pAlignment) ((pAlignment)->core.l_qseq)
#define SR_SetQueryLen(pAlignment, queryLen) ((pAlignment)->core.l_qseq = queryLen)

//=========================================================
// function: 
//      get the strand direction of an alignment
//
// args:
//      1. pAlignment: a pointer to an alignment structure
//
//return:
//      the strand of the alignment
//=========================================================
static inline SR_Strand SR_GetStrand(const bam1_t* const pAlignment)
{
    if (bam1_strand(pAlignment))
        return SR_REVERSE_COMP;
    else
        return SR_FORWARD;
}

//==========================================================
// function:
//      set the strand direction of an alignment
//
// args:
//      1. pAlignment: a pointer to an alignment structure
//      2. strand    : strand direction
//==========================================================
static inline void SR_SetStrand(bam1_t* pAlignment, SR_Strand strand)
{
    if (strand == SR_REVERSE_COMP)
        pAlignment->core.flag |= BAM_FREVERSE;
    else
        pAlignment->core.flag &= BAM_FFORWD;
}


//======================
// Interface functions
//======================

//================================================================
// function:
//      load the unique-orphan pair into the query region 
//      object
//
// args:
//      1. pQueryRegion: a pointer to a query region object
//      2. iter: a iterator to an alignment buffer
//================================================================ 
SR_Status SR_QueryRegionLoadPair(SR_QueryRegion* pQueryRegion, SR_BamInStreamIter* pIter);

//==============================================================
// function:
//      transfer the DNA sequence in bam structure from the
//      4-bits representation to the ascii representation and
//
// args:
//      1. pQueryRegion: a pointer to an query region structure
//==============================================================
void SR_QueryRegionLoadSeq(SR_QueryRegion* pQueryRegion);

//==============================================================
// function:
//      change the orphan sequence according to the action 
//      argument
//
// args:
//      1. pQueryRegion: a pointer to an query region structure
//      2. action      : action applied on the DNA sequence
//
// discussion:
//      the strand of the orphan mate will not be set.
//      you have to set it through the "SR_SetStrand"
//==============================================================
void SR_QueryRegionChangeSeq(SR_QueryRegion* pQueryRegion, SR_SeqAction action);


//==============================================================
// function:
//      set the search region for the split aligner according
//      to the position of anchor mate and user specified
//      parameters
//
// args:
//      1. pQueryRegion: a pointer to an query region structure
//      2. pSearchArgs : a pointer to an search arguments
//                      structure containing the user specific
//                      search arguments
//      3. refLen     : length of the current chromosome
//      4. direction  : search direction relative to the anchor
//                      position, upstream or downstream
// return:
//      if the region is successfully set, return TRUE; if the
//      region can not be set due to the range limit of chr it 
//      will return FALSE
//==============================================================
SR_Bool SR_QueryRegionSetRange(SR_QueryRegion* pQueryRegion, const SR_SearchArgs* pSearchArgs, uint32_t refLen, SR_Direction direction);

//==============================================================
// function:
//      set the search region on special reference for the 
//      split aligner
//
// args:
//      1. pQueryRegion: a pointer to an query region structure
//      2. specialRefLen : length of the special references
//==============================================================
static inline void SR_QueryRegionSetRangeSpecial(SR_QueryRegion* pQueryRegion, unsigned int specialRefLen) 
{
    pQueryRegion->closeRefBegin = 0;
    pQueryRegion->closeRefEnd = specialRefLen - 1;

    pQueryRegion->farRefBegin = 0;
    pQueryRegion->farRefEnd = specialRefLen - 1;
}

#endif  /*SR_QUERYREGION_H*/


