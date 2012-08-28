/*
 * =====================================================================================
 *
 *       Filename:  TGM_BamHeader.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  09/15/2011 04:21:11 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef  TGM_BAMHEADER_H
#define  TGM_BAMHEADER_H

#include "bam.h"

//===============================
// Type and constant definition
//===============================

// stucture holding header information
typedef struct TGM_BamHeader
{
    bam_header_t* pOrigHeader;

    const char** pMD5s;

}TGM_BamHeader;


//===============================
// Constructors and Destructors
//===============================

TGM_BamHeader* TGM_BamHeaderAlloc(void);

void TGM_BamHeaderFree(TGM_BamHeader* pBamHeader);


//======================
// Inline functions
//======================

//===============================================================
// function:
//      get the number of references stored in the bam header
//
// args:
//      1. pBamHeader: a pointer to the header structure
// 
// return:
//      number of references (chromosomes)
//=============================================================== 
inline int32_t TGM_BamHeaderGetRefNum(const TGM_BamHeader* pBamHeader)
{
    return (pBamHeader->pOrigHeader->n_targets);
}

//===============================================================
// function:
//      get the reference ID to reference name dictionary
//
// args:
//      1. pBamHeader: a pointer to the header structure
// 
// return:
//      the dictionary of reference ID to reference name
//=============================================================== 
inline const char** TGM_BamHeaderGetRefNames(const TGM_BamHeader* pBamHeader)
{
    return (const char**) pBamHeader->pOrigHeader->target_name;
}

//===============================================================
// function:
//      get the array of reference length
//
// args:
//      1. pBamHeader: a pointer to the header structure
// 
// return:
//      an array contains the length of each chromosome
//=============================================================== 
inline const uint32_t* TGM_BamHeaderGetRefLens(const TGM_BamHeader* pBamHeader)
{
    return pBamHeader->pOrigHeader->target_len;
}


#endif  /*TGM_BAMHEADER_H*/
