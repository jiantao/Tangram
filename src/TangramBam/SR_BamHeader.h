/*
 * =====================================================================================
 *
 *       Filename:  SR_BamHeader.h
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

#ifndef  SR_BAMHEADER_H
#define  SR_BAMHEADER_H


#include "../OutSources/samtools/bam.h"


//===============================
// Type and constant definition
//===============================

// stucture holding header information
typedef struct SR_BamHeader
{
    bam_header_t* pOrigHeader;

    const char** pMD5s;

}SR_BamHeader;


//===============================
// Constructors and Destructors
//===============================

SR_BamHeader* SR_BamHeaderAlloc(void);

void SR_BamHeaderFree(SR_BamHeader* pBamHeader);


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
inline int32_t SR_BamHeaderGetRefNum(const SR_BamHeader* pBamHeader);

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
inline const char** SR_BamHeaderGetRefNames(const SR_BamHeader* pBamHeader);

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
inline const uint32_t* SR_BamHeaderGetRefLens(const SR_BamHeader* pBamHeader);

#endif  /*SR_BAMHEADER_H*/
