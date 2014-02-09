/*
 * =====================================================================================
 *
 *       Filename:  TGM_BamPairAux.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  09/17/2011 09:43:53 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef  TGM_BAMPAIRAUX_H
#define  TGM_BAMPAIRAUX_H


#include "bam.h"
#include "TGM_Types.h"
#include "TGM_Utilities.h"
#include "TGM_LibInfo.h"


//===============================
// Type and constant definition
//===============================

static const unsigned int TGM_UNIQUE_ORPHAN_FMASK = (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP);

static const unsigned int TGM_NORMAL_FMASK = (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP | BAM_FMUNMAP);

static const unsigned int TGM_READ_PAIR_FMASK = (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP | BAM_FMUNMAP);

/*

typedef struct TGM_FilterDataRP
{
    bam1_t* pUpAlgn;

    bam1_t* pDownAlgn;

    bamFile pBamInput;

    bam_index_t* pBamIndex;

    const TGM_AnchorInfo* pAnchorInfo;

    uint32_t binLen;

    TGM_Bool loadCross;

    TGM_Bool isFilled;

}TGM_FilterDataRP;


TGM_FilterDataRP* TGM_FilterDataRPAlloc(const TGM_AnchorInfo* pAnchorInfo, uint32_t binLen);

void TGM_FilterDataRPFree(TGM_FilterDataRP* pFilterData);
*/

typedef struct TGM_FilterDataNoZA
{
    const TGM_LibInfoTable* pLibTable;

    TGM_Bool keepNormal;

}TGM_FilterDataNoZA;

//======================
// Interface functions
//======================

static inline TGM_StreamCode TGM_CommonFilter(const bam1_t* pAlignment, void* pFilterData)
{
    if ((pAlignment->core.flag & BAM_FPAIRED) == 0
        || pAlignment->core.tid < 0
        || pAlignment->core.mtid < 0
        || strcmp(bam1_qname(pAlignment), "*") == 0
        || strcmp(bam1_qname(pAlignment), "0") == 0
        || (pAlignment->core.flag & TGM_UNIQUE_ORPHAN_FMASK) != 0
        || (pAlignment->core.flag | (BAM_FUNMAP | BAM_FMUNMAP)) == pAlignment->core.flag)
    {
        return STREAM_PASS;
    }
    
    return STREAM_KEEP;
}

static inline TGM_StreamCode TGM_NormalFilter(const bam1_t* pAlignment, void* pFilterData)
{
    if ((pAlignment->core.flag & BAM_FPAIRED) == 0
        || pAlignment->core.tid < 0
        || pAlignment->core.mtid < 0
        || strcmp(bam1_qname(pAlignment), "*") == 0
        || strcmp(bam1_qname(pAlignment), "0") == 0
        || (pAlignment->core.flag & TGM_NORMAL_FMASK) != 0
        || (pAlignment->core.isize == 0))
    {
        return STREAM_PASS;
    }

    // any reads aligned to different chromosome will be kept as SV candidates
    if (pAlignment->core.tid != pAlignment->core.mtid)
        return STREAM_PASS;

    return STREAM_KEEP;
}

static inline TGM_StreamCode TGM_ReadPairFilter(const bam1_t* pAlignment, void* pFilterData)
{
    TGM_Bool* pLoadCross = (TGM_Bool*) pFilterData; 

    if ((pAlignment->core.flag & BAM_FPAIRED) == 0
        || pAlignment->core.tid < 0
        || pAlignment->core.mtid < 0
        || strcmp(bam1_qname(pAlignment), "0") == 0
        || strcmp(bam1_qname(pAlignment), "*") == 0
        || (pAlignment->core.flag & TGM_READ_PAIR_FMASK) != 0)
    {
        return STREAM_PASS;
    }

    // any reads aligned to different chromosome will be kept as SV candidates
    if (!(*pLoadCross) && pAlignment->core.tid != pAlignment->core.mtid) 
        return STREAM_PASS;

    return STREAM_KEEP;
}

TGM_StreamCode TGM_ReadPairNoZAFilter(const bam1_t* pAlignment, void* pFilterData);

static inline TGM_StreamCode TGM_SplitFilter(const bam1_t* pAlignment, void* pFilterData)
{
    if (strcmp(bam1_qname(pAlignment), "*") == 0
        || strcmp(bam1_qname(pAlignment), "0") == 0
        || pAlignment->core.tid < 0
        || pAlignment->core.mtid < 0
        || (pAlignment->core.flag & TGM_UNIQUE_ORPHAN_FMASK) != 0
        || (pAlignment->core.flag | (BAM_FUNMAP | BAM_FMUNMAP)) == pAlignment->core.flag)
    {
        return STREAM_PASS;
    }
    
    return STREAM_KEEP;
}

// void TGM_FilterDataRPInit(TGM_FilterDataRP* pFilterData, const char* bamInputFile);

#define TGM_FilterDataRPTurnOffCross(pFilterData) (pFilterData)->loadCross = FALSE

#define TGM_FilterDataRPTurnOnCross(pFilterData) (pFilterData)->loadCross = TRUE



    
#endif  /*TGM_BAMPAIRAUX_H*/
