/*
 * =====================================================================================
 *
 *       Filename:  TGM_BamHeader.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  09/17/2011 02:11:18 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include "TGM_Error.h"
#include "TGM_BamHeader.h"

TGM_BamHeader* TGM_BamHeaderAlloc(void)
{
    TGM_BamHeader* pNewHeader = (TGM_BamHeader*) calloc(1, sizeof(TGM_BamHeader));
    if (pNewHeader == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for a bam header object");

    return pNewHeader;
}

void TGM_BamHeaderFree(TGM_BamHeader* pBamHeader)
{
    if (pBamHeader != NULL)
    {
        free(pBamHeader->pMD5s);
        bam_header_destroy(pBamHeader->pOrigHeader);

        free(pBamHeader);
    }
}

