/*
 * =====================================================================================
 *
 *       Filename:  SR_BamHeader.c
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

#include "SR_Error.h"
#include "SR_BamHeader.h"

SR_BamHeader* SR_BamHeaderAlloc(void)
{
    SR_BamHeader* pNewHeader = (SR_BamHeader*) calloc(1, sizeof(SR_BamHeader));
    if (pNewHeader == NULL)
        SR_ErrQuit("ERROR: Not enough memory for a bam header object");

    return pNewHeader;
}

void SR_BamHeaderFree(SR_BamHeader* pBamHeader)
{
    if (pBamHeader != NULL)
    {
        free(pBamHeader->pMD5s);
        bam_header_destroy(pBamHeader->pOrigHeader);

        free(pBamHeader);
    }
}

int32_t SR_BamHeaderGetRefNum(const SR_BamHeader* pBamHeader)
{
    return (pBamHeader->pOrigHeader->n_targets);
}

const char** SR_BamHeaderGetRefNames(const SR_BamHeader* pBamHeader)
{
    return (const char**) pBamHeader->pOrigHeader->target_name;
}

const uint32_t* SR_BamHeaderGetRefLens(const SR_BamHeader* pBamHeader)
{
    return pBamHeader->pOrigHeader->target_len;
}

