/*
 * =====================================================================================
 *
 *       Filename:  TGM_LibInfo.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  02/04/2012 08:37:12 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include <ctype.h>
#include <dirent.h>
#include <sys/stat.h>

#include "khash.h"
#include "TGM_Error.h"
#include "TGM_Types.h"
#include "TGM_LibInfo.h"
#include "TGM_Utilities.h"

// index of the mode count for invalid pair mode 
#define INVALID_PAIR_MODE_SET_INDEX 1

#define NUM_ZA_FIELD 7

static const char* TGM_LibTableFileName = "lib_table.dat";

static const char* TGM_HistFileName = "hist.dat";

// sample name hash
KHASH_MAP_INIT_STR(name, uint32_t);

KHASH_SET_INIT_INT(key);

void TGM_GetNumMismatchFromZA(int16_t* pNumMM, int32_t* pLen, const char* cigarStr, unsigned int cigarLen, const char* mdStr, unsigned int mdLen)
{
    if (cigarStr == NULL)
        return;

    int numMM = 0;
    int len = 0;
    for (const char* cigarEnd = strpbrk(cigarStr, "DIM"); cigarEnd != NULL && cigarEnd - cigarStr < cigarLen; cigarEnd = strpbrk(cigarEnd + 1, "DIM"))
    {
        if (*cigarEnd != 'M')
        {
            const char* currPos = cigarEnd - 1;
            while (isdigit(*currPos))
                --currPos;

            int mm = atoi(currPos + 1);
            numMM += mm;
            
            if (*cigarEnd == 'D')
                len += mm;
        }
        else
        {
            const char* currPos = cigarEnd - 1;
            while (isdigit(*currPos))
                --currPos;
            
            len += atoi(currPos + 1);
        }
    }

    const char* mdFieldPos = mdStr;
    while (mdFieldPos != NULL && mdFieldPos - mdStr < mdLen)
    {
        if (isdigit(*mdFieldPos))
        {
            ++mdFieldPos;
            continue;
        }

        const char* mdFieldEnd = mdFieldPos + 1;
        while (!isdigit(*mdFieldEnd) && *mdFieldEnd != '\0')
            ++mdFieldEnd;

        if (*mdFieldPos != '^')
            numMM += mdFieldEnd - mdFieldPos;

        mdFieldPos = mdFieldEnd;
    }

    *pNumMM = numMM;
    *pLen = len;
}

static TGM_PairMode TGM_GetPairMode(const bam1_t* pAlignment)
{
    int8_t upMode = 0;
    int8_t downMode = 0;

    if ((pAlignment->core.flag & BAM_FREVERSE) != 0)
    {
        upMode |= 1;
    }

    if ((pAlignment->core.flag & BAM_FMREVERSE) != 0)
    {
        downMode |= 1;
    }

    if ((pAlignment->core.flag & BAM_FREAD1) != 0 && (pAlignment->core.flag & BAM_FREAD2) != 0)
        return TGM_BAD_PAIR_MODE;
    else if ((pAlignment->core.flag & BAM_FREAD2) != 0)
        upMode |= (1 << 1);
    else if ((pAlignment->core.flag & BAM_FREAD1) != 0)
        downMode |= (1 << 1);
    else
        return TGM_BAD_PAIR_MODE;

    if (pAlignment->core.tid > pAlignment->core.mtid)
        TGM_SWAP(upMode, downMode, unsigned int);
    else if(pAlignment->core.tid == pAlignment->core.mtid)
    {
        if (pAlignment->core.pos > pAlignment->core.mpos)
            TGM_SWAP(upMode, downMode, unsigned int);
    }

    return (TGM_PairMode) TGM_PairModeMap[((upMode << 2) | downMode)];
}

/*
static TGM_Status TGM_LibInfoTableAddAnchor(TGM_LibInfoTable* pTable, const char* tagPos)
{
    const TGM_Bool isLoaded = pTable->pAnchorInfo->size == 0 ? FALSE : TRUE;
    TGM_AnchorInfo* pAnchorInfo = pTable->pAnchorInfo;

    char buff[255];
    while ((tagPos = strstr(tagPos, "@SQ")) != NULL)
    {
        const char* lineEnd = strpbrk(tagPos, "\n");
        const char* refNamePos = strstr(tagPos, "SN:");

        if (refNamePos == NULL || refNamePos > lineEnd)
        {
            TGM_ErrMsg("ERROR: Reference name in the bam header is not specified\n");
            return TGM_ERR;
        }

        refNamePos += 3;
        const char* refNameEnd = strpbrk(refNamePos, " \t\n");
        if (refNameEnd == NULL || refNameEnd > lineEnd)
        {
            TGM_ErrMsg("ERROR: Reference name in the bam header is not specified\n");
            return TGM_ERR;
        }

        int refNameLen = refNameEnd - refNamePos;
        if (refNameLen <= 0)
        {
            TGM_ErrMsg("ERROR: Reference name in the bam header is not specified\n");
            return TGM_ERR;
        }

        int ret = 0;
        khiter_t khIter = 0;
        unsigned int refIndex = 0;
        TGM_Bool isUsedRef = TRUE;

        if (isLoaded)
        {
            buff[refNameLen] = '\0';
            strncpy(buff, refNamePos, refNameLen);

            khIter = kh_put(name, pAnchorInfo->pAnchorHash, buff, &ret);
            if (ret != 0)
            {
                TGM_ErrMsg("ERROR: Found a reference name that is not recorded in the previous bam files.\n");
                return TGM_ERR;
            }

            refIndex = kh_value((khash_t(name)*) pAnchorInfo->pAnchorHash, khIter);
            if (pAnchorInfo->pLength[refIndex] < 0)
                isUsedRef = FALSE;
        }
        else
        {
            if (pAnchorInfo->size == pAnchorInfo->capacity)
            {
                pAnchorInfo->capacity *= 2;
                pAnchorInfo->pAnchors = (char**) realloc(pAnchorInfo->pAnchors, sizeof(char*) * pAnchorInfo->capacity);
                if (pAnchorInfo->pAnchors == NULL)
                    TGM_ErrQuit("ERROR: Not enough memory for the storage of the anchor names.\n");

                pAnchorInfo->pLength = (int32_t*) realloc(pAnchorInfo->pLength, sizeof(int32_t) * pAnchorInfo->capacity);
                if (pAnchorInfo->pLength == NULL)
                    TGM_ErrQuit("ERROR: Not enough memory for the storage of the anchor length.\n");

                pAnchorInfo->pMd5s = (char*) realloc(pAnchorInfo->pMd5s, sizeof(char) * MD5_STR_LEN * pAnchorInfo->capacity);
                if (pAnchorInfo->pMd5s == NULL)
                    TGM_ErrQuit("ERROR: Not enough memory for the storage of the md5 strings.\n");
            }

            pAnchorInfo->pAnchors[pAnchorInfo->size] = (char*) malloc((refNameLen + 1) * sizeof(char));
            if (pAnchorInfo->pAnchors[pAnchorInfo->size] == NULL)
                TGM_ErrQuit("ERROR: Not enough memory for the storage of reference name.\n");

            pAnchorInfo->pAnchors[pAnchorInfo->size][refNameLen] = '\0';
            strncpy(pAnchorInfo->pAnchors[pAnchorInfo->size], refNamePos, refNameLen);

            khIter = kh_put(name, pAnchorInfo->pAnchorHash, pAnchorInfo->pAnchors[pAnchorInfo->size], &ret);
            if (ret == 0)
            {
                TGM_ErrMsg("ERROR: Found a duplicated reference name in a bam file.\n");
                return TGM_ERR;
            }

            refIndex = pAnchorInfo->size;
            kh_value((khash_t(name)*) pAnchorInfo->pAnchorHash, khIter) = refIndex;

            if (strncmp("GL0", pAnchorInfo->pAnchors[refIndex], 3) == 0
                || strncmp("NC_", pAnchorInfo->pAnchors[refIndex], 3) == 0
                || strncmp("NT_", pAnchorInfo->pAnchors[refIndex], 3) == 0)
            {
                isUsedRef = FALSE;
                pAnchorInfo->pLength[refIndex] = -1;
            }
        }

        if (isUsedRef)
        {
            const char* refLenPos = strstr(tagPos, "LN:");
            if (refLenPos == NULL || refLenPos > lineEnd)
            {
                TGM_ErrMsg("ERROR: Cannot find the length of the reference.\n");
                return TGM_ERR;
            }

            refLenPos += 3;
            const char* refLenEnd = strpbrk(refLenPos, " \t\n");
            int refLenLen = refLenEnd - refLenPos;
            if (refLenEnd == NULL || refLenLen <= 0 || refLenEnd > lineEnd)
            {
                TGM_ErrMsg("ERROR: Cannot find the length of the reference.\n");
                return TGM_ERR;
            }

            buff[refLenLen] = '\0';
            strncpy(buff, refLenPos, refLenLen);

            int refLen = atoi(buff);
            if (refLen == 0)
            {
                TGM_ErrMsg("ERROR: Cannot find the length of the reference.\n");
                return TGM_ERR;
            }

            if (isLoaded)
            {
                if (pAnchorInfo->pLength[refIndex] > 0 && refLen != pAnchorInfo->pLength[refIndex])
                {
                    TGM_ErrMsg("ERROR: The length of the reference found in this bam file is inconsistent with that in the previous bam file.\n");
                    return TGM_ERR;
                }
            }
            else
                pAnchorInfo->pLength[refIndex] = refLen;

            const char* md5Pos = strstr(tagPos, "M5:");
            if (md5Pos == NULL || md5Pos > lineEnd)
            {
                TGM_ErrMsg("ERROR: Cannot find the md5 string of the reference.\n");
                return TGM_ERR;
            }

            md5Pos += 3;
            if (isLoaded)
            {
                strncpy(buff, md5Pos, MD5_STR_LEN);
                if (strncmp(buff, pAnchorInfo->pMd5s + refIndex * MD5_STR_LEN, MD5_STR_LEN) != 0)
                {
                    TGM_ErrMsg("ERROR: The md5 string of the reference found in this bam file is inconsistent with that in the previous bam file.\n");
                    return TGM_ERR;
                }
            }
            else
                strncpy(pAnchorInfo->pMd5s + refIndex * MD5_STR_LEN, md5Pos, MD5_STR_LEN);

        }
        else if (!isLoaded)
            memset(pAnchorInfo->pMd5s + refIndex * MD5_STR_LEN, '-', MD5_STR_LEN);

        if (!isLoaded)
            ++(pAnchorInfo->size);

        ++tagPos;
    }

    return TGM_OK;
}

*/

static TGM_Status TGM_LibInfoTableAddAnchor(TGM_LibInfoTable* pTable, const TGM_BamHeader* pBamHeader)
{
    const TGM_Bool isLoaded = pTable->pAnchorInfo->size == 0 ? FALSE : TRUE;

    int ret = 0;
    khiter_t khIter = 0;

    if (isLoaded)
    {
        if (pBamHeader->pOrigHeader->n_targets != pTable->pAnchorInfo->size)
        {
            TGM_ErrMsg("ERROR: The number of reference sequences in this bam file is inconsistent with that in previous bam files.\n");
            return TGM_ERR;
        }

        unsigned int refIndex = 0;
        for (unsigned int i = 0; i != pBamHeader->pOrigHeader->n_targets; ++i)
        {
            khIter = kh_put(name, pTable->pAnchorInfo->pAnchorHash, pBamHeader->pOrigHeader->target_name[i], &ret);
            if (ret == 0)
            {
                khash_t(name)* pAnchorHash = pTable->pAnchorInfo->pAnchorHash;
                refIndex = kh_value(pAnchorHash, khIter);
                if (refIndex != i)
                {
                    TGM_ErrMsg("ERROR: Reference ID in this bam file is inconsistent with that in the previous bam files.\n");
                    return TGM_ERR;
                }

                if (pTable->pAnchorInfo->pLength[i] > 0 && pTable->pAnchorInfo->pLength[i] != pBamHeader->pOrigHeader->target_len[i])
                {
                    TGM_ErrMsg("ERROR: The length of the reference sequence in this bam file is inconsistent with that in previous bam files.\n");
                    return TGM_ERR;
                }
            }
            else
            {
                TGM_ErrMsg("ERROR: Found a reference sequence that is not in the previous bam headers.\n");
                return TGM_ERR;
            }

            if (strncmp(pTable->pAnchorInfo->pMd5s + MD5_STR_LEN * i, pBamHeader->pMD5s[i], MD5_STR_LEN) != 0)
            {
                TGM_ErrMsg("ERROR: MD5 string in this bam file is inconsistent with that in previous bam files.\n");
                return TGM_ERR;
            }
        }
    }
    else
    {
        if (pBamHeader->pOrigHeader->n_targets > pTable->pAnchorInfo->capacity)
        {
            const char* pSpecialPrefix = pTable->pAnchorInfo->pSpecialPrefix;

            TGM_AnchorInfoFree(pTable->pAnchorInfo);
            pTable->pAnchorInfo = TGM_AnchorInfoAlloc(pBamHeader->pOrigHeader->n_targets, pSpecialPrefix);
        }

        pTable->pAnchorInfo->size = pBamHeader->pOrigHeader->n_targets;

        for (unsigned int i = 0; i != pBamHeader->pOrigHeader->n_targets; ++i)
        {
            unsigned int nameLen = strlen(pBamHeader->pOrigHeader->target_name[i]);

            pTable->pAnchorInfo->pAnchors[i] = (char*) calloc(nameLen + 1, sizeof(char));
            if (pTable->pAnchorInfo->pAnchors[i] == NULL)
                TGM_ErrQuit("ERROR: Not enough memory for the storage of reference name.\n");

            strncpy(pTable->pAnchorInfo->pAnchors[i], pBamHeader->pOrigHeader->target_name[i], nameLen);
            khIter = kh_put(name, pTable->pAnchorInfo->pAnchorHash, pTable->pAnchorInfo->pAnchors[i], &ret);

            khash_t(name)* pAnchorHash = pTable->pAnchorInfo->pAnchorHash;
            kh_value(pAnchorHash, khIter) = i;

            // set the length of those unused references to -1 so that
            // we will ignore any reads aligned to them
            pTable->pAnchorInfo->pLength[i] = pBamHeader->pOrigHeader->target_len[i];
            if (strncmp("GL0", pTable->pAnchorInfo->pAnchors[i], 3) == 0
                || strncmp("NC_", pTable->pAnchorInfo->pAnchors[i], 3) == 0
                || strncmp("NT_", pTable->pAnchorInfo->pAnchors[i], 3) == 0
                || strncmp("hs", pTable->pAnchorInfo->pAnchors[i], 2) == 0)
            {
                pTable->pAnchorInfo->pLength[i] = -1;
            }

            strncpy(pTable->pAnchorInfo->pMd5s + MD5_STR_LEN * i, pBamHeader->pMD5s[i], MD5_STR_LEN);
        }
    }

    return TGM_OK;
}

static TGM_Status TGM_LibInfoTableAddAnchorSplit(TGM_LibInfoTable* pTable, const TGM_BamHeader* pBamHeader, const char* specialPrefix, unsigned int prefixLen)
{
    const TGM_Bool isLoaded = pTable->pAnchorInfo->size == 0 ? FALSE : TRUE;

    int ret = 0;
    khiter_t khIter = 0;

    if (isLoaded)
    {
        unsigned int refIndex = 0;
        for (unsigned int i = 0; i != pBamHeader->pOrigHeader->n_targets; ++i)
        {
            khIter = kh_put(name, pTable->pAnchorInfo->pAnchorHash, pBamHeader->pOrigHeader->target_name[i], &ret);
            if (ret == 0)
            {
                khash_t(name)* pAnchorHash = pTable->pAnchorInfo->pAnchorHash;
                refIndex = kh_value(pAnchorHash, khIter);
                if (refIndex != i)
                {
                    TGM_ErrMsg("ERROR: Reference ID in this bam file is inconsistent with that in the previous bam files.\n");
                    return TGM_ERR;
                }

                if (pTable->pAnchorInfo->pLength[i] > 0 && pTable->pAnchorInfo->pLength[i] != pBamHeader->pOrigHeader->target_len[i])
                {
                    TGM_ErrMsg("ERROR: The length of the reference sequence in this bam file is inconsistent with that in previous bam files.\n");
                    return TGM_ERR;
                }

                if (strncmp(pTable->pAnchorInfo->pMd5s + MD5_STR_LEN * i, pBamHeader->pMD5s[i], MD5_STR_LEN) != 0)
                {
                    TGM_ErrMsg("ERROR: MD5 string in this bam file is inconsistent with that in previous bam files.\n");
                    return TGM_ERR;
                }
            }
            else
            {
                // if there are already reference records in the anchor info table
                // then only special reference is allowed to be added into the table
                if (strncmp(pBamHeader->pOrigHeader->target_name[i], specialPrefix, prefixLen) != 0)
                {
                    TGM_ErrMsg("ERROR: Found a reference sequence that is not in the previous bam headers.\n");
                    return TGM_ERR;
                }

                unsigned int nameLen = strlen(pBamHeader->pOrigHeader->target_name[i]);

                pTable->pAnchorInfo->pAnchors[i] = (char*) calloc(nameLen + 1, sizeof(char));
                if (pTable->pAnchorInfo->pAnchors[i] == NULL)
                    TGM_ErrQuit("ERROR: Not enough memory for the storage of reference name.\n");

                strncpy(pTable->pAnchorInfo->pAnchors[i], pBamHeader->pOrigHeader->target_name[i], nameLen);

                khash_t(name)* pAnchorHash = pTable->pAnchorInfo->pAnchorHash;
                kh_value(pAnchorHash, khIter) = i;

                // set the length of those unused references to -1 so that
                // we will ignore any reads aligned to them
                pTable->pAnchorInfo->pLength[i] = pBamHeader->pOrigHeader->target_len[i];
                if (strncmp("GL0", pTable->pAnchorInfo->pAnchors[i], 3) == 0
                    || strncmp("NC_", pTable->pAnchorInfo->pAnchors[i], 3) == 0
                    || strncmp("NT_", pTable->pAnchorInfo->pAnchors[i], 3) == 0
                    || strncmp("hs", pTable->pAnchorInfo->pAnchors[i], 2) == 0)
                {
                    pTable->pAnchorInfo->pLength[i] = -1;
                }

                strncpy(pTable->pAnchorInfo->pMd5s + MD5_STR_LEN * i, pBamHeader->pMD5s[i], MD5_STR_LEN);
            }
        }
    }
    else
    {
        if (pBamHeader->pOrigHeader->n_targets > pTable->pAnchorInfo->capacity)
        {
            const char* pSpecialPrefix = pTable->pAnchorInfo->pSpecialPrefix;

            TGM_AnchorInfoFree(pTable->pAnchorInfo);
            pTable->pAnchorInfo = TGM_AnchorInfoAlloc(pBamHeader->pOrigHeader->n_targets, pSpecialPrefix);
        }

        pTable->pAnchorInfo->size = pBamHeader->pOrigHeader->n_targets;

        for (unsigned int i = 0; i != pBamHeader->pOrigHeader->n_targets; ++i)
        {
            unsigned int nameLen = strlen(pBamHeader->pOrigHeader->target_name[i]);

            pTable->pAnchorInfo->pAnchors[i] = (char*) calloc(nameLen + 1, sizeof(char));
            if (pTable->pAnchorInfo->pAnchors[i] == NULL)
                TGM_ErrQuit("ERROR: Not enough memory for the storage of reference name.\n");

            strncpy(pTable->pAnchorInfo->pAnchors[i], pBamHeader->pOrigHeader->target_name[i], nameLen);
            khIter = kh_put(name, pTable->pAnchorInfo->pAnchorHash, pTable->pAnchorInfo->pAnchors[i], &ret);

            khash_t(name)* pAnchorHash = pTable->pAnchorInfo->pAnchorHash;
            kh_value(pAnchorHash, khIter) = i;

            // set the length of those unused references to -1 so that
            // we will ignore any reads aligned to them
            pTable->pAnchorInfo->pLength[i] = pBamHeader->pOrigHeader->target_len[i];
            if (strncmp("GL0", pTable->pAnchorInfo->pAnchors[i], 3) == 0
                || strncmp("NC_", pTable->pAnchorInfo->pAnchors[i], 3) == 0
                || strncmp("NT_", pTable->pAnchorInfo->pAnchors[i], 3) == 0
                || strncmp("hs", pTable->pAnchorInfo->pAnchors[i], 2) == 0)
            {
                pTable->pAnchorInfo->pLength[i] = -1;
            }

            strncpy(pTable->pAnchorInfo->pMd5s + MD5_STR_LEN * i, pBamHeader->pMD5s[i], MD5_STR_LEN);
        }
    }

    return TGM_OK;
}

static TGM_Status TGM_AnchorInfoMerge(TGM_AnchorInfo* pDstAnchor, const TGM_AnchorInfo* pSrcAnchor)
{
    int ret = 0;
    khiter_t khIter = 0;
    unsigned int refIndex = 0;

    for (unsigned int i = 0; i != pSrcAnchor->size; ++i)
    {
        khIter = kh_put(name, pDstAnchor->pAnchorHash, pSrcAnchor->pAnchors[i], &ret);
        if (ret == 0)
        {
            khash_t(name)* pAnchorHash = pDstAnchor->pAnchorHash;
            refIndex = kh_value(pAnchorHash, khIter);
            if (refIndex != i)
            {
                TGM_ErrMsg("ERROR: Reference ID in this bam file is inconsistent with that in the previous bam files.\n");
                return TGM_ERR;
            }

            if (pDstAnchor->pLength[i] > 0 && pDstAnchor->pLength[i] != pSrcAnchor->pLength[i])
            {
                TGM_ErrMsg("ERROR: The length of the reference sequence in this bam file is inconsistent with that in previous bam files.\n");
                return TGM_ERR;
            }

            if (strncmp(pDstAnchor->pMd5s + MD5_STR_LEN * i, pSrcAnchor->pMd5s + MD5_STR_LEN * i, MD5_STR_LEN) != 0)
            {
                TGM_ErrMsg("ERROR: MD5 string in this bam file is inconsistent with that in previous bam files.\n");
                return TGM_ERR;
            }
        }
        else
            return TGM_ERR;
    }

    return TGM_OK;
}

TGM_Status TGM_SampleInfoMerge(TGM_SampleInfo* pDstSampleInfo, unsigned int* pIndexMap, TGM_SampleInfo* pSrcSampleInfo)
{
    for (unsigned int i = 0; i != pSrcSampleInfo->size; ++i)
    {
        int ret = 0;

        khash_t(name)* pSampleHash = pDstSampleInfo->pSampleHash;
        khiter_t khIter = kh_put(name, pSampleHash, pSrcSampleInfo->pSamples[i], &ret);

        if (ret != 0)
        {
            if (pDstSampleInfo->size == pDstSampleInfo->capacity)
            {
                pDstSampleInfo->capacity *= 2;
                pDstSampleInfo->pSamples = (char**) realloc(pDstSampleInfo->pSamples, pDstSampleInfo->capacity * sizeof(char*));
                if (pDstSampleInfo->pSamples == NULL)
                    TGM_ErrQuit("ERROR: Not enought memory for the sample names in the fragment length distribution object.\n");
            }

            kh_value(pSampleHash, khIter) = pDstSampleInfo->size;
            pDstSampleInfo->pSamples[pDstSampleInfo->size] = pSrcSampleInfo->pSamples[i];
            pSrcSampleInfo->pSamples[i] = NULL;
            ++(pDstSampleInfo->size);
        }

        pIndexMap[i] = kh_value(pSampleHash, khIter);
    }

    return TGM_OK;
}

static TGM_Status TGM_LibInfoTableAddSample(int* pSampleID, TGM_LibInfoTable* pTable, const char* tagPos, const char* lineEnd)
{
    const char* sampleNamePos = strstr(tagPos, "SM:");
    if (sampleNamePos != NULL && sampleNamePos < lineEnd)
        sampleNamePos += 3;

    const char* sampleNameEnd = strpbrk(sampleNamePos, " \t\n");
    int sampleNameLen = sampleNameEnd - sampleNamePos;

    TGM_SampleInfo* pSampleInfo = pTable->pSampleInfo;

    if (sampleNameLen > 0)
    {
        char* buff = (char*) malloc((sampleNameLen + 1) * sizeof(char));
        if (buff == NULL)
            TGM_ErrQuit("ERROR: Not enought memory for the read group ID in the fragment length distribution object.\n");

        buff[sampleNameLen] = '\0';
        memcpy(buff, sampleNamePos, sampleNameLen);

        int ret = 0;

        khash_t(name)* pSampleHash = pSampleInfo->pSampleHash;
        khiter_t khIter = kh_put(name, pSampleHash, buff, &ret);

        if (ret == 0)
        {
            free(buff);
            *pSampleID = kh_value(pSampleHash, khIter);
        }
        else
        {
            if (pSampleInfo->size == pSampleInfo->capacity)
            {
                pSampleInfo->capacity *= 2;
                pSampleInfo->pSamples = (char**) realloc(pSampleInfo->pSamples, pSampleInfo->capacity * sizeof(char*));
                if (pSampleInfo->pSamples == NULL)
                    TGM_ErrQuit("ERROR: Not enought memory for the sample names in the fragment length distribution object.\n");
            }

            *pSampleID = pSampleInfo->size;
            kh_value(pSampleHash, khIter) = pSampleInfo->size;
            pSampleInfo->pSamples[pSampleInfo->size] = buff;
            ++(pSampleInfo->size);
        }

        return TGM_OK;
    }

    return TGM_ERR;
}

static TGM_Status TGM_LibInfoTableAddSeqTech(TGM_LibInfoTable* pTable, const char* platformPos)
{
    char buff[30];

    const char* platformEnd = strpbrk(platformPos, " \t\n");
    unsigned int platformLen = platformEnd - platformPos;
    if (platformLen > 29)
        return TGM_ERR;

    strncpy(buff, platformPos, platformLen);
    buff[platformLen] = '\0';

    StrToUpper(buff);
    unsigned int lastIndex = pTable->size;
    if (strstr(buff, "LONG") != NULL)
        pTable->pSeqTech[lastIndex] = ST_ILLUMINA_LONG;
    else if (strcmp(buff, "ILLUMINA") == 0)
        pTable->pSeqTech[lastIndex] = ST_ILLUMINA;
    else if (strcmp(buff, "LS454") == 0)
        pTable->pSeqTech[lastIndex] = ST_454;
    else if (strcmp(buff, "SOLID") == 0)
        pTable->pSeqTech[lastIndex] = ST_SOLID;
    else
        return TGM_ERR;

    return TGM_OK;
}

static TGM_Status TGM_LibInfoTableAddReadGrp(TGM_LibInfoTable* pTable, const char* tagPos, const char* lineEnd, int sampleID)
{
    // get the name of the current read group
    const char* readGrpNamePos = strstr(tagPos, "ID:");
    if (readGrpNamePos != NULL && readGrpNamePos < lineEnd)
        readGrpNamePos += 3;
    else
        return TGM_ERR;

    const char* platformPos = strstr(tagPos, "PL:");
    if (platformPos != NULL && platformPos < lineEnd)
        platformPos += 3;
    else
        return TGM_ERR;

    // expand the array if necessary
    if (pTable->size == pTable->capacity)
    {
        pTable->capacity *= 2;
        pTable->pReadGrps = (char**) realloc(pTable->pReadGrps, pTable->capacity * sizeof(char*));
        if (pTable->pReadGrps == NULL)
            TGM_ErrQuit("ERROR: Not enough memory for the storage of read group names in the library information object.\n");

        pTable->pSampleMap = (int32_t*) realloc(pTable->pSampleMap, pTable->capacity * sizeof(int32_t));
        if (pTable->pSampleMap == NULL)
            TGM_ErrQuit("ERROR: Not enough memory for the storage of the sample ID map in the library information object.\n");

        pTable->pSeqTech = (int8_t*) realloc(pTable->pSeqTech, pTable->capacity * sizeof(int8_t));
        if (pTable->pSeqTech == NULL)
            TGM_ErrQuit("ERROR: Not enough memory for the storage of the sequencing technology in the library information  object.\n");
    }

    TGM_Status status = TGM_LibInfoTableAddSeqTech(pTable, platformPos);
    if (status != TGM_OK)
        return TGM_ERR;

    const char* readGrpNameEnd = strpbrk(readGrpNamePos, " \t\n\0");
    size_t readGrpNameLen = readGrpNameEnd - readGrpNamePos;
    if (readGrpNameLen > 0)
    {
        pTable->pReadGrps[pTable->size] = (char*) calloc(readGrpNameLen + 1, sizeof(char));
        if (pTable->pReadGrps[pTable->size] == NULL)
            TGM_ErrQuit("ERROR: Not enough memory for the storage of read group name in the fragment length distribution object.\n");

        memcpy(pTable->pReadGrps[pTable->size], readGrpNamePos, readGrpNameLen);

        int ret = 0;

        khash_t(name)* pReadGrpHash = pTable->pReadGrpHash;
        khiter_t khIter = kh_put(name, pReadGrpHash, pTable->pReadGrps[pTable->size], &ret);

        if (ret != 0)
        {
            pTable->pSampleMap[pTable->size] = sampleID;
            kh_value(pReadGrpHash, khIter) = pTable->size;
            ++(pTable->size);
        }
        else
        {
            free(pTable->pReadGrps[pTable->size]);
            pTable->pReadGrps[pTable->size] = NULL;
        }

        return TGM_OK;
    }

    return TGM_ERR;
}

static void TGM_AnchorInfoWrite(const TGM_AnchorInfo* pInfo, FILE* libFile)
{
    unsigned int writeSize = 0;

    writeSize = fwrite(pInfo->pLength, sizeof(int32_t), pInfo->size, libFile);
    if (writeSize != pInfo->size)
        TGM_ErrQuit("ERROR: Cannot write the length of anchors into the information file.\n");

    writeSize = fwrite(pInfo->pMd5s, sizeof(char), pInfo->size * MD5_STR_LEN, libFile);
    if (writeSize != pInfo->size * MD5_STR_LEN)
        TGM_ErrQuit("ERROR: Cannot write the md5 string into the information file.\n");

    for (unsigned int i = 0; i != pInfo->size; ++i)
    {
        uint32_t anchorNameLen = strlen(pInfo->pAnchors[i]);
        writeSize = fwrite(&(anchorNameLen), sizeof(uint32_t), 1, libFile);
        if (writeSize != 1)
            TGM_ErrQuit("ERROR: Cannot write the anchor name length into the information file.\n");

        writeSize = fwrite(pInfo->pAnchors[i], sizeof(char), anchorNameLen, libFile);
        if (writeSize != anchorNameLen)
            TGM_ErrQuit("ERROR: Cannot write the sample name into the information file.\n");
    }

    writeSize = fwrite(&(pInfo->specialPrefixLen), sizeof(uint32_t), 1, libFile);
    if (writeSize != 1)
        TGM_ErrQuit("ERROR: Cannot write the length of the special prefix into the information file.\n");

    if (pInfo->specialPrefixLen > 0)
    {
        writeSize = fwrite(pInfo->pSpecialPrefix, sizeof(char), pInfo->specialPrefixLen, libFile);
        if (writeSize != pInfo->specialPrefixLen)
            TGM_ErrQuit("ERROR: Cannot write the special prefix into the information file.\n");
    }
}

static void TGM_SampleInfoWrite(const TGM_SampleInfo* pInfo, FILE* libFile)
{
    unsigned int writeSize = 0;

    for (unsigned int i = 0; i != pInfo->size; ++i)
    {
        uint32_t sampleNameLen = strlen(pInfo->pSamples[i]);
        writeSize = fwrite(&(sampleNameLen), sizeof(uint32_t), 1, libFile);
        if (writeSize != 1)
            TGM_ErrQuit("ERROR: Cannot write the sample name length into the information file.\n");

        writeSize = fwrite(pInfo->pSamples[i], sizeof(char), sampleNameLen, libFile);
        if (writeSize != sampleNameLen)
            TGM_ErrQuit("ERROR: Cannot write the sample name into the information file.\n");
    }
}

static void TGM_AnchorInfoRead(TGM_AnchorInfo* pInfo, FILE* libFile)
{
    unsigned int readSize = 0;
    
    readSize = fread(pInfo->pLength, sizeof(int32_t), pInfo->size, libFile);
    if (readSize != pInfo->size)
        TGM_ErrQuit("ERROR: Cannot read the length of anchors from the library file.\n");

    readSize = fread(pInfo->pMd5s, sizeof(char), pInfo->size * MD5_STR_LEN, libFile);
    if (readSize != pInfo->size * MD5_STR_LEN)
        TGM_ErrQuit("ERROR: Cannot read the md5 strings from the library file.\n");

    uint32_t anchorNameLen = 0;
    for (unsigned int i = 0; i != pInfo->size; ++i)
    {
        readSize = fread(&anchorNameLen, sizeof(uint32_t), 1, libFile);
        if (readSize != 1)
            TGM_ErrQuit("ERROR: Cannot read the length of anchor name from the library file.\n");

        pInfo->pAnchors[i] = (char*) malloc((anchorNameLen + 1) * sizeof(char));
        if (pInfo->pAnchors[i] == NULL)
            TGM_ErrQuit("ERROR: Not enough memory for the storage of the anchor name.\n");

        pInfo->pAnchors[i][anchorNameLen] = '\0';
        
        readSize = fread(pInfo->pAnchors[i], sizeof(char), anchorNameLen, libFile);
        if (readSize != anchorNameLen)
            TGM_ErrQuit("ERROR: Cannot read the anchor name from the library file.\n");

        int ret = 0;
        khiter_t khIter = kh_put(name, pInfo->pAnchorHash, pInfo->pAnchors[i], &ret);
        kh_value((khash_t(name)*) pInfo->pAnchorHash, khIter) = i;
    }

    readSize = fread(&(pInfo->specialPrefixLen), sizeof(uint32_t), 1, libFile);
    if (readSize != 1)
        TGM_ErrQuit("ERROR: Cannot read the length of the special prefix from the library file.\n");

    if (pInfo->specialPrefixLen > 0)
    {
        if (pInfo->pSpecialPrefix == NULL)
        {
            pInfo->pSpecialPrefix = (char*) malloc(sizeof(char) * (pInfo->specialPrefixLen + 1));
            if (pInfo->pSpecialPrefix == NULL)
                TGM_ErrQuit("ERROR: Not enough memory for the special reference prefix.\n");

            pInfo->pSpecialPrefix[pInfo->specialPrefixLen] = '\0';
        }

        readSize = fread(pInfo->pSpecialPrefix, sizeof(char), pInfo->specialPrefixLen, libFile);
        if (readSize != pInfo->specialPrefixLen)
            TGM_ErrQuit("ERROR: Cannot read the special prefix from the library file.\n");
    }
}

static void TGM_SampleInfoRead(TGM_SampleInfo* pInfo, FILE* libFile)
{
    unsigned int readSize = 0;
    uint32_t sampleNameLen = 0;

    for (unsigned int i = 0; i != pInfo->size; ++i)
    {
        readSize = fread(&sampleNameLen, sizeof(uint32_t), 1, libFile);
        if (readSize != 1)
            TGM_ErrQuit("ERROR: Cannot read the length of sample name from the library file.\n");

        pInfo->pSamples[i] = (char*) malloc((sampleNameLen + 1) * sizeof(char));
        if (pInfo->pSamples[i] == NULL)
            TGM_ErrQuit("ERROR: Not enough memory for the storage of the sample name.\n");

        pInfo->pSamples[i][sampleNameLen] = '\0';
        
        readSize = fread(pInfo->pSamples[i], sizeof(char), sampleNameLen, libFile);
        if (readSize != sampleNameLen)
            TGM_ErrQuit("ERROR: Cannot read the sample name from the library file.\n");

        int ret = 0;
        khiter_t khIter = kh_put(name, pInfo->pSampleHash, pInfo->pSamples[i], &ret);
        kh_value((khash_t(name)*) pInfo->pSampleHash, khIter) = i;
    }
}


TGM_AnchorInfo* TGM_AnchorInfoAlloc(uint32_t capacity, const char* pSpecialPrefix)
{
    TGM_AnchorInfo* pNewInfo = (TGM_AnchorInfo*) malloc(sizeof(TGM_AnchorInfo));
    if (pNewInfo == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for an anchor information object.\n");

    pNewInfo->pAnchors = (char**) malloc(capacity * sizeof(char*));
    if (pNewInfo->pAnchors == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for the storage of anchor names in the anchor information object.\n");

    pNewInfo->pLength = (int32_t*) malloc(capacity * sizeof(int32_t));
    if (pNewInfo->pLength == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for the storage of anchor length in the anchor information object.\n");

    pNewInfo->pMd5s = (char*) malloc(capacity * MD5_STR_LEN * sizeof(char));
    if (pNewInfo->pMd5s == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for the storage of md5 string in the anchor information object.\n");

    pNewInfo->size = 0;
    pNewInfo->capacity = capacity;

    if (pSpecialPrefix != NULL)
    {
        pNewInfo->specialPrefixLen = strlen(pSpecialPrefix);
        pNewInfo->pSpecialPrefix = (char*) malloc(sizeof(char) * (pNewInfo->specialPrefixLen + 1));
        pNewInfo->pSpecialPrefix[pNewInfo->specialPrefixLen] = '\0';

        strcpy(pNewInfo->pSpecialPrefix, pSpecialPrefix);
    }
    else
    {
        pNewInfo->pSpecialPrefix = NULL;
        pNewInfo->specialPrefixLen = 0;
    }

    pNewInfo->pAnchorHash = kh_init(name);
    kh_resize(name, pNewInfo->pAnchorHash, 2 * capacity);

    return pNewInfo;
}

void TGM_AnchorInfoFree(TGM_AnchorInfo* pInfo)
{
    if (pInfo != NULL)
    {
        for (unsigned int i = 0; i != pInfo->size; ++i)
            free(pInfo->pAnchors[i]);

        kh_destroy(name, pInfo->pAnchorHash);

        free(pInfo->pAnchors);
        free(pInfo->pLength);
        free(pInfo->pMd5s);
        free(pInfo->pSpecialPrefix);
        free(pInfo);
    }
}

TGM_SampleInfo* TGM_SampleInfoAlloc(uint32_t capacity)
{
    TGM_SampleInfo* pNewInfo = (TGM_SampleInfo*) malloc(sizeof(TGM_SampleInfo));
    if (pNewInfo == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for a sample information object.\n");

    pNewInfo->pSamples = (char**) malloc(capacity * sizeof(char*));
    if (pNewInfo->pSamples == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for the storage of sample names in the sample information object.\n");

    pNewInfo->pReadFraction = (double*) malloc(capacity * sizeof(double));
    if (pNewInfo->pReadFraction == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for the storage of read fraction in the sample information object.\n");

    pNewInfo->pSampleHash = kh_init(name);
    kh_resize(name, pNewInfo->pSampleHash, 2 * capacity);

    pNewInfo->size = 0;
    pNewInfo->capacity = capacity;

    return pNewInfo;
}

void TGM_SampleInfoFree(TGM_SampleInfo* pInfo)
{
    if (pInfo != NULL)
    {
        if (pInfo->pSamples != NULL)
        {
            for (unsigned int i = 0; i != pInfo->size; ++i)
                free(pInfo->pSamples[i]);

            free(pInfo->pSamples);
        }

        kh_destroy(name, pInfo->pSampleHash);
        free(pInfo->pReadFraction);
        free(pInfo);
    }
}

TGM_SpecialID* TGM_SpecialIDAlloc(unsigned int capacity)
{
    TGM_SpecialID* pSpecialID = (TGM_SpecialID*) malloc(sizeof(TGM_SpecialID));
    if (pSpecialID == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for the special ID object.\n");

    pSpecialID->names = (char (*)[3]) malloc(capacity * 3 * sizeof(char));
    if (pSpecialID->names == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for the special ID names.\n");

    pSpecialID->pHash = kh_init(name);
    pSpecialID->size = 0;
    pSpecialID->capacity = capacity;

    return pSpecialID;
}

void TGM_SpecialIDFree(TGM_SpecialID* pSpecialID)
{
    if (pSpecialID != NULL)
    {
        free(pSpecialID->names);
        kh_destroy(name, pSpecialID->pHash);

        free(pSpecialID);
    }
}


TGM_LibInfoTable* TGM_LibInfoTableAlloc(uint32_t capAnchor, uint32_t capSample, uint32_t capReadGrp, const char* pSpecialPrefix)
{
    TGM_LibInfoTable* pNewTable = (TGM_LibInfoTable*) malloc(sizeof(TGM_LibInfoTable));
    if (pNewTable == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for a library information table object.\n");

    pNewTable->pSampleInfo = TGM_SampleInfoAlloc(capSample);

    pNewTable->pLibInfo = (TGM_LibInfo*) malloc(capReadGrp * sizeof(TGM_LibInfo));
    if (pNewTable->pLibInfo == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for the storage of library information in an library table object.\n");

    pNewTable->pReadGrps = (char**) malloc(capReadGrp * sizeof(char*));
    if (pNewTable->pLibInfo == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for the storage of read group names in an library table object.\n");

    pNewTable->pSampleMap = (int32_t*) malloc(capReadGrp * sizeof(int32_t));
    if (pNewTable->pSampleMap == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for the storage of read-group-to-sample map in an library table object.\n");

    pNewTable->pSeqTech = (int8_t*) malloc(capReadGrp * sizeof(int8_t));
    if (pNewTable->pSeqTech == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for the storage of sequencing technologies in an library table object.\n");

    pNewTable->pAnchorInfo = TGM_AnchorInfoAlloc(capAnchor, pSpecialPrefix);

    pNewTable->pReadGrpHash = kh_init(name);
    kh_resize(name, pNewTable->pReadGrpHash, 2 * capReadGrp);

    pNewTable->size = 0;
    pNewTable->capacity = capReadGrp;
    pNewTable->fragLenMax = 0;
    pNewTable->cutoff = 0.0;
    pNewTable->trimRate = 0.0;

    return pNewTable;
}

void TGM_LibInfoTableFree(TGM_LibInfoTable* pTable)
{
    if (pTable != NULL)
    {
        if (pTable->pReadGrps != NULL)
        {
            for (unsigned int i = 0; i != pTable->size; ++i)
                free(pTable->pReadGrps[i]);

            free(pTable->pReadGrps);
        }

        TGM_SampleInfoFree(pTable->pSampleInfo);
        TGM_AnchorInfoFree(pTable->pAnchorInfo);

        kh_destroy(name, pTable->pReadGrpHash);

        free(pTable->pLibInfo);
        free(pTable->pSampleMap);
        free(pTable->pSeqTech);

        free(pTable);
    }
}

TGM_Status TGM_LoadZAtag(TGM_ZAtag* pZAtag, const bam1_t* pUpAlgn)
{
    uint8_t* zaPos = bam_aux_get(pUpAlgn, "ZA");
    const char* zaStr = NULL;
    if (zaPos != NULL)
        zaStr = bam_aux2Z(zaPos);
    else
        return TGM_NOT_FOUND;

    // set the first character of the special reference name to be ' '
    // we can later check if the pair hit any special references
    pZAtag->spRef[0][0] = ' ';
    pZAtag->spRef[1][0] = ' ';

    pZAtag->spRef[0][2] = '\0';
    pZAtag->spRef[1][2] = '\0';

    unsigned int whichMate = 0;
    int numMappings = 0;
    const char* currFieldPos = zaStr + 1;
    if (*currFieldPos == '&')
        whichMate = 1;

    currFieldPos += 2;
    unsigned int currFieldNum = 1;
    const char* currFieldEnd = NULL;

    const char* cigarStr = NULL;
    unsigned int cigarLen = 0;

    const char* mdStr = NULL;
    unsigned int mdLen = 0;

    do
    {
        ++currFieldNum;
        currFieldEnd = strpbrk(currFieldPos, ";>");

        if (currFieldEnd - currFieldPos == 0)
        {
            ++currFieldPos;
            continue;
        }

        switch (currFieldNum)
        {
            case 2:
                pZAtag->bestMQ[whichMate] = atoi(currFieldPos);
                break;
            case 3:
                pZAtag->secMQ[whichMate] = atoi(currFieldPos);
                break;
            case 4:
                pZAtag->spRef[whichMate][0] = *currFieldPos;
                pZAtag->spRef[whichMate][1] = *(currFieldPos + 1);
                break;
            case 5:
                numMappings = atoi(currFieldPos);
                pZAtag->numMappings[whichMate] = numMappings;
                if (numMappings > UINT16_MAX)
                    pZAtag->numMappings[whichMate] = UINT16_MAX;
                break;
            case 6:
                cigarStr = currFieldPos;
                cigarLen = currFieldEnd - currFieldPos;
                break;
            case 7:
                mdStr = currFieldPos;
                mdLen = currFieldEnd - currFieldPos;
                break;
            default:
                break;
        }

        currFieldPos = currFieldEnd + 1;

    }while(*currFieldEnd != '>');

    if (cigarStr != NULL && mdStr != NULL)
    {
        TGM_GetNumMismatchFromZA(&(pZAtag->numMM[whichMate]), &(pZAtag->end[whichMate]), cigarStr, cigarLen, mdStr, mdLen);
        pZAtag->end[whichMate] = pZAtag->end[whichMate] + pUpAlgn->core.mpos - 1;
    }
    else
    {
        pZAtag->numMM[whichMate] = TGM_GetNumMismatchFromBam(pUpAlgn);
        pZAtag->end[whichMate] = bam_calend(&(pUpAlgn->core), bam1_cigar(pUpAlgn));
    }


    // handle the sencond ZA tag
    currFieldPos += 3;
    currFieldNum = 1;
    currFieldEnd = NULL;

    cigarStr = NULL;
    cigarLen = 0;

    mdStr = NULL;
    mdLen = 0;

    whichMate ^= 1;

    do
    {
        ++currFieldNum;
        currFieldEnd = strpbrk(currFieldPos, ";>");

        if (currFieldEnd - currFieldPos == 0)
        {
            ++currFieldPos;
            continue;
        }

        switch (currFieldNum)
        {
            case 2:
                pZAtag->bestMQ[whichMate] = atoi(currFieldPos);
                break;
            case 3:
                pZAtag->secMQ[whichMate] = atoi(currFieldPos);
                break;
            case 4:
                pZAtag->spRef[whichMate][0] = *currFieldPos;
                pZAtag->spRef[whichMate][1] = *(currFieldPos + 1);
                break;
            case 5:
                numMappings = atoi(currFieldPos);
                pZAtag->numMappings[whichMate] = numMappings;
                if (numMappings > UINT16_MAX)
                    pZAtag->numMappings[whichMate] = UINT16_MAX;
                break;
            case 6:
                cigarStr = currFieldPos;
                cigarLen = currFieldEnd - currFieldPos;
                break;
            case 7:
                mdStr = currFieldPos;
                mdLen = currFieldEnd - currFieldPos;
            default:
                break;
        }

        currFieldPos = currFieldEnd + 1;

    }while(*currFieldEnd != '>');

    if (cigarStr != NULL && mdStr != NULL)
    {
        TGM_GetNumMismatchFromZA(&(pZAtag->numMM[whichMate]), &(pZAtag->end[whichMate]), cigarStr, cigarLen, mdStr, mdLen);
        pZAtag->end[whichMate] = pZAtag->end[whichMate] + pUpAlgn->core.mpos - 1;
    }
    else
    {
        pZAtag->numMM[whichMate] = TGM_GetNumMismatchFromBam(pUpAlgn);
        pZAtag->end[whichMate] = bam_calend(&(pUpAlgn->core), bam1_cigar(pUpAlgn));
    }

    return TGM_OK;
}

TGM_Status TGM_LoadPairStats(TGM_PairStats* pPairStats, const bam1_t* pAlignment, const TGM_LibInfoTable* pTable)
{
    pPairStats->pairMode = TGM_GetPairMode(pAlignment);
    if (pPairStats->pairMode == TGM_BAD_PAIR_MODE)
        return TGM_ERR;

    pPairStats->fragLen = abs(pAlignment->core.isize);

    if (pAlignment->core.tid != pAlignment->core.mtid)
        pPairStats->fragLen = -1;

    static const char tagRG[2] = {'R', 'G'};
    uint8_t* rgPos = bam_aux_get(pAlignment, tagRG);
    if (rgPos != NULL)
    {
        const char* RG = bam_aux2Z(rgPos);
        TGM_Status status = TGM_LibInfoTableGetRGIndex(&(pPairStats->readGrpID), pTable, RG);
        if (status != TGM_OK)
            return TGM_ERR;

        TGM_SeqTech seqTech = pTable->pSeqTech[pPairStats->readGrpID];
        pPairStats->pairMode = TGM_SeqTechMap[seqTech][pPairStats->pairMode];
    }
    else
        return TGM_ERR;

    return TGM_OK;
}

TGM_Bool TGM_IsNormalPair(TGM_PairStats* pPairStats, TGM_ZAtag* pZAtag, TGM_Status* pZAstatus, const TGM_MateInfo* pMateInfo, unsigned int* pBackHistIndex, 
                        const bam1_t* pAlgns[3], int retNum, const TGM_LibInfoTable* pTable, unsigned short minMQ)
{
    const bam1_t* pUpAlgn = NULL;
    const bam1_t* pDownAlgn = NULL;

    *pZAstatus = TGM_ERR;

    // check the mapping quality of the read pair 
    // to make sure that both mates are unique
    if (retNum == 1)
    {
        // for sorted by coordinate no za, this is acutally the downstream mate
        // but it does not matter for this funcion
        pUpAlgn = pAlgns[0];
        if (pMateInfo != NULL)
        {
            if (pMateInfo->mapQ < minMQ || pUpAlgn->core.qual < minMQ)
                return FALSE;
        }
        else
        {
            TGM_Status zaStatus = TGM_LoadZAtag(pZAtag, pUpAlgn);
            if (zaStatus != TGM_OK)
                return FALSE;

            *pZAstatus = TGM_OK;

            if (pZAtag->bestMQ[0] < minMQ || pZAtag->bestMQ[1] < minMQ)
                return FALSE;
        }
    }
    else if (retNum == 2)
    {
        pUpAlgn = pAlgns[0];
        pDownAlgn = pAlgns[1];

        if (pUpAlgn->core.qual < minMQ || pDownAlgn->core.qual < minMQ)
            return FALSE;
    }
    else
        return FALSE;

    TGM_Status status = TGM_LoadPairStats(pPairStats, pUpAlgn, pTable);
    if (status != TGM_OK)
        return FALSE;

    *pBackHistIndex = pTable->size - pPairStats->readGrpID;

    // check the orientation of the reads
    if (SV_ReadPairTypeMap[0][pPairStats->pairMode] != SV_NORMAL && SV_ReadPairTypeMap[1][pPairStats->pairMode] != SV_NORMAL)
        return FALSE;

    return TRUE;
}

TGM_Status TGM_LibInfoTableSetRG(TGM_LibInfoTable* pTable, unsigned int* oldSize, const TGM_BamHeader* pBamHeader)
{
    TGM_Status status = TGM_OK;
    *oldSize = pTable->size;
    unsigned int oldCapacity = pTable->capacity;

    const char* tagPos = pBamHeader->pOrigHeader->text;
    status = TGM_LibInfoTableAddAnchor(pTable, pBamHeader);
    if (status != TGM_OK)
        return TGM_ERR;

    while ((tagPos = strstr(tagPos, "@RG")) != NULL)
    {
        const char* lineEnd = strpbrk(tagPos, "\n");
        int32_t sampleID = 0;

        status = TGM_LibInfoTableAddSample(&sampleID, pTable, tagPos, lineEnd);
        if (status != TGM_OK)
        {
            TGM_ErrMsg("ERROR: the \"SM\" field is not found under the read group tag in the bam header.\n");

            ++tagPos;
            continue;
        }

        status = TGM_LibInfoTableAddReadGrp(pTable, tagPos, lineEnd, sampleID);
        if (status != TGM_OK)
            TGM_ErrMsg("ERROR: the \"ID\" field is required under the read group tag.\n");

        ++tagPos;
    }

    if (pTable->capacity > oldCapacity)
    {
        pTable->pLibInfo = (TGM_LibInfo*) realloc(pTable->pLibInfo, sizeof(TGM_LibInfo) * pTable->capacity);
        if (pTable->pLibInfo == NULL)
            TGM_ErrQuit("ERROR: Not enough memory for the storage of library summary in the library information object.\n");
    }

    return TGM_OK;
}

TGM_Status TGM_LibInfoTableSetRGSplit(TGM_LibInfoTable* pTable, const TGM_BamHeader* pBamHeader, const char* pSpecialPrefix, unsigned int prefixLen)
{
    const char* tagPos = pBamHeader->pOrigHeader->text;

    TGM_Status status = TGM_OK;
    status = TGM_LibInfoTableAddAnchorSplit(pTable, pBamHeader, pSpecialPrefix, prefixLen);
    if (status != TGM_OK)
        return TGM_ERR;

    while ((tagPos = strstr(tagPos, "@RG")) != NULL)
    {
        const char* lineEnd = strpbrk(tagPos, "\n");
        int32_t sampleID = 0;

        status = TGM_LibInfoTableAddSample(&sampleID, pTable, tagPos, lineEnd);
        if (status != TGM_OK)
        {
            TGM_ErrMsg("ERROR: the \"SM\" field is not found under the read group tag in the bam header.\n");

            ++tagPos;
            continue;
        }

        status = TGM_LibInfoTableAddReadGrp(pTable, tagPos, lineEnd, sampleID);
        if (status != TGM_OK)
            TGM_ErrMsg("ERROR: the \"ID\" field is required under the read group tag.\n");

        ++tagPos;
    }

    return TGM_OK;
}

TGM_Status TGM_LibInfoTableGetRGIndex(int32_t* pReadGrpIndex, const TGM_LibInfoTable* pTable, const char* pReadGrpName)
{
    *pReadGrpIndex = 0;

    khash_t(name)* pRgHash = pTable->pReadGrpHash;
    khiter_t khIter = kh_get(name, pRgHash, pReadGrpName);
    if (khIter != kh_end(pRgHash))
    {
        *pReadGrpIndex = kh_value(pRgHash, khIter);
    }
    else
    {
        TGM_ErrMsg("ERROR: Found a read group name that is not recorded in the library information table.\n");
        return TGM_ERR;
    }

    return TGM_OK;
}


void TGM_LibInfoTableUpdate(TGM_LibInfoTable* pTable, const TGM_FragLenHistArray* pHistArray, unsigned int oldSize)
{
    for (unsigned int i = 0; i != pHistArray->size; ++i)
    {
        const TGM_FragLenHist* pHist = pHistArray->data + i;
        if (pHist->size < MIN_FRAGLEN_HIST_SIZE || pHist->modeCount[0] < MIN_FRAGLEN_HIST_FREQ)
        {
            pTable->pLibInfo[oldSize].fragLenMedian = 0;
            pTable->pLibInfo[oldSize].fragLenLow = 0;
            pTable->pLibInfo[oldSize].fragLenHigh = 0;

            ++oldSize;
            continue;
        }

        pTable->pLibInfo[oldSize].fragLenMedian = pHist->median;

        double oneSideFlow = pTable->trimRate / 2.0 * pHist->modeCount[0];
        double total = pHist->modeCount[0] * (1.0 - pTable->trimRate);

        double cumFreq = 0.0;
        double oneSideCutoff = pTable->cutoff / 2;
        for (unsigned int i = 0; i != pHist->size; ++i)
        {
            cumFreq += pHist->freq[i];
            if (((cumFreq - oneSideFlow) / total) > oneSideCutoff)
            {
                pTable->pLibInfo[oldSize].fragLenLow = pHist->fragLen[i];
                break;
            }
        }

        cumFreq = 0.0;
        for (int i = pHist->size - 1; i != -1; --i)
        {
            cumFreq += pHist->freq[i];
            if (((cumFreq - oneSideFlow) / total) > oneSideCutoff)
            {
                pTable->pLibInfo[oldSize].fragLenHigh = pHist->fragLen[i];
                break;
            }
        }

        if (pTable->pLibInfo[oldSize].fragLenHigh > pTable->fragLenMax)
            pTable->fragLenMax = pTable->pLibInfo[oldSize].fragLenHigh;

        ++oldSize;
    }
}

void TGM_LibInfoTableWrite(const TGM_LibInfoTable* pTable, TGM_Bool writeInfo, FILE* libFile)
{
    unsigned int writeSize = 0;

    writeSize = fwrite(&(pTable->pAnchorInfo->size), sizeof(uint32_t), 1, libFile);
    if (writeSize != 1)
        TGM_ErrQuit("ERROR: Cannot write the number of anchors into the information file.\n");

    writeSize = fwrite(&(pTable->pSampleInfo->size), sizeof(uint32_t), 1, libFile);
    if (writeSize != 1)
        TGM_ErrQuit("ERROR: Cannot write the number of samples into the information file.\n");

    writeSize = fwrite(&(pTable->size), sizeof(uint32_t), 1, libFile);
    if (writeSize != 1)
        TGM_ErrQuit("ERROR: Cannot write the number of read groups into the information file.\n");

    TGM_AnchorInfoWrite(pTable->pAnchorInfo, libFile);
    TGM_SampleInfoWrite(pTable->pSampleInfo, libFile);

    writeSize = fwrite(pTable->pSampleMap, sizeof(int32_t), pTable->size, libFile);
    if (writeSize != pTable->size)
        TGM_ErrQuit("ERROR: Cannot write the read-group-to-sample map into the information file.\n");

    for (unsigned int i = 0; i != pTable->size; ++i)
    {
        uint32_t readGrpNameLen = strlen(pTable->pReadGrps[i]);
        writeSize = fwrite(&(readGrpNameLen), sizeof(uint32_t), 1, libFile);
        if (writeSize != 1)
            TGM_ErrQuit("ERROR: Cannot write the readGrp name length into the information file.\n");

        writeSize = fwrite(pTable->pReadGrps[i], sizeof(char), readGrpNameLen, libFile);
        if (writeSize != readGrpNameLen)
            TGM_ErrQuit("ERROR: Cannot write the readGrp name into the information file.\n");
    }

    if (writeInfo)
    {
        writeSize = fwrite(&(pTable->fragLenMax), sizeof(uint32_t), 1, libFile);
        if (writeSize != 1)
            TGM_ErrQuit("ERROR: Cannot write the maximum fragment length into the library file.\n");

        writeSize = fwrite(&(pTable->cutoff), sizeof(double), 1, libFile);
        if (writeSize != 1)
            TGM_ErrQuit("ERROR: Cannot write the cutoff into the library information file.\n");

        writeSize = fwrite(&(pTable->trimRate), sizeof(double), 1, libFile);
        if (writeSize != 1)
            TGM_ErrQuit("ERROR: Cannot write the trim rate into the library information file.\n");

        writeSize = fwrite(pTable->pSeqTech, sizeof(int8_t), pTable->size, libFile);
        if (writeSize != pTable->size)
            TGM_ErrQuit("ERROR: Cannot write the sequencing technology from the library information file.\n");

        writeSize = fwrite(pTable->pLibInfo, sizeof(TGM_LibInfo), pTable->size, libFile);
        if (writeSize != pTable->size)
            TGM_ErrQuit("ERROR: Cannot write the library information into the output file.\n");
    }
    else
    {
        uint32_t zero = 0;
        writeSize = fwrite(&zero, sizeof(uint32_t), 1, libFile);
        if (writeSize != 1)
            TGM_ErrQuit("ERROR: Cannot write the maximum fragment length into the library file.\n");
    }

    fflush(libFile);
}

void TGM_SpecialIDWrite(const TGM_SpecialID* pSpecialID, FILE* pLibOutput)
{
    unsigned int writeSize = 0;

    writeSize = fwrite(&(pSpecialID->size), sizeof(uint32_t), 1, pLibOutput);
    if (writeSize != 1)
        TGM_ErrQuit("ERROR: Cannot write the number of special reference ID into the library information file.\n");

    for (unsigned int i = 0; i != pSpecialID->size; ++i)
    {
        writeSize = fwrite(pSpecialID->names[i], sizeof(char), 2, pLibOutput);
        if (writeSize != 2)
            TGM_ErrQuit("ERROR: Cannot write the special reference ID into the library information file.\n");
    }
}

void TGM_SpecialIDRead(TGM_SpecialID* pSpecialID, FILE* pLibInput)
{
    unsigned int readSize = fread(&(pSpecialID->size), sizeof(uint32_t), 1, pLibInput);
    if (readSize != 1)
        TGM_ErrQuit("ERROR: Cannot read the number of special reference ID from the library information file.\n");

    if (pSpecialID->size > pSpecialID->capacity)
    {
        pSpecialID->capacity = pSpecialID->size * 2;
        free(pSpecialID->names);
        pSpecialID->names = (char (*)[3]) malloc(sizeof(char) * 3 * pSpecialID->capacity);
        if (pSpecialID->names == NULL)
            TGM_ErrQuit("ERROR: Not enough memory for the special ID names.\n");
    }

    kh_clear(name, pSpecialID->pHash);

    int ret = 0;
    khiter_t khIter = 0;

    for (unsigned int i = 0; i != pSpecialID->size; ++i)
    {
        readSize = fread(pSpecialID->names[i], sizeof(char), 2, pLibInput);
        if (readSize != 2)
            TGM_ErrQuit("ERROR: Cannot read the special reference ID into the library information file.\n");

        pSpecialID->names[i][2] = '\0';

        khIter = kh_put(name, pSpecialID->pHash, pSpecialID->names[i], &ret);
        kh_value((khash_t(name)*) pSpecialID->pHash, khIter) = i;
    }
}

void TGM_SpecialIDMergeRead(TGM_SpecialID* pSpecialID, FILE* pLibOutput)
{
    uint32_t numSpecialID = 0;
    unsigned int readSize = fread(&numSpecialID, sizeof(uint32_t), 1, pLibOutput);
    if (readSize != 1)
        TGM_ErrQuit("ERROR: Cannot read the number of special reference ID from the library information file.\n");

    if (numSpecialID > pSpecialID->capacity)
    {
        pSpecialID->capacity = numSpecialID * 2;
        pSpecialID->names = (char (*)[3]) realloc(pSpecialID->names, sizeof(char) * 3 * pSpecialID->capacity);
        if (pSpecialID->names == NULL)
            TGM_ErrQuit("ERROR: Not enough memory for the special ID names.\n");
    }

    int ret = 0;
    khiter_t khIter = 0;

    for (unsigned int i = 0; i != numSpecialID; ++i)
    {
        readSize = fread(pSpecialID->names[pSpecialID->size], sizeof(char), 2, pLibOutput);
        if (readSize != 2)
            TGM_ErrQuit("ERROR: Cannot read the special reference ID into the library information file.\n");

        pSpecialID->names[pSpecialID->size][2] = '\0';

        khIter = kh_put(name, pSpecialID->pHash, pSpecialID->names[pSpecialID->size], &ret);
        if (ret != 0)
        {
            kh_value((khash_t(name)*) pSpecialID->pHash, khIter) = pSpecialID->size;
            ++(pSpecialID->size);
        }
    }
}

TGM_LibInfoTable* TGM_LibInfoTableRead(FILE* libFile)
{
    unsigned int readSize = 0;
    uint32_t sizeAC = 0;
    uint32_t sizeSM = 0;
    uint32_t sizeRG = 0;

    readSize = fread(&sizeAC, sizeof(uint32_t), 1, libFile);
    if (readSize != 1)
        TGM_ErrQuit("ERROR: Cannot read the number of anchors from the library file.\n");

    readSize = fread(&sizeSM, sizeof(uint32_t), 1, libFile);
    if (readSize != 1)
        TGM_ErrQuit("ERROR: Cannot read the number of sample from the library file.\n");

    readSize = fread(&sizeRG, sizeof(uint32_t), 1, libFile);
    if (readSize != 1)
        TGM_ErrQuit("ERROR: Cannot read the number of read group from the library file.\n");

    TGM_LibInfoTable* pTable = TGM_LibInfoTableAlloc(sizeAC, sizeSM, sizeRG, NULL);

    pTable->pAnchorInfo->size = sizeAC;
    pTable->pSampleInfo->size = sizeSM;
    pTable->size = sizeRG;

    pTable->pAnchorInfo->capacity = sizeAC;
    pTable->pSampleInfo->capacity = sizeSM;
    pTable->capacity = sizeRG;

    TGM_AnchorInfoRead(pTable->pAnchorInfo, libFile);
    TGM_SampleInfoRead(pTable->pSampleInfo, libFile);

    readSize = fread(pTable->pSampleMap, sizeof(int32_t), sizeRG, libFile);
    if (readSize != sizeRG)
        TGM_ErrQuit("ERROR: Cannot read the read-group-to-sample map from the library file.\n");

    uint32_t readGrpNameLen = 0;

    for (unsigned int i = 0; i != sizeRG; ++i)
    {
        readSize = fread(&readGrpNameLen, sizeof(uint32_t), 1, libFile);
        if (readSize != 1)
            TGM_ErrQuit("ERROR: Cannot read the length of read group name from the library file.\n");

        pTable->pReadGrps[i] = (char*) malloc((readGrpNameLen + 1) * sizeof(char));
        if (pTable->pReadGrps[i] == NULL)
            TGM_ErrQuit("ERROR: Not enough memory for the storage of the sample name.\n");

        pTable->pReadGrps[i][readGrpNameLen] = '\0';

        readSize = fread(pTable->pReadGrps[i], sizeof(char), readGrpNameLen, libFile);
        if (readSize != readGrpNameLen)
            TGM_ErrQuit("ERROR: Cannot read the read group name from the library file.\n");

        int ret = 0;
        khiter_t khIter = kh_put(name, pTable->pReadGrpHash, pTable->pReadGrps[i], &ret);
        kh_value((khash_t(name)*) pTable->pReadGrpHash, khIter) = i;
    }

    readSize = fread(&(pTable->fragLenMax), sizeof(uint32_t), 1, libFile);
    if (readSize != 1)
        TGM_ErrQuit("ERROR: Cannot read the maximum fragment length from the library file.\n");

    if (pTable->fragLenMax > 0)
    {
        readSize = fread(&(pTable->cutoff), sizeof(double), 1, libFile);
        if (readSize != 1)
            TGM_ErrQuit("ERROR: Cannot read the cutoff from the library file.\n");

        readSize = fread(&(pTable->trimRate), sizeof(double), 1, libFile);
        if (readSize != 1)
            TGM_ErrQuit("ERROR: Cannot read the trim rate from the library file.\n");

        readSize = fread(pTable->pSeqTech, sizeof(int8_t), pTable->size, libFile);
        if (readSize != pTable->size)
            TGM_ErrQuit("ERROR: Cannot read the sequencing technology from the library file.\n");

        readSize = fread(pTable->pLibInfo, sizeof(TGM_LibInfo), pTable->size, libFile);
        if (readSize != pTable->size)
            TGM_ErrQuit("ERROR: Cannot read the library information from the library file.\n");
    }


    return pTable;
}

uint32_t TGM_LibInfoTableCountNormalChr(const TGM_LibInfoTable* pLibTable)
{
    uint32_t numChr = 0;

    for (unsigned int i = 0; i != pLibTable->pAnchorInfo->size; ++i)
    {
        if (pLibTable->pAnchorInfo->pLength[i] > 0 
            && strncmp(pLibTable->pAnchorInfo->pAnchors[i], pLibTable->pAnchorInfo->pSpecialPrefix, pLibTable->pAnchorInfo->specialPrefixLen) != 0)
        {
            ++numChr;
        }
    }

    return numChr;
}

void TGM_LibInfoTableMerge(const char* workingDir)
{
    DIR* pDir = opendir(workingDir);
    if (pDir != NULL)
    {
        struct dirent* pDirRecord = NULL;
        char libFile[TGM_MAX_LINE];
        char histFile[TGM_MAX_LINE];

        unsigned int count = 0;
        TGM_LibInfoTable* pDstLibTable = NULL;
        TGM_LibInfoTable* pSrcLibTable = NULL;

        sprintf(histFile, "%s%s/", workingDir, "merged");
        TGM_Status status = TGM_CheckWorkingDir(histFile);
        if (status != TGM_OK)
            TGM_ErrQuit("ERROR: Cannot create merged directory.\n");

        sprintf(histFile, "%s%s/%s", workingDir, "merged", TGM_HistFileName);
        FILE* pHistOutput = fopen(histFile, "wb");
        if (pHistOutput == NULL)
            TGM_ErrQuit("ERROR: Cannot open the fragment length histogram file \"%s\" for writing.\n", histFile);

        TGM_FragLenHistLite* pHistLite = TGM_FragLenHistLiteAlloc(200);
        TGM_FragLenHistArrayWriteHeader(0, pHistOutput);

        TGM_SpecialID* pSpecialID = TGM_SpecialIDAlloc(DEFAULT_NUM_SPECIAL_REF);

        while ((pDirRecord = readdir(pDir)) != NULL)
        {
            if (strcmp(".", pDirRecord->d_name) != 0 && strcmp("..", pDirRecord->d_name) != 0 && strcmp("merged", pDirRecord->d_name))
            {
                sprintf(libFile, "%s%s/%s", workingDir, pDirRecord->d_name, TGM_LibTableFileName);
                sprintf(histFile, "%s%s/%s", workingDir, pDirRecord->d_name, TGM_HistFileName);
                
                FILE* pLibTableInput = fopen(libFile, "rb");
                if (pLibTableInput == NULL)
                    TGM_ErrQuit("ERROR: Cannot open the library table file \"%s\".\n", libFile);

                FILE* pHistInput = fopen(histFile, "rb");
                if (pHistInput == NULL)
                    TGM_ErrQuit("ERROR: Cannot open the fragment length histogram file \"%s\".\n", histFile);

                pSrcLibTable = TGM_LibInfoTableRead(pLibTableInput);
                if (pSrcLibTable->fragLenMax != 0)
                {
                    TGM_FragLenHistTransfer(pHistLite, pHistInput, pHistOutput);

                    if (count != 0)
                    {
                        TGM_Status mergeStatus = TGM_LibInfoTableDoMerge(pDstLibTable, pSrcLibTable);
                        if (mergeStatus != TGM_OK)
                            TGM_ErrQuit("ERROR: Cannot merge the library table file.\n");

                        TGM_LibInfoTableFree(pSrcLibTable);
                        pSrcLibTable = NULL;
                        
                        TGM_SpecialIDMergeRead(pSpecialID, pLibTableInput);
                    }
                    else
                    {
                        TGM_SWAP(pDstLibTable, pSrcLibTable, TGM_LibInfoTable*);
                        TGM_SpecialIDRead(pSpecialID, pLibTableInput);
                    }

                    ++count;
                }
                else
                {
                    TGM_LibInfoTableFree(pSrcLibTable);
                    pSrcLibTable = NULL;
                }

                fclose(pLibTableInput);
                fclose(pHistInput);
            }
        }

        TGM_FragLenHistArrayWriteHeader(pDstLibTable->size, pHistOutput);

        sprintf(libFile, "%s%s/%s", workingDir, "merged", TGM_LibTableFileName);
        FILE* pLibTableOutput = fopen(libFile, "wb");
        if (pLibTableOutput == NULL)
            TGM_ErrQuit("ERROR: Cannot open the library information table file \"%s\" for writing.\n", libFile);

        TGM_LibInfoTableWrite(pDstLibTable, TRUE, pLibTableOutput);
        TGM_SpecialIDWrite(pSpecialID, pLibTableOutput);

        TGM_FragLenHistLiteFree(pHistLite);
        TGM_LibInfoTableFree(pDstLibTable);
        TGM_SpecialIDFree(pSpecialID);

        fclose(pHistOutput);
        fclose(pLibTableOutput);
        closedir(pDir);
    }
    else
        TGM_ErrQuit("ERROR: Cannot open the working directory \"%s\"\n", workingDir);
}

TGM_Status TGM_LibInfoTableDoMerge(TGM_LibInfoTable* pDstLibTable, TGM_LibInfoTable* pSrcLibTable)
{
    TGM_Status mergeStatus = TGM_OK;

    if (pDstLibTable->cutoff != pSrcLibTable->cutoff
        || pDstLibTable->trimRate != pDstLibTable->trimRate)
    {
        return TGM_ERR;
    }

    if (pDstLibTable->fragLenMax < pSrcLibTable->fragLenMax)
        pDstLibTable->fragLenMax = pSrcLibTable->fragLenMax;

    mergeStatus = TGM_AnchorInfoMerge(pDstLibTable->pAnchorInfo, pSrcLibTable->pAnchorInfo);
    if (mergeStatus != TGM_OK)
        return mergeStatus;

    unsigned int* pIndexMap = (unsigned int*) malloc(sizeof(unsigned int) * pSrcLibTable->pSampleInfo->size);
    if (pIndexMap == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for the sample index map.\n");

    mergeStatus = TGM_SampleInfoMerge(pDstLibTable->pSampleInfo, pIndexMap, pSrcLibTable->pSampleInfo);
    if (mergeStatus != TGM_OK)
        return mergeStatus;

    if (pDstLibTable->capacity < (pDstLibTable->size + pSrcLibTable->size))
    {
        pDstLibTable->capacity = (pDstLibTable->size + pSrcLibTable->size) * 2;

        pDstLibTable->pReadGrps = (char**) realloc(pDstLibTable->pReadGrps, sizeof(char*) * pDstLibTable->capacity);
        if (pDstLibTable->pReadGrps == NULL)
            TGM_ErrQuit("ERROR: Not enough memory for the read group names.\n");

        pDstLibTable->pLibInfo = (TGM_LibInfo*) realloc(pDstLibTable->pLibInfo, sizeof(TGM_LibInfo) * pDstLibTable->capacity);
        if (pDstLibTable->pLibInfo == NULL)
            TGM_ErrQuit("ERROR: Not enough memory for the read group names.\n");

        pDstLibTable->pSampleMap = (int32_t*) realloc(pDstLibTable->pSampleMap, sizeof(int32_t) * pDstLibTable->capacity);
        if (pDstLibTable->pSampleMap == NULL)
            TGM_ErrQuit("ERROR: Not enough memory for the sample map.\n");

        pDstLibTable->pSeqTech = (int8_t*) realloc(pDstLibTable->pSeqTech, sizeof(int8_t) * pDstLibTable->capacity);
        if (pDstLibTable->pSeqTech == NULL)
            TGM_ErrQuit("ERROR: Not enough memory for the sequencing technology array.\n");
    }

    unsigned int oldSize = pDstLibTable->size;

    int ret = 0;
    khiter_t khIter = 0;

    for (unsigned int i = 0; i != pSrcLibTable->size; ++i)
    {
        khIter = kh_put(name, pDstLibTable->pReadGrpHash, pSrcLibTable->pReadGrps[i], &ret);
        if (ret != 0)
        {
            kh_value((khash_t(name)*)pDstLibTable->pReadGrpHash, khIter) = pDstLibTable->size;
            pDstLibTable->pReadGrps[pDstLibTable->size] = pSrcLibTable->pReadGrps[i];
            pSrcLibTable->pReadGrps[i] = NULL;

            pDstLibTable->pSampleMap[pDstLibTable->size] = pIndexMap[pSrcLibTable->pSampleMap[i]];

            ++(pDstLibTable->size);
        }
        else
            return TGM_ERR;
    }

    unsigned int delta = pDstLibTable->size - oldSize;
    memcpy(pDstLibTable->pLibInfo + oldSize, pSrcLibTable->pLibInfo, sizeof(TGM_LibInfo) * delta);
    memcpy(pDstLibTable->pSeqTech + oldSize, pSrcLibTable->pSeqTech, sizeof(int8_t) * delta);

    free(pIndexMap);
    return mergeStatus;
}
