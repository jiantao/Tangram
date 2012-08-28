/*
 * =====================================================================================
 *
 *       Filename:  TGM_Utilities.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  02/10/2012 06:29:54 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <errno.h>
#include <sys/stat.h>

#include "TGM_Error.h"
#include "TGM_Utilities.h"


int FindKthSmallestInt(int array[], int size, int k)
{
    register int i;
    register int j;
    register int l;
    register int m;

    register int x;

    l = 0;
    m = size - 1;

    while (l < m)
    {
        x = array[k];
        i = l;
        j = m;

        do
        {
            while (array[i] < x)
                ++i;

            while (x < array[j])
                --j;

            if (i <= j)
            {
                TGM_SWAP(array[i], array[j], int);
                ++i;
                --j;
            }

        }while (i <= j);

        if (j < k)
            l = i;

        if (k < i)
            m = j;
    }

    return array[k];
}

int FindMedianInt(int array[], int size)
{
    int median = FindKthSmallestInt(array, size, size / 2);
    if (size % 2 == 0)
    {
        int lowerMedian = 0;
        for (int i = 0; i != size / 2; ++i)
        {
            if (array[i] > lowerMedian)
                lowerMedian = array[i];
        }

        median = (int) ((lowerMedian + median) / 2.0);
    }

    return median;
}

unsigned int FindKthSmallestUint(unsigned int array[], unsigned int size, unsigned int k)
{
    register unsigned int i;
    register unsigned int j;
    register unsigned int l;
    register unsigned int m;

    register unsigned int x;

    l = 0;
    m = size - 1;

    while (l < m)
    {
        x = array[k];
        i = l;
        j = m;

        do
        {
            while (array[i] < x)
                ++i;

            while (x < array[j])
                --j;

            if (i <= j)
            {
                TGM_SWAP(array[i], array[j], int);
                ++i;
                --j;
            }

        }while (i <= j);

        if (j < k)
            l = i;

        if (k < i)
            m = j;
    }

    return array[k];
}

unsigned int FindMedianUint(unsigned int array[], unsigned int size)
{
    int median = FindKthSmallestUint(array, size, size / 2);
    if (size % 2 == 0)
    {
        unsigned int lowerMedian = 0;
        for (unsigned int i = 0; i != size / 2; ++i)
        {
            if (array[i] > lowerMedian)
                lowerMedian = array[i];
        }

        median = (unsigned int) ((lowerMedian + median) / 2.0);
    }

    return median;
}

char* TGM_CreateFileName(const char* workingDir, const char* fileName)
{
    unsigned int dirLen = strlen(workingDir);
    unsigned int nameLen = strlen(fileName);

    if (dirLen < 1 || nameLen < 1)
        return NULL;

    char* fullName = (char*) calloc(sizeof(char), dirLen + nameLen + 2);
    if (fullName == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for the file name.\n");

    strncpy(fullName, workingDir, dirLen);
    if (workingDir[dirLen - 1] == '/')
        strncpy(fullName + dirLen, fileName, nameLen);
    else
    {
        fullName[dirLen] = '/';
        strncpy(fullName + dirLen + 1, fileName, nameLen);
    }

    return fullName;
}

unsigned int TGM_TrimSpaces(char* str)
{
    const char* start = str;
    const char* end = str + strlen(str) - 1;

    while (*start == ' ' || *start == '\t' || *start == '\n')
        ++start;

    while (*end == ' ' || *end == '\t' || *end == '\n')
        --end;

    if (end <= start)
        return 0;

    unsigned strLen = end - start + 1;
    memmove(str, start, strLen);

    str[strLen] = '\0';
    return strLen;
}

TGM_Status TGM_GetNextLine(char* buff, unsigned int buffSize, FILE* input)
{
    while (fgets(buff, buffSize, input) != NULL)
    {
        unsigned int len = TGM_TrimSpaces(buff);
        if (len > 0)
            return TGM_OK;
        else
            continue;
    }

    return TGM_EOF;
}

int TGM_GetNumMismatchFromBam(const bam1_t* pAlgn)
{
    int numMM = 0;
    uint32_t* cigar = bam1_cigar(pAlgn);
    for (unsigned i = 0; i != pAlgn->core.n_cigar; ++i)
    {
        int type = (cigar[i] & BAM_CIGAR_MASK);
        if (type == BAM_CINS || type == BAM_CDEL)
            numMM += (cigar[i] >> BAM_CIGAR_SHIFT);
    }

    uint8_t* mdPos = bam_aux_get(pAlgn, "MD");
    if (mdPos != NULL)
    {
        const char* mdStr = bam_aux2Z(mdPos);
        const char* mdFieldPos = mdStr;
        while (mdFieldPos != NULL && *mdFieldPos != '\0')
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
    }

    return numMM;
}

static TGM_Status TGM_DoDir(const char* path, mode_t mode)
{
    struct stat     st;
    TGM_Status status = TGM_OK;

    if (stat(path, &st) != 0)
    {
        /*  Directory does not exist */
        if (mkdir(path, mode) != 0)
            status = TGM_ERR;
    }
    else if (!S_ISDIR(st.st_mode))
    {
        errno = ENOTDIR;
        status = TGM_ERR;
    }

    return status;
}

/* *
 * ** mkpath - ensure all directories in path exist
 * ** Algorithm takes the pessimistic view and works top-down to ensure
 * ** each directory in path exists, rather than optimistically creating
 * ** the last element and working backwards.
 * */
static TGM_Status TGM_MakePath(const char* path, mode_t mode)
{
    char           *pp;
    char           *sp;
    TGM_Status      status;
    char           *copypath = strdup(path);

    status = 0;
    pp = copypath;

    while (status == TGM_OK && (sp = strchr(pp, '/')) != 0)
    {
        if (sp != pp)
        {
            /*  Neither root nor double slash in path */
            *sp = '\0';
            status = TGM_DoDir(copypath, mode);
            *sp = '/';
        }

        pp = sp + 1;
    }

    if (status == 0)
        status = TGM_DoDir(path, mode);

    free(copypath);
    return status;
}

TGM_Bool TGM_IsDir(const char* path)
{
    struct stat st;

    if (stat(path, &st) != 0)
        return FALSE;
    
    if (!S_ISDIR(st.st_mode))
        return FALSE;

    return TRUE;
}

// make sure the working dir is empty
TGM_Status TGM_CheckWorkingDir(const char* workingDir)
{
    TGM_Status status = TGM_MakePath(workingDir, (S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH));

    if (status == TGM_OK)
    {
        DIR* pDir = opendir(workingDir);
        if (pDir != NULL)
        {
            status = TGM_ERR;
            struct dirent* pDirRecord = NULL;

            while ((pDirRecord = readdir(pDir)) != NULL)
            {
                status = TGM_OK;
                if (strcmp(pDirRecord->d_name, ".") != 0 && strcmp(pDirRecord->d_name, "..") != 0)
                {
                    status = TGM_ERR;
                    TGM_ErrMsg("ERROR: Working directory \"%s\" is not empty.\n", workingDir);
                    break;
                }
            }

            closedir(pDir);
        }
        else
        {
            status = TGM_ERR;
            TGM_ErrMsg("ERROR: Cannot open working directory: %s.\n", workingDir);
        }
    }

    return status;
}

TGM_Status TGM_TransferFile(char* buffer, int buffSize, int64_t fileSize, FILE* input, FILE* output)
{
    unsigned int readSize = 0;
    unsigned int writeSize = 0;
    while (fileSize >= buffSize)
    {
        readSize = fread(buffer, sizeof(char), buffSize, input);
        if (readSize != buffSize)
            return TGM_ERR;

        writeSize = fwrite(buffer, sizeof(char), buffSize, output);
        if (writeSize != buffSize)
            return TGM_ERR;

        fileSize -= buffSize;
    }

    if (fileSize > 0)
    {
        readSize = fread(buffer, sizeof(char), fileSize, input);
        if (readSize != fileSize)
            return TGM_ERR;

        writeSize = fwrite(buffer, sizeof(char), fileSize, output);
        if (writeSize != fileSize)
            return TGM_ERR;
    }

    return TGM_OK;
}
