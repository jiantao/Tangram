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
    int length = strlen(str);

    const char* start = str;
    const char* end = str + length - 1;

    int count = 0;

    while (count != length)
    {
        if (*start == ' ' || *start == '\t' || *start == '\n')
            ++start;
        else
            break;

        ++count;
    }

    count = 0;
    while (count != length)
    {
        if (*end == ' ' || *end == '\t' || *end == '\n')
            --end;
        else
            break;

        ++count;
    }

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
