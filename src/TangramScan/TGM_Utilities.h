/*
 * =====================================================================================
 *
 *       Filename:  TGM_Utilities->h
 *
 *    Description:  
 *
 *        Version:  1->0
 *        Created:  07/10/2011 03:01:14 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef  TGM_UTILITIES_H
#define  TGM_UTILITIES_H

#include <math.h>
#include <ctype.h>
#include <stdio.h>
#include <stdint.h>

#include "bam.h"
#include "TGM_Types.h"



#define TGM_TO_STRING_PRE(obj) #obj

#define TGM_TO_STRING(obj) TGM_TO_STRING_PRE(obj)

#define TGM_SWAP(objectOne, objectTwo, type)      \
    do                                           \
    {                                            \
        type temp;                               \
        temp = (objectOne);                      \
        (objectOne) = (objectTwo);               \
        (objectTwo) = temp;                      \
                                                 \
    }while(0)


//=======================
// Array utilities
//=======================

#define TGM_ARRAY_GET(pArray, i) ((pArray)->data[(i)])

#define TGM_ARRAY_GET_PT(pArray, i) ((pArray)->data + (i))

#define TGM_ARRAY_GET_SIZE(pArray) ((pArray)->size)

#define TGM_ARRAY_GET_FIRST(pArray) ((pArray)->data[0])

#define TGM_ARRAY_GET_FIRST_PT(pArray) ((pArray)->data)

#define TGM_ARRAY_GET_LAST(pArray) ((pArray)->data[(pArray)->size - 1])

#define TGM_ARRAY_GET_LAST_PT(pArray) ((pArray)->data + (pArray)->size - 1)

#define TGM_ARRAY_ALLOC(pArray, pArrayCap, objType, dataType)              \
    do                                                                    \
    {                                                                     \
        (pArray) = (objType*) malloc(sizeof(objType));                    \
        if ((pArray) == NULL)                                             \
            TGM_ErrSys("ERROR: not enough memory for an object->\n");      \
                                                                          \
        TGM_ARRAY_INIT((pArray), pArrayCap, dataType);                     \
                                                                          \
    }while(0)

#define TGM_ARRAY_INIT(pArray, pArrayCap, dataType)                           \
    do                                                                       \
    {                                                                        \
        (pArray)->size = TGM_EMPTY;                                           \
        (pArray)->capacity = pArrayCap;                                      \
                                                                             \
        (pArray)->data = (dataType*) calloc(pArrayCap, sizeof(dataType));    \
        if ((pArray)->data == NULL)                                          \
            TGM_ErrSys("ERROR: not enough memory for the data storage->\n");  \
                                                                             \
    }while(0)


#define TGM_ARRAY_PUSH(pArray, pNewElt, dataType)                                                              \
    do                                                                                                        \
    {                                                                                                         \
        if ((pArray)->size == (pArray)->capacity)                                                             \
        {                                                                                                     \
            (pArray)->capacity *= 2;                                                                          \
            (pArray)->data = (dataType*) realloc((pArray)->data, sizeof(dataType) * (pArray)->capacity);      \
            if ((pArray)->data == NULL)                                                                       \
                TGM_ErrSys("ERROR: not enough memory to expand the storage of data in an array.\n");           \
        }                                                                                                     \
                                                                                                              \
        TGM_ARRAY_GET(pArray, (pArray)->size) = *(pNewElt);                                                    \
        ++((pArray)->size);                                                                                   \
                                                                                                              \
    }while(0)

#define TGM_ARRAY_RESIZE(pArray, newCapacity, dataType)                                                        \
    do                                                                                                        \
    {                                                                                                         \
        if ((newCapacity) > (pArray)->capacity)                                                               \
        {                                                                                                     \
            (pArray)->capacity = (newCapacity);                                                               \
            (pArray)->data = (dataType*) realloc((pArray)->data, sizeof(dataType) * (pArray)->capacity);      \
            if ((pArray)->data == NULL)                                                                       \
                TGM_ErrSys("ERROR: not enough memory to expand the storage of data in an array.\n");           \
                                                                                                              \
        }                                                                                                     \
                                                                                                              \
    }while(0)

#define TGM_ARRAY_RESIZE_NO_COPY(pArray, newCapacity, dataType)                                                \
    do                                                                                                        \
    {                                                                                                         \
        if ((newCapacity) > (pArray)->capacity)                                                               \
        {                                                                                                     \
            free(pArray->data);                                                                               \
            (pArray)->capacity = newCapacity;                                                                 \
            (pArray)->data = (dataType*) calloc(sizeof(dataType), (pArray)->capacity);                        \
            if ((pArray)->data == NULL)                                                                       \
                TGM_ErrSys("ERROR: not enough memory to expand the storage of data in an array.\n");           \
                                                                                                              \
        }                                                                                                     \
                                                                                                              \
        (pArray)->size = 0;                                                                                   \
                                                                                                              \
    }while(0)

#define TGM_ARRAY_POP(pArray)   \
    do                         \
    {                          \
        (pArray)->size -= 1;   \
                               \
    }while(0)


#define TGM_ARRAY_RESET(pArray) (pArray)->size = TGM_EMPTY


#define TGM_ARRAY_FREE(pArray, freeObj)   \
    do                                   \
    {                                    \
        if ((pArray) != NULL)            \
        {                                \
            free((pArray)->data);        \
            if (freeObj)                 \
                free(pArray);            \
            else                         \
                (pArray)->data = NULL;   \
        }                                \
                                         \
    }while(0)


#define TGM_ARRAY_IS_FULL(pArray) ((pArray)->size == (pArray)->capacity)

int FindKthSmallestInt(int array[], int size, int k);

int FindMedianInt(int array[], int size);

unsigned int FindKthSmallestUint(unsigned int array[], unsigned int size, unsigned int k);

unsigned int FindMedianUint(unsigned int array[], unsigned int size);

static inline int DoubleRoundToInt(double number)
{
    return ((int) floor(number + 0.5));
}

static inline void StrToUpper(char* str)
{
    for (unsigned int i = 0; str[i] != '\0'; ++i)
        str[i] = toupper(str[i]);
}

char* TGM_CreateFileName(const char* workingDir, const char* fileName);

TGM_Status TGM_GetNextLine(char* buff, unsigned int buffSize, FILE* input);

int TGM_GetNumMismatchFromBam(const bam1_t* pAlgn);

static inline int TGM_CompareUint(const void* a, const void* b)
{
    const uint32_t* first = (const uint32_t*) a;
    const uint32_t* second = (const uint32_t*) b;

    if (*first < *second)
        return -1;
    else if (*first > *second)
        return 1;
    else
        return 0;
}

TGM_Bool TGM_IsDir(const char* path);


//=================================================================
// function:
//      write the read pair table into the read pair files
//
// args:
//      1. pReadPairTable: a pointer to a read pair table
//      2. pFileHash: a pointer to a hash table of read pair files
//=================================================================
TGM_Status TGM_CheckWorkingDir(const char* workingDir);

TGM_Status TGM_TransferFile(char* buffer, int buffSize, int64_t fileSize, FILE* input, FILE* output);

#endif  /*TGM_UTILITIES_H*/
