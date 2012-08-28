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


#ifdef __cplusplus      //added by Jim Howard so that these functions can be called from c++
extern "C"
{
#endif

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

TGM_Status TGM_GetNextLine(char* buff, unsigned int buffSize, FILE* input);

#ifdef __cplusplus
}
#endif

#endif  /*TGM_UTILITIES_H*/
