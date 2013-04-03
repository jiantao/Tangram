/*
 * =====================================================================================
 *
 *       Filename:  TGM_FragLenHist.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  02/08/2012 02:24:19 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */
#ifndef  TGM_FRAGLENHIST_H
#define  TGM_FRAGLENHIST_H

#include <stdint.h>
#include <stdio.h>

#include "TGM_Types.h"

#define INVALID_FRAG_LEN_QUAL 255

#define MIN_FRAGLEN_HIST_SIZE 10

// the object used to hold the fragment length histogram of a given read group
typedef struct TGM_FragLenHist
{
    void* rawHist;                  // raw fragment length histogram. this is a hash table for histogram building

    uint32_t* fragLen;              // array of the fragment length

    uint64_t* freq;                 // frequency of each fragment length

    double mean;                    // mean of the histogram

    double median;                  // median of the histogram

    double stdev;                   // standard deviation of the histogram

    uint32_t size;                  // number of unique fragment length

    uint32_t capacity;              // capacity of the histogram

    uint64_t modeCount[2];          // total counts of a histogram

}TGM_FragLenHist;

typedef struct TGM_FragLenHistArray
{
    TGM_FragLenHist* data;

    uint32_t size;

    uint32_t capacity;

}TGM_FragLenHistArray;

typedef struct TGM_FragLenHistLite
{
    uint32_t* fragLen;              // array of the fragment length

    uint64_t* freq;                 // frequency of each fragment length

    uint32_t size;

    uint32_t capacity;

}TGM_FragLenHistLite;

typedef struct TGM_FragLenHistLiteArray
{
    TGM_FragLenHistLite* data;

    uint32_t size;

    uint32_t capacity;

}TGM_FragLenHistLiteArray;

TGM_FragLenHistArray* TGM_FragLenHistArrayAlloc(unsigned int capacity);

void TGM_FragLenHistArrayFree(TGM_FragLenHistArray* pHistArray);

TGM_FragLenHistLite* TGM_FragLenHistLiteAlloc(unsigned int capacity);

void TGM_FragLenHistLiteFree(TGM_FragLenHistLite* pHistLite);

void TGM_FragLenHistLiteArrayFree(TGM_FragLenHistLiteArray* pHistLiteArray);

void TGM_FragLenHistArrayClear(TGM_FragLenHistArray* pHistArray);

void TGM_FragLenHistArrayInit(TGM_FragLenHistArray* pHistArray, unsigned int size);

TGM_Status TGM_FragLenHistArrayUpdate(TGM_FragLenHistArray* pHistArray, unsigned int backHistIndex, uint32_t fragLen);

int TGM_FragLenHistLiteArrayGetFragLenQual(const TGM_FragLenHistLiteArray* pHistArray, uint32_t refID, uint32_t fragLen, uint32_t median);

void TGM_FragLenHistArrayFinalize(TGM_FragLenHistArray* pHistArray);

void TGM_FragLenHistArrayWriteHeader(uint32_t size, FILE* output);

void TGM_FragLenHistArrayWrite(const TGM_FragLenHistArray* pHistArray, FILE* output);

TGM_FragLenHistLiteArray* TGM_FragLenHistLiteArrayRead(FILE* pHistArrayInput);

void TGM_FragLenHistTransfer(TGM_FragLenHistLite* pHistLite, FILE* input, FILE* output);

void TGM_FragLenHistLiteInit(TGM_FragLenHistLite* pHistLite, uint32_t newSize);

void TGM_FragLenHistLiteRead(TGM_FragLenHistLite* pHistLite, FILE* input);

void TGM_FragLenHistLiteWrite(const TGM_FragLenHistLite* pHistLite, FILE* output);

#endif  /*TGM_FRAGLENHIST_H*/
