/*
 * =====================================================================================
 *
 *       Filename:  SR_InHashTable.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/12/2011 22:32:41
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */
#ifndef  SR_INHASHTABLE_H
#define  SR_INHASHTABLE_H

#include <stdio.h>
#include <stdint.h>

#include "SR_Types.h"
#include "SR_Reference.h"


//===============================
// Type and constant definition
//===============================

// generate a mask to clear the highest 2 bits in a hash key (the leftmost base pair)
#define GET_HIGH_END_MASK(hashSize) ((uint32_t) 0xffffffff >> (34 - (2 * (hashSize))))

typedef struct HashPosView
{
    const uint32_t* data;     // the address where a certain hash start in the "hashPos" array in the "SR_InHashTable" object 

    unsigned int size;        // total number of hash position found in the reference for a certain hash 

}HashPosView;

// input format of reference hash table
typedef struct SR_InHashTable
{
    int32_t id;                    // chromosome of the current reference hash table

    unsigned char hashSize;        // size of hash

    uint32_t* hashPos;             // positions of hashes found in the reference sequence

    uint32_t* indices;             // index of a given hash in the "hashPos" array

    uint32_t  highEndMask;         // a mask to clar the highest 2 bits in a hash key

    uint32_t  numPos;              // total number of hash positions found in reference

    uint32_t  numHashes;           // total number of different hashes

}SR_InHashTable;


//===============================
// Constructors and Destructors
//===============================

SR_InHashTable* SR_InHashTableAlloc(unsigned char hashSize);

void SR_InHashTableFree(SR_InHashTable* pHashTable);


//===============================
// Interface functions
//===============================

//================================================================
// function:
//      read the start part of the hash table file, including the
//      reference header position and the hash size
//
// args:
//      1. pHashSize: a pointer to the hash size
//      2. htInput: a file pointer to the hash table input file
// 
// return:
//      the reference header position (secrete code for 
//      compatibility check with the reference file)
//================================================================ 
int64_t SR_InHashTableReadStart(unsigned char* pHashSize, FILE* htInput);

//============================================================================
// function:
//      read the hash positions of the special sequence into the 
//      hash table structure
//
// args:
//      1. pSpecialHashTable: a pointer to the special hash table structure
//      2. pRefHeader: a pointer to a reference header structure
//      3. htInput: a file pointer to the hash table input file
// 
// return:
//      SR_OK: read successfully
//      SR_ERR: found an error during reading
//============================================================================
SR_Status SR_InHashTableReadSpecial(SR_InHashTable* pSpecialHashTable, const SR_RefHeader* pRefHeader, FILE* htInput);

//==================================================================
// function:
//      jump to a certain chromosome in the hash table file given
//      the reference ID
//
// args:
//      1. htInput: a file pointer to the hash table input file
//      2. pRefHeader: a pointer to the reference header structure
//      3. refID: the reference ID
// 
// return:
//      SR_OK: jump successfully
//      SR_ERR: found an error during jump
//==================================================================
SR_Status SR_InHashTableJump(FILE* htInput, const SR_RefHeader* pRefHeader, int32_t refID);

//==================================================================
// function:
//      read the hash positions of a chromosome into the hash table
//      structure
//
// args:
//      1. pHashTable: a pointer to the hash table structure
//      2. htInput: a file pointer to the hash table input file
// 
// return:
//      SR_OK: read successfully
//      SR_ERR: found an error during reading
//==================================================================
SR_Status SR_InHashTableRead(SR_InHashTable* pHashTable, FILE* htInput);

//======================================================================
// function:
//      get the hash position array of a given hash key
//
// args:
//      1. pHashPosView: a pointer to the hash position view structure
//      2. pHashTable: a pointer to the hash table structure
//      3. hashKey: hash key
// 
// return:
//      if the hash key is found in the hash table, hash position
//      structure will be loaded and TRUE will be returned; otherwise
//      FALSE is returned.
//======================================================================
SR_Bool SR_InHashTableSearch(HashPosView* pHashPosView, const SR_InHashTable* pHashTable, uint32_t hashKey);


#endif  /*SR_INHASHTABLE_H*/
