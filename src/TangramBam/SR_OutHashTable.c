/*
 * =====================================================================================
 *
 *       Filename:  SR_OutHashTable.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  07/01/2011 08:35:32 PM
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
#include <assert.h>
#include <ctype.h>

#include "SR_Utilities.h"
#include "SR_OutHashTable.h"

#define DEFAULT_HASH_SIZE 7
#define DEFAULT_POS_ARR_CAPACITY 2000

// calculate the hash key for a hash in a read that starts at a certain position
static SR_Bool GetNextHashKey(uint32_t* hashKey, uint32_t* pos, const char* query, uint32_t queryLen, uint32_t mask, unsigned char hashSize)
{
    // table use to translate a nucleotide into its corresponding 2-bit representation
    static const char translation[26] = { 0, -1, 1, -1, -1, -1, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 3, -1, -1, -1, -1, -1, -1 };

    unsigned int endPos = *pos + hashSize;
    unsigned int startPos = *pos == 0 ? *pos : endPos - 1;


    while (endPos <= queryLen)
    {
        // remove the highest 2 bits in the previous hash key
        *hashKey &= mask;

        char tValue = 0;
        for (unsigned int i = startPos; i != endPos; ++i)
        {
            char cur_base = 0;
            if (islower(query[i])) cur_base = toupper(query[i]);
            else cur_base = query[i];
            tValue = translation[cur_base - 'A'];
            if (tValue < 0)
            {
                *hashKey = 0;

                *pos = i + 1;
                startPos = *pos;
                endPos = startPos + hashSize;
                break;
            }

            *hashKey = *hashKey << 2 | tValue;
        }

        if (tValue >= 0)
            return TRUE;
    }

    return FALSE;
}


SR_HashPosArray* SR_HashPosArrayAlloc(unsigned int capacity)
{
    SR_HashPosArray* newArray = (SR_HashPosArray*) malloc(sizeof(SR_HashPosArray));
    if (newArray == NULL)
        SR_ErrSys("ERROR: Not enough memory for a hash position array object.\n");

    if (capacity == 0)
        capacity = DEFAULT_POS_ARR_CAPACITY;

    newArray->data = (uint32_t*) malloc(sizeof(uint32_t) * capacity);
    if (newArray->data == NULL)
        SR_ErrSys("ERROR: Not enough memory for the storage of hash positions in a hash position array object.\n");

    newArray->size = 0;
    newArray->capacity = capacity;

    return newArray;
}

void SR_HashPosArrayInit(SR_HashPosArray* hashPosArray, unsigned int capacity)
{
    assert(hashPosArray != NULL);

    if (capacity == 0)
        capacity = DEFAULT_POS_ARR_CAPACITY;

    hashPosArray->data = (uint32_t*) malloc(sizeof(uint32_t) * capacity);
    if (hashPosArray->data == NULL)
        SR_ErrSys("ERROR: Not enough memory for the storage of hash positions in a hash position array object.\n");

    hashPosArray->size = 0;
    hashPosArray->capacity = capacity;
}

void SR_HashPosArrayFree(SR_HashPosArray* hashPosArray)
{
    if (hashPosArray != NULL)
    {
        free(hashPosArray->data);

        free(hashPosArray);
    }
}

void SR_HashPosArrayPushBack(SR_HashPosArray* hashPosArray, uint32_t pos)
{
    if (hashPosArray->size >= hashPosArray->capacity)
    {
        hashPosArray->capacity *= 2;
        hashPosArray->data = (uint32_t*) realloc(hashPosArray->data, sizeof(uint32_t) * hashPosArray->capacity);
        if (hashPosArray->data == NULL)
            SR_ErrSys("ERROR: Not enough memory for the storage of hash positions in a hash position array object.\n");
    }

    SR_ARRAY_GET(hashPosArray, hashPosArray->size) = pos;
    ++(hashPosArray->size);
}

SR_OutHashTable* SR_OutHashTableAlloc(unsigned char hashSize)
{
    if (hashSize > MAX_HASH_SIZE)
        SR_ErrQuit("ERROR: Hash size can not be greater than %d\n", MAX_HASH_SIZE);
    else if (hashSize == 0)
    {
        SR_ErrMsg("WARNING: Hash size should be greater than zero. Default hash size %d will be used.\n", DEFAULT_HASH_SIZE);
        hashSize = DEFAULT_HASH_SIZE;
    }

    SR_OutHashTable* newTable = (SR_OutHashTable*) malloc(sizeof(SR_OutHashTable));
    if (newTable == NULL)
        SR_ErrSys("ERROR: Not enough memory for a reference hash table object.\n");
    
    newTable->id = 0;
    newTable->hashSize = hashSize;
    newTable->numPos = 0;
    newTable->numHashes = (uint32_t) 1 << (2 * hashSize);

    newTable->hashPosTable = (SR_HashPosArray*) malloc(sizeof(SR_HashPosArray) * newTable->numHashes);
    if (newTable->hashPosTable == NULL)
        SR_ErrSys("ERROR: Not enough memory for the storage of the hash position table in a reference hash table object.\n");

    for (unsigned int i = 0; i != newTable->numHashes; ++i)
        SR_HashPosArrayInit(&((newTable->hashPosTable)[i]), DEFAULT_POS_ARR_CAPACITY);

    return newTable;
}


void SR_OutHashTableFree(SR_OutHashTable* pHashTable)
{
    if (pHashTable != NULL)
    {
        if (pHashTable->hashPosTable != NULL)
        {
            for (unsigned int i = 0; i != pHashTable->numHashes; ++i)
                free((pHashTable->hashPosTable[i]).data);
        }

        free(pHashTable->hashPosTable);
        free(pHashTable);
    }
}

void SR_OutHashTableLoad(SR_OutHashTable* pHashTable, const char* refSeq, uint32_t refLen, int32_t id)
{
    uint32_t mask = GetHighEndMask(pHashTable->hashSize);
    uint32_t hashKey = 0;
    uint32_t pos = 0;

    pHashTable->id = id;
    while (GetNextHashKey(&hashKey, &pos, refSeq, refLen, mask, pHashTable->hashSize))
    {
        SR_HashPosArrayPushBack(&((pHashTable->hashPosTable)[hashKey]), pos);
        ++(pHashTable->numPos);
        ++pos;
    }
}

int64_t SR_OutHashTableWrite(const SR_OutHashTable* pHashTable, FILE* htOutput)
{
    int64_t fileOffset = ftello(htOutput);
    if (fileOffset == -1)
        SR_ErrSys("ERROR: Cannot get the offset of current file.\n");

    size_t writeSize = 0;

    writeSize = fwrite(&(pHashTable->id), sizeof(int32_t), 1, htOutput);
    if (writeSize != 1)
        SR_ErrSys("ERROR: Cannot write the chromosome ID to the hash table file.\n");

    uint32_t index = 0;
    for (unsigned int i = 0; i != pHashTable->numHashes; ++i)
    {
        writeSize = fwrite(&index, sizeof(uint32_t), 1, htOutput);
        if (writeSize != 1)
            SR_ErrSys("ERROR: Cannot write hash position index to the hash table file.\n");

        uint32_t hashPosSize = (pHashTable->hashPosTable)[i].size;
        index += hashPosSize;
    }

    writeSize = fwrite(&(pHashTable->numPos), sizeof(uint32_t), 1, htOutput);
    if (writeSize != 1)
        SR_ErrSys("ERROR: Cannot write the total number of hash positions to the hash table file.\n");

    for (unsigned int i = 0; i != pHashTable->numHashes; ++i)
    {
        const uint32_t* hashPos = (pHashTable->hashPosTable)[i].data;
        uint32_t hashPosSize = (pHashTable->hashPosTable)[i].size;

        writeSize = fwrite(hashPos, sizeof(uint32_t), hashPosSize, htOutput);
        if (writeSize != hashPosSize)
            SR_ErrSys("ERROR: Cannot write hash position to the hash table file.\n");

    }

    fflush(htOutput);

    return fileOffset;
}

void SR_OutHashTableWriteStart(unsigned char hashSize, FILE* htOutput)
{
    size_t writeSize = 0;                                                                
    int64_t emptyOffset = 0;

    writeSize = fwrite(&emptyOffset, sizeof(emptyOffset), 1, htOutput);
    if (writeSize != 1)                                                                  
        SR_ErrQuit("ERROR: Cannot write the offset of reference header into hash table file.\n");

    writeSize = fwrite(&(hashSize), sizeof(unsigned char), 1, htOutput);           
    if (writeSize != 1)                                                                  
        SR_ErrQuit("ERROR: Cannot write the hash size into hash table file.\n");

    fflush(htOutput);
}

void SR_OutHashTableSetStart(int64_t refHeaderPos, FILE* htOutput)
{
    size_t writeSize = 0;                                                                

    if (fseeko(htOutput, 0, SEEK_SET) != 0)
        SR_ErrQuit("ERROR: Cannot seek in the hash table file.\n");

    writeSize = fwrite(&refHeaderPos, sizeof(refHeaderPos), 1, htOutput);
    if (writeSize != 1)                                                                  
        SR_ErrQuit("ERROR: Cannot write the offset of reference header into hash table file.\n");

    fflush(htOutput);
}


void SR_OutHashTableReset(SR_OutHashTable* pHashTable)
{
    pHashTable->id = 0;
    pHashTable->numPos = 0;

    for (unsigned int i = 0; i != pHashTable->numHashes; ++i)
    {
        SR_HashPosArray* currPosArray = &((pHashTable->hashPosTable)[i]);
        SR_ARRAY_RESET(currPosArray);
    }
}
