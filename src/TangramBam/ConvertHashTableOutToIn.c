#include <stdlib.h>
#include <string.h>
#include "SR_InHashTable.h"
#include "SR_OutHashTable.h"

void ConvertHashTableOutToIn(const SR_OutHashTable* out, SR_InHashTable* in) {
  in->id = out->id;

  uint32_t index = 0;
  for (uint32_t i = 0; i != out->numHashes; ++i) {
    in->indices[i] = index;
    index += (out->hashPosTable)[i].size;
  }

  in->numPos = out->numPos;

  free(in->hashPos);
  in->hashPos = (uint32_t*) malloc(sizeof(uint32_t) * in->numPos);
  uint32_t* ptr = in->hashPos;
  for (uint32_t i = 0; i != out->numHashes; ++i) {
    const uint32_t* hashPos = (out->hashPosTable)[i].data;
    const uint32_t hashPosSize = (out->hashPosTable)[i].size;
    memcpy(ptr, hashPos, sizeof(uint32_t) * hashPosSize);
    ptr += hashPosSize;
  }
}
