#include "hashes_collection.h"

#include <limits.h>
#include <stdio.h>
#include <algorithm>

#include "SR_HashRegionTable.h"

bool OperatorLength(const BestRegion* r1, const BestRegion* r2) {
  return (r1->length) < (r2->length);
}

bool OperatorQueryBegin(const BestRegion* r1, const BestRegion* r2) {
  return (r1->queryBegin) < (r2->queryBegin);
}

bool OperatorQueryBeginDescending(const BestRegion* r1, const BestRegion* r2) {
  return (r2->queryBegin) < (r1->queryBegin);
}

namespace Scissors {

void HashesCollection::Init(const BestRegionArray& array) {
  hash_regions_.clear();
  hash_regions_.resize(array.size);
  for (unsigned int i = 0; i < array.size; ++i) 
    hash_regions_[i] = &array.data[i];
}

void HashesCollection::SortByLength() {
  sort(hash_regions_.begin(), hash_regions_.end(), OperatorLength);
}

bool HashesCollection::GetBestCoverPair(HashesCollection* hc, 
                                        unsigned int* best1, 
					unsigned int* best2) {
  sort(hash_regions_.begin(), hash_regions_.end(), OperatorQueryBegin);
  sort(hc->hash_regions_.begin(), hc->hash_regions_.end(), OperatorQueryBeginDescending);

  int best_cover = INT_MIN;
  bool found = false;
  for (unsigned int i = 0; i < hash_regions_.size(); ++i) {
    if (hash_regions_[i]->length == 0) continue;;
    int end_i = hash_regions_[i]->queryBegin + hash_regions_[i]->length - 1;
    for (unsigned int j = 0; j < hc->hash_regions_.size(); ++j) {
      if (hc->hash_regions_[j]->length == 0) continue;
      if (hash_regions_[i]->queryBegin > hc->hash_regions_[j]->queryBegin) break;

      int overlap = end_i - hc->hash_regions_[j]->queryBegin + 1;
      int cover   = hash_regions_[i]->length + hc->hash_regions_[j]->length;
      if (overlap > 0) cover -= overlap;

      // update the best cover
      if (best_cover < cover) {
        best_cover = cover;
	*best1 = i;
	*best2 = j;
	found = true;
      } else {
        // nothing
      }
    }
  }

  return found;
}

// @description:
//   Get a pair of BestRegions which cover the read longest.
//   In other words, get a pair of the longest summed lengths.
// @params:
//   The ids of the alignments of the pair are assigned in best1 and best2.
// @returns:
//   Returns true when a pair is found; otherwise false.
bool HashesCollection::GetBestCoverPair(unsigned int* best1, unsigned int* best2) {
  // first, sort by query begins
  sort(hash_regions_.begin(), hash_regions_.end(), OperatorQueryBegin);

  int best_cover = INT_MIN;
  bool found = false;
  for (unsigned int i = 0; i < hash_regions_.size(); ++i) {
    int end_i = hash_regions_[i]->queryBegin + hash_regions_[i]->length - 1;
    for (unsigned int j = hash_regions_.size() - 1; j > i; --j) {
      int overlap = end_i - hash_regions_[j]->queryBegin + 1;
      int cover   = hash_regions_[i]->length + hash_regions_[j]->length;
      if (overlap > 0) cover -= overlap;
      
      // update the best cover
      if (best_cover < cover) {
        best_cover = cover;
	*best1 = i;
	*best2 = j;
	found = true;
      } else {
        // nothing
      } // end if-else
    } // end for
  } // end for

  return found;
}

const BestRegion* HashesCollection::Get (const unsigned int& index) const {
  if (index >= hash_regions_.size()) return NULL;

  return hash_regions_[index];
}

void HashesCollection::Print(void) const {
  for (vector<BestRegion*>::const_iterator ite = hash_regions_.begin();
      ite != hash_regions_.end(); ++ite) {
    fprintf(stderr, "%u %u %u %u %u %u\n", (*ite)->refBegins[0], (*ite)->refBegins[1], (*ite)->refBegins[2], (*ite)->queryBegin, (*ite)->length, (*ite)->numPos);
  }
}
} //namespace Scissors
