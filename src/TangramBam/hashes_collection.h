#ifndef UTILITIES_MISCELLANEOUS_HASH_REGION_COLLECTION_H_ 
#define UTILITIES_MISCELLANEOUS_HASH_REGION_COLLECTION_H_

#include <vector>

//#include "utilities/hashTable/SR_HashRegionTable.h"
struct BestRegionArray;
struct BestRegion;

using std::vector;

namespace Scissors {
class HashesCollection {
 public:
  HashesCollection(){};
  ~HashesCollection(){};
  void Init(const BestRegionArray& array);
  void SortByLength(void);
  void Print(void)const;
  const BestRegion* Get (const unsigned int& index) const;
  inline const int GetSize() const;
  bool GetBestCoverPair(unsigned int* best1, unsigned int* best2);
  bool GetBestCoverPair(HashesCollection* hc, unsigned int* best1, unsigned int* best2);

 private:
  vector<BestRegion*> hash_regions_;
  HashesCollection (const HashesCollection&);
  HashesCollection& operator= (const HashesCollection&);
}; // end HashRegionCollection

inline const int HashesCollection::GetSize() const {
  return hash_regions_.size();
}
} //namespace Scissors
#endif  // UTILITIES_HASH_REGION_COLLECTION_H_
