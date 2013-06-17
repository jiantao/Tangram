#ifndef UTILITIES_HASHTABLE_REFERENCE_HASHER_H_
#define UTILITIES_HASHTABLE_REFERENCE_HASHER_H_

#include <string>
extern "C" {
#include "SR_InHashTable.h"
#include "SR_Reference.h"
}

class ReferenceHasher {
 public:
  ReferenceHasher(void);
  ReferenceHasher(const char* sequence);
  ~ReferenceHasher(void);

  // @function: Setting sequence
  //            Notice that before Load(), the fasta filename should be set.
  //            If the filename is already given in the constructor,
  //            then you don't have to use this function.
  // @param:    fasta: fasta filename
  bool SetSequence(const char* fasta);

  // @function: Setting hash size.
  //            Notice that 1) Default hash size is 7; 
  //            2) before Load(), the hash size should be set.
  //            The size also can be given in the constructor.
  void SetHashSize(const int& hash_size) {hash_size_ = hash_size;};

  // @function: Loading special references from the fasta file 
  //            and hashing them.
  bool Load(void);

  void Clear(void);

  const SR_Reference* GetReference(void) const {return(is_loaded_ ? references_ : NULL);};
  //const SR_RefHeader* GetReferenceHeader(void) const {return(is_loaded ? reference_header_ : NULL);};
  const SR_InHashTable* GetHashTable(void) const {return(is_loaded_ ? hash_table_ : NULL);};

 private:
  //std::string fasta_;
  //SR_RefHeader* reference_header_;
  SR_Reference* references_;
  SR_InHashTable* hash_table_;
  int hash_size_;
  bool is_loaded_;

  void Init(void);
  ReferenceHasher (const ReferenceHasher&);
  ReferenceHasher& operator= (const ReferenceHasher&);
};

#endif //UTILITIES_HASHTABLE_SPECIAL_HASHER_H_
