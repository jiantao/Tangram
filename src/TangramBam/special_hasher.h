#ifndef UTILITIES_HASHTABLE_SPECIAL_HASHER_H_
#define UTILITIES_HASHTABLE_SPECIAL_HASHER_H_

#include <string>
extern "C" {
#include "SR_InHashTable.h"
#include "SR_Reference.h"
}

class SpecialHasher {
 public:
  SpecialHasher(void);
  SpecialHasher(const char* fasta, 
                const int& hash_size = 7, 
		const int& ref_id_start_no = 0);
  ~SpecialHasher(void);

  // @function: Setting fasta filename.
  //            Notice that before Load(), the fasta filename should be set.
  //            If the filename is already given in the constructor,
  //            then you don't have to use this function.
  // @param:    fasta: fasta filename
  void SetFastaName(const char* fasta) {
    if (!is_loaded_) fasta_ = fasta;
  };

  // @function: Setting hash size.
  //            Notice that 1) Default hash size is 7; 
  //            2) after Load(), the hash size cannot be reset.
  //            The size also can be given in the constructor.
  void SetHashSize(const int& hash_size) {
    if (!is_loaded_) hash_size_ = hash_size;
  };

  // @function: Setting the start number of special references.
  //            Since special references are attached after the original references
  //            in the bam header, the start number of special references is 
  //            recommended to be set if the output format is bam.
  void SetRefIdStartNo(const int& start_no) {
    ref_id_start_no_ = start_no;
    if (reference_header_->pSpecialRefInfo) 
      reference_header_->pSpecialRefInfo->ref_id_start_no = ref_id_start_no_;
  };

  // @function: Loading special references from the fasta file 
  //            and hashing them.
  bool Load(void);

  const SR_Reference* GetReference(void) const {return(is_loaded_ ? references_ : NULL);};
  const SR_RefHeader* GetReferenceHeader(void) const {return(is_loaded_ ? reference_header_ : NULL);};
  const SR_InHashTable* GetHashTable(void) const {return(is_loaded_ ? hash_table_ : NULL);};

 private:
  std::string fasta_;
  SR_RefHeader* reference_header_;
  SR_Reference* references_;
  SR_InHashTable* hash_table_;
  int hash_size_;
  bool is_loaded_;
  int ref_id_start_no_;

  void Init(void);
  SpecialHasher (const SpecialHasher&);
  SpecialHasher& operator= (const SpecialHasher&);
};

#endif //UTILITIES_HASHTABLE_SPECIAL_HASHER_H_
