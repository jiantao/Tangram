#include "special_hasher.h"

#include <stdlib.h>
#include <stdint.h>
#include <string.h>

extern "C" {
#include "SR_InHashTable.h"
#include "SR_OutHashTable.h"
#include "SR_Reference.h"
#include "ConvertHashTableOutToIn.h"
}

using std::string;

SpecialHasher::SpecialHasher(void)
    : fasta_("")
    , reference_header_(NULL)
    , references_(NULL)
    , hash_table_(NULL)
    , hash_size_(7)
    , is_loaded_(false)
    , ref_id_start_no_(0){
  Init();
}

SpecialHasher::SpecialHasher(const char* fasta, const int& hash_size, const int& ref_id_start_no)
    : fasta_(fasta)
    , reference_header_(NULL)
    , references_(NULL)
    , hash_table_(NULL)
    , hash_size_(hash_size)
    , is_loaded_(false)
    , ref_id_start_no_(ref_id_start_no){
  Init();
}

SpecialHasher::~SpecialHasher(void) {
  // SR_RefHeaderFree also frees memory 
  //   that is allocated by SR_SpecialRefInfoAlloc
  SR_RefHeaderFree(reference_header_);
  SR_ReferenceFree(references_);
  SR_InHashTableFree(hash_table_);
}

void SpecialHasher::Init(void) {
  // allocate memory for the header and special info
  reference_header_ = SR_RefHeaderAlloc(1,0);
  const int initial_reference_number = 30;
  reference_header_->pSpecialRefInfo = 
    SR_SpecialRefInfoAlloc(initial_reference_number);

  if (reference_header_->pSpecialRefInfo)
    reference_header_->pSpecialRefInfo->ref_id_start_no = ref_id_start_no_;

  references_ = SR_ReferenceAlloc();
}

bool SpecialHasher::Load(void) {
  if (fasta_.empty()) {
    fprintf(stderr, "ERROR: Please set fasta filename before loading.\n");
    return false;
  }
  // load the special references
  FILE* input = fopen(fasta_.c_str(), "r");
  if (input == NULL) {
    fprintf(stderr, "ERROR: The file (%s) cannot be opened.\n", fasta_.c_str());
  } else {
    SR_SpecialRefLoad(references_, reference_header_, input);
  }
  fclose(input);

  // index every possible hash position in the current chromosome
  // and write the results into hash position index file and hash position file
  SR_OutHashTable* out_hash_table = SR_OutHashTableAlloc(hash_size_);
  SR_OutHashTableLoad(out_hash_table, references_->sequence, 
                      references_->seqLen, references_->id);

  hash_table_ = SR_InHashTableAlloc(hash_size_);
  ConvertHashTableOutToIn(out_hash_table, hash_table_);

  SR_OutHashTableFree(out_hash_table);

  reference_header_->pSpecialRefInfo->ref_id_start_no = ref_id_start_no_;

  is_loaded_ = true;
  return true;
}
