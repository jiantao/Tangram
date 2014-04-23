#ifndef TANGRAM_BAM_H_
#define TANGRAM_BAM_H_

#include <string>

#include "api/BamAlignment.h"
#include "../OutSources/stripedSW/ssw_cpp.h"

using namespace std;

struct Alignment {
  BamTools::BamAlignment bam_alignment;
  bool hit_insertion;
  string ins_prefix;

  Alignment()
      : bam_alignment()
      , hit_insertion(false)
      , ins_prefix()
  {}

  void Clear() {
    hit_insertion = false;
    ins_prefix.clear();
  }
};

struct Param {
  string in_bam; // -i
  string out_bam; // -o
  string ref_fasta; // -r
  string command_line;
  string target_ref_name; // -t, the target chromosome
  int required_match; // -m

  Param()
      : in_bam("stdin")
      , out_bam("stdout")
      , ref_fasta()
      , command_line()
      , target_ref_name("-1")
      , required_match(50)
  {}
};

struct SpecialReference {
  string concatnated;
  int concatnated_len;
  vector<int> ref_lens;
  vector<string> ref_names;

  SpecialReference()
      : concatnated()
      , concatnated_len(0)
      , ref_lens()
      , ref_names()
  {}
};

#endif
