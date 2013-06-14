#ifndef TANGRAM_BAM_H_
#define TANGRAM_BAM_H_

#include <string>

#include "api/BamAlignment.h"
#include "../OutSources/stripedSW/ssw_cpp.h"

using namespace std;

const StripedSmithWaterman::Filter kFilter(true, false, 0, 32467);
const int kAlignmentMapSize = 10000;
const float kSoftClipRate = 0.15; // the max ratio of allowed soft clips

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
  string in_bam;
  string out_bam;
  string ref_fasta;
  string command_line;
  string target_ref_name; // the target chromosome

  Param()
      : in_bam()
      , out_bam()
      , ref_fasta()
      , command_line()
      , target_ref_name("-1")
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
