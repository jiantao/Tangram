#include "seq_converter.h"
#include "../OutSources/samtools/bam.h"

static const char COMPLEMENT_BASE_MAP[16] = {0x0, 0x8, 0x4, 0xf, 0x2,0xf,0xf,0xf,0x1,0xf,0xf,0xf,0xf,0xf,0xf,0xf};

void GetReverseComplementSequence(const uint8_t* original, 
                                  const int length, 
              	 		  uint8_t* reverse) 
{
  int length_even = ((length % 2) == 1) ? length - 1 : length;
  for (int i = 0; i < length_even; i=i+2) {
    char upper = COMPLEMENT_BASE_MAP[bam1_seqi(original, length - i - 1)];
    char lower = COMPLEMENT_BASE_MAP[bam1_seqi(original, length - i - 2)];
    *reverse = (upper << 4) | lower;
    ++reverse;
  }

  if ((length % 2) == 1) {
    char upper = COMPLEMENT_BASE_MAP[bam1_seqi(original, 0)];
    *reverse = (upper << 4) & 0xf0;
    ++reverse;
  }
}

void GetComplementSequence(const uint8_t* original,
                           const int length,
			   uint8_t* reverse)
{
  for (int i = 0; i < length; i=i+2) {
    char upper = COMPLEMENT_BASE_MAP[bam1_seqi(original, i)];
    char lower = COMPLEMENT_BASE_MAP[bam1_seqi(original, i+1)];
    *reverse = (upper << 4) | lower;
    ++reverse;
  }
}

void GetInverseSequence(const uint8_t* original,
                        const int length,
			uint8_t* reverse)
{
  int length_even = ((length % 2) == 1) ? length - 1 : length;
  for (int i = 0; i < length_even; i=i+2) {
    char upper = bam1_seqi(original, length - i - 1);
    char lower = bam1_seqi(original, length - i - 2);
    *reverse = (upper << 4) | lower;
    ++reverse;
  }

  if ((length % 2) == 1) {
    char upper = bam1_seqi(original, 0);
    *reverse = (upper << 4) & 0xf0;
    ++reverse;
  }
}

void GetInverseQual(uint8_t* original,
                    const int length) {
  for (int i = 0, j = length - 1; i < j; ++i, --j) {
    uint8_t tmp = *(original + i);
    *(original + i) = *(original + j);
    *(original + j) = tmp;
  }
}
