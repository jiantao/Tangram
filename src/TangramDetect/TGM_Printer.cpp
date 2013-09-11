/*
 * =====================================================================================
 *
 *       Filename:  TGM_Printer.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  08/16/2012 03:35:16 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include <time.h>
#include <string>
#include <iomanip>
#include "TGM_Printer.h"

using namespace std;
using namespace Tangram;


Printer::Printer(const Detector* pDetector, const DetectPars& detectPars, const Aligner* pAligner, const Reference* pRef, 
                 const LibTable& libTable, const BamPairTable& bamPairTable, const GenotypePars& genotypePars, Genotype& genotype)
                : pDetector(pDetector), detectPars(detectPars), pAligner(pAligner), pRef(pRef), libTable(libTable), 
                  bamPairTable(bamPairTable), genotypePars(genotypePars), genotype(genotype)
{
  #ifdef TD_VERBOSE_DEBUG
  fprintf(stderr, "familyMap:\n");
  for (unsigned int i = 0; i < pRef->familyMap.Size(); ++i) {
    fprintf(stderr, "%d: %d\n", i, pRef->familyMap[i]);
  }
  fprintf(stderr, "familyName:\n");
  for (unsigned int i = 0; i < pRef->familyName.size(); ++i) {
    fprintf(stderr, "%d: %s\n", i, pRef->familyName[i].c_str());
  }
  #endif
}

Printer::~Printer()
{
    CloseOutputGrp();
}

void Printer::Init(void)
{
    InitOutputGrp();
    
    InitPrintSubset();

    InitPrintElmnts();
}

void Printer::Print(void)
{
    PrintHeader();

    PrintElmnt element;
    PrintElmnt temp;

    while (printElmnts.size() > 0)
    {
        element = *(printElmnts.begin());
        bool hasGenotype = false;
        switch(element.svType)
        {
            case SV_SPECIAL:
                PrintSpecial(element);
                hasGenotype = genotype.Special(element.pRpSpecial, element.pSplitEvent);

                if (genotypePars.doGenotype)
                    PrintGenotype(genotype, hasGenotype);
                else
                    PrintSampleInfo(genotype);
                break;
            case SV_INVERSION:
                break;
            default:
                break;
        }

        int road = element.subsetIdx;
        printElmnts.erase(printElmnts.begin());

        if (GetNextPrintElmnt(temp, road))
            printElmnts.insert(temp);
    }
}

void Printer::InitPrintSubset(void)
{
    unsigned int numSp = libTable.GetNumSpecialRef();
    unsigned int numFamily = 0;

    if (pAligner != NULL)
        numFamily = pRef->familyName.size();

    unsigned int numMax = numSp + numFamily;
    printSubset.Init(numMax);
     
    unsigned int subsetIdx = 0;

    for (unsigned int i = 0; i != numSp; ++i)
    {
        unsigned int numEvents = pDetector->pSpecialEventsTable[i].Size();
        const SpecialEvent* pRpSpecials = NULL;

        if (numEvents > 0)
            pRpSpecials = pDetector->pSpecialEventsTable[i].GetPointer(0);

        const SplitEvent* pSplitSpecials = NULL;
        unsigned int splitLen = 0;

        if (pAligner != NULL)
            pSplitSpecials = pAligner->GetSpecialStartFromZA(splitLen, i);

        if (numEvents == 0 && splitLen == 0)
            continue;

        printSubset[subsetIdx].pRpSpecials = pRpSpecials;
        printSubset[subsetIdx].rpSize = numEvents;
        printSubset[subsetIdx].rpIdx = 0;

        printSubset[subsetIdx].pSplitEvents = pSplitSpecials;
        printSubset[subsetIdx].srSize = splitLen;
        printSubset[subsetIdx].srIdx = 0;

        printSubset[subsetIdx].svType = SV_SPECIAL;

        ++subsetIdx;
        printSubset.Increment();
    }

    if (pAligner != NULL)
    {
        unsigned int numFamily = pRef->familyName.size();

        for (unsigned int i = 0; i != numFamily; ++i)
        {
            int zaID = pRef->familyToZA[i];
            if (zaID < 0)
            {
                unsigned int len = 0;
                const SplitEvent* pSplitSpecials = pAligner->GetSpecialStartFromFamily(len, i);

                if (len == 0)
                    continue;

                printSubset[subsetIdx].pRpSpecials = NULL;
                printSubset[subsetIdx].rpSize = 0;
                printSubset[subsetIdx].rpIdx = 0;

                printSubset[subsetIdx].pSplitEvents = pSplitSpecials;
                printSubset[subsetIdx].srSize = len;
                printSubset[subsetIdx].srIdx = 0;

                printSubset[subsetIdx].svType = SV_SPECIAL;

                ++subsetIdx;
                printSubset.Increment();
            }
        }
    }
}

void Printer::InitPrintElmnts(void)
{
    unsigned int subsetSize = printSubset.Size();

    PrintElmnt element;
    for (unsigned int i = 0; i != subsetSize; ++i)
    {
        if (GetNextPrintElmnt(element, i));
            printElmnts.insert(element);
    }
}

bool Printer::GetNextPrintElmnt(PrintElmnt& element, int subsetIdx)
{
    bool ret = false;
    switch (printSubset[subsetIdx].svType)
    {
        case SV_SPECIAL:
            ret = GetNextSpecial(element, subsetIdx);
            break;
        case SV_INVERSION:
            break;
        default:
            break;
    }

    return ret;
}

bool Printer::GetNextSpecial(PrintElmnt& element, int subsetIdx)
{
    PrintSubset& subset = printSubset[subsetIdx];

    // find the first read-pair event that is not merged with any split-read event
    element.pRpSpecial = NULL;
    while (subset.rpIdx < subset.rpSize)
    {
        if (subset.pRpSpecials[subset.rpIdx].splitIdx >= 0)
            ++(subset.rpIdx);
        else
        {
            element.pRpSpecial = subset.pRpSpecials + subset.rpIdx;
            break;
        }
    }

    if (subset.srIdx == subset.srSize)
        element.pSplitEvent = NULL;
    else
        element.pSplitEvent = subset.pSplitEvents + subset.srIdx;

    if (element.pRpSpecial == NULL && element.pSplitEvent == NULL)
        return false;
    if (element.pRpSpecial != NULL && element.pSplitEvent != NULL)
    {
        if (element.pRpSpecial->pos <= (unsigned int) element.pSplitEvent->pos)
            element.pSplitEvent = NULL;
        else
            element.pRpSpecial = NULL;
    }

    if (element.pRpSpecial != NULL)
    {
        element.pos = element.pRpSpecial->pos;
        ++(subset.rpIdx);
    }

    if (element.pSplitEvent != NULL)
    {
        element.pos = element.pSplitEvent->pos;

        int rpIdx = element.pSplitEvent->rpIdx;
        if (rpIdx >= 0)
        {
            element.pRpSpecial = subset.pRpSpecials + rpIdx;
        }

        ++(subset.srIdx);
    }

    element.subsetIdx = subsetIdx;
    element.svType = SV_SPECIAL;

    return true;
}

void Printer::PrintHeader(void)
{
    for (unsigned int i = SV_DELETION; i <= SV_INTER_CHR_TRNSLCTN; ++i)
    {
        unsigned int shiftSize = i - 1;
        string outputFile;

        switch(i)
        {
            case SV_DELETION:
                break;
            case SV_TANDEM_DUP:
                break;
            case SV_INVERSION:
                break;
            case SV_SPECIAL:
                if ((detectPars.detectSet & (1 << shiftSize)) != 0)
                {
                    if (detectPars.outputPrefix != NULL)
                    {
                        string outputFile(detectPars.outputPrefix);
                        outputFile += ".mei.vcf";
                        outputGrp.fpSpecial = fopen(outputFile.c_str(), "w");
                        if (outputGrp.fpSpecial == NULL)
                            TGM_ErrQuit("Error: Cannot open the MEI VCF file: %s\n", outputFile.c_str());
                    }

                    PrintSpecialHeader();
                }
                break;
            case SV_INTER_CHR_TRNSLCTN:
                break;
            default:
                break;
        }
    }
}

void Printer::PrintSpecialHeader(void)
{
    // get the current time
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%Y%m%d", &tstruct);

    // output to the standard input
    if (outputGrp.fpSpecial == NULL)
    {

        printf("##fileformat=VCFv4.1\n"
               "##fileDate=%s\n"
               "##source=Tangram\n"
               "##ALT=<ID=INS:ME:AL,Description=\"Insertion of ALU element\">\n"
               "##ALT=<ID=INS:ME:L1,Description=\"Insertion of L1 element\">\n"
               "##ALT=<ID=INS:ME:SV,Description=\"Insertion of SVA element\">\n"
               "##ALT=<ID=INS:ME:HE,Description=\"Insertion of HERVK element\">\n"
               "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">\n"
               "##INFO=<ID=STRAND,Number=1,Type=String,Description=\"Orientation of the inserted mobile elements. '/' means the strand information is not available\">\n"
               "##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants. Only presents if the 'IMPRECISE' flag is set\">\n"
               "##INFO=<ID=MEILEN,Number=1,Type=Integer,Description=\"Inserted length of MEI. -1 means the inserted length is not available.\">\n"
	       "##INFO=<ID=RP5,Number=1,Type=Integer,Description=\"supportive read-pair fragments from 5-prime\">\n"
	       "##INFO=<ID=RP3,Number=1,Type=Integer,Description=\"supportive read-pair fragments from 3-prime\">\n"
	       "##INFO=<ID=SR5,Number=1,Type=Integer,Description=\"supportive split-read fragments from 5-prime\">\n"
	       "##INFO=<ID=SR3,Number=1,Type=Integer,Description=\"supportive split-read fragments from 3-prime\">\n"
               "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
               "##FORMAT=<ID=GL,Number=3,Type=Float,Description=\"Genotype likelihood\">\n"
               "##FORMAT=<ID=RO,Number=1,Type=Integer,Description=\"Reference allele observation count\">\n"
	       "##FORMAT=<ID=R5,Number=1,Type=Integer,Description=\"supportive read-pair fragments from 5-prime\">\n"
	       "##FORMAT=<ID=R3,Number=1,Type=Integer,Description=\"supportive read-pair fragments from 3-prime\">\n"
	       "##FORMAT=<ID=S5,Number=1,Type=Integer,Description=\"supportive split-read fragments from 5-prime\">\n"
	       "##FORMAT=<ID=S3,Number=1,Type=Integer,Description=\"supportive split-read fragments from 3-prime\">\n",
               buf);

        printf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    }
    else // output to a file
    {

        fprintf(outputGrp.fpSpecial,
	       "##fileformat=VCFv4.1\n"
               "##fileDate=%s\n"
               "##source=Tangram\n"
               "##ALT=<ID=INS:ME:AL,Description=\"Insertion of ALU element\">\n"
               "##ALT=<ID=INS:ME:L1,Description=\"Insertion of L1 element\">\n"
               "##ALT=<ID=INS:ME:SV,Description=\"Insertion of SVA element\">\n"
               "##ALT=<ID=INS:ME:HE,Description=\"Insertion of HERVK element\">\n"
               "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">\n"
               "##INFO=<ID=STRAND,Number=1,Type=String,Description=\"Orientation of the inserted mobile elements. '/' means the strand information is not available\">\n"
               "##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants. Only presents if the 'IMPRECISE' flag is set\">\n"
               "##INFO=<ID=MEILEN,Number=1,Type=Integer,Description=\"Inserted length of MEI. -1 means the inserted length is not available.\">\n"
	       "##INFO=<ID=RP5,Number=1,Type=Integer,Description=\"supportive read-pair fragments from 5-prime\">\n"
	       "##INFO=<ID=RP3,Number=1,Type=Integer,Description=\"supportive read-pair fragments from 3-prime\">\n"
	       "##INFO=<ID=SR5,Number=1,Type=Integer,Description=\"supportive split-read fragments from 5-prime\">\n"
	       "##INFO=<ID=SR3,Number=1,Type=Integer,Description=\"supportive split-read fragments from 3-prime\">\n"
               "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
               "##FORMAT=<ID=GL,Number=3,Type=Float,Description=\"Genotype likelihood\">\n"
               "##FORMAT=<ID=RO,Number=1,Type=Integer,Description=\"Reference allele observation count\">\n"
	       "##FORMAT=<ID=R5,Number=1,Type=Integer,Description=\"supportive read-pair fragments from 5-prime\">\n"
	       "##FORMAT=<ID=R3,Number=1,Type=Integer,Description=\"supportive read-pair fragments from 3-prime\">\n"
	       "##FORMAT=<ID=S5,Number=1,Type=Integer,Description=\"supportive split-read fragments from 5-prime\">\n"
	       "##FORMAT=<ID=S3,Number=1,Type=Integer,Description=\"supportive split-read fragments from 3-prime\">\n",
               buf);

        fprintf(outputGrp.fpSpecial, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    }

    // print ou the sample names
    const Array<char*>& sampleNames = libTable.GetSampleNames();

    unsigned int numSamples = sampleNames.Size();
    for (unsigned int i = 0; i != numSamples; ++i)
    {
        if (outputGrp.fpSpecial == NULL)
            printf("\t%s", sampleNames[i]);
        else
            fprintf(outputGrp.fpSpecial, "\t%s", sampleNames[i]);

    }

    if (outputGrp.fpSpecial == NULL)
        printf("\n");
    else
        fprintf(outputGrp.fpSpecial, "\n");
}

void Printer::PrintSpecial(const PrintElmnt& element)
{
    InitFeatures();
    SetSpecialFeatures(element);

    int insertedLen = -1;
    string insertedSeq;
    if (element.pSplitEvent != NULL)
    {
        if (element.pSplitEvent->pSpecialData->end >= 0 && element.pSplitEvent->pSpecialData->pos >=0) {
            insertedLen = element.pSplitEvent->pSpecialData->end - element.pSplitEvent->pSpecialData->pos + 1;
	    /*
	    int concatenated_pos = pRef->GetSpRefConcatenatedPos(
	                             element.pSplitEvent->pSpecialData->spRefID, 
	                             element.pSplitEvent->pSpecialData->pos);
            const int8_t *sr_seq_ptr = pRef->spRefSeq.GetPointer(concatenated_pos);
	    if (sr_seq_ptr != NULL) {
              if (features.strand == '-')  
	        ReverseSeqToString(insertedSeq, sr_seq_ptr, insertedLen);
              else
	        SeqToString(insertedSeq, sr_seq_ptr, insertedLen);
	    }
	    */
	}
    }

    int splitFrag = features.splitFrag[0] + features.splitFrag[1];
    char refChar = char_table[pRef->refSeq[features.pos - pRef->pos]];

    formatted.clear();
    formatted.str("");

    if (splitFrag == 0)
    {
        int ciPos1 = 0;
        int ciPos2 = 0;
        if (features.pos5[1] < features.pos3[0])
        {
            ciPos1 = features.pos5[0] - features.pos5[1];
            ciPos2 = features.pos3[1] - features.pos5[1];
        }
        else
        {
            ciPos1 = features.pos5[0] - features.pos3[0];
            ciPos2 = features.pos3[1] - features.pos3[0];
        }

        formatted << "CIPOS=" << ciPos1 << "," << ciPos2 << ";";
    }

    if (outputGrp.fpSpecial == NULL)
    {
        string alt;
	if (insertedSeq.empty()) {
	  alt = "<INS:ME:";
	  alt += features.spRefName;
	  alt += ">";
	} else {
	  alt = refChar;
	  alt += insertedSeq;
	}
	printf("%s\t"
               "%d\t"
               ".\t"
               "%c\t"
               "%s\t"
               ".\t"
               ".\t"
               "%s%s",
               features.anchorName,
               features.pos + 1,
               refChar,
               alt.c_str(),
	       splitFrag > 0 ? "" : "IMPRECISE;",
               formatted.str().c_str()
              );
	#ifdef TD_VERBOSE_DEBUG
	if (element.pSplitEvent != NULL)
	  fprintf(stderr, "%s\t%d\t%s\t%d\t%d\t%d\t%d\n", features.anchorName, features.pos + 1, features.spRefName, 
	                element.pSplitEvent->pSpecialData->spRefID,
			element.pSplitEvent->pSpecialData->familyID,
			element.pSplitEvent->pSpecialData->pos,
			element.pSplitEvent->pSpecialData->end);
	#endif
    }
    else
    {
        string alt;
	if (alt.empty()) {
	  alt = "<INS:ME:";
	  alt += features.spRefName;
	  alt += ">";
	} else {
	  alt = refChar;
	  alt += insertedSeq;
	}
	fprintf(outputGrp.fpSpecial, "%s\t"
               "%d\t"
               ".\t"
               "%c\t"
               "%s\t"
               ".\t"
               ".\t"
               "%s%s",
               features.anchorName,
               features.pos + 1,
	       refChar,
	       alt.c_str(),
               splitFrag > 0 ? "" : "IMPRECISE;",
               formatted.str().c_str()
              );
	#ifdef TD_VERBOSE_DEBUG
	if (element.pSplitEvent != NULL)
	  fprintf(stderr, "%s\t%d\t%s\t%d\t%d\t%d\t%d\n", features.anchorName, features.pos + 1, features.spRefName, 
	                element.pSplitEvent->pSpecialData->spRefID,
			element.pSplitEvent->pSpecialData->familyID,
			element.pSplitEvent->pSpecialData->pos,
			element.pSplitEvent->pSpecialData->end);
	#endif
    }

    formatted.clear();
    formatted.str("");

    formatted << "RP5=" << features.rpFrag[0] << ";"
              << "RP3=" << features.rpFrag[1] << ";" 
	      << "SR5=" << features.splitFrag[0] << ";" 
	      << "SR3=" << features.splitFrag[1] << ";";

    if (outputGrp.fpSpecial == NULL)
    {
        printf("%sSTRAND=%c;MEILEN=%d",
               formatted.str().c_str(),
               features.strand,
               insertedLen
              );
    }
    else
    {
        fprintf(outputGrp.fpSpecial, "%sSTRAND=%c;MEILEN=%d",
                formatted.str().c_str(),
                features.strand,
                insertedLen
               );
    }
}

void Printer::SetSpecialFeatures(const PrintElmnt& element)
{
    if (element.pSplitEvent != NULL)
        SetSpecialFeaturesFromSplit(*(element.pSplitEvent));

    if (element.pRpSpecial != NULL)
        SetSpecialFeaturesFromRp(*(element.pRpSpecial));
}

void Printer::SetSpecialFeaturesFromSplit(const SplitEvent& splitEvent)
{
    const Array<char*>* pSpecialRefs = libTable.GetSpecialRefNames();

    features.splitFrag[0] = splitEvent.success3; 
    features.splitFrag[1] = splitEvent.success5;

    features.pos = splitEvent.pos;
    features.len = splitEvent.len;

    features.strand = splitEvent.strand == 0 ? '+' : '-';

    int zaID = pRef->familyToZA[splitEvent.pSpecialData->familyID];
    if (zaID < 0)
        features.spRefName = pRef->familyName[splitEvent.pSpecialData->familyID].c_str();
    else
        features.spRefName = (*pSpecialRefs)[zaID];

    features.anchorName = pRef->refHeader.names[splitEvent.refID];

    features.pos5[0] = splitEvent.pos5[0];
    features.pos5[1] = splitEvent.pos5[1];

    features.pos3[0] = splitEvent.pos3[0];
    features.pos3[1] = splitEvent.pos3[1];
}

void Printer::SetSpecialFeaturesFromRp(const SpecialEvent& rpSpecial)
{
    const Array<char*>& anchorNames = libTable.GetAnchorNames();
    const Array<char*>* pSpecialRefs = libTable.GetSpecialRefNames();

    features.rpFrag[0] = rpSpecial.numFrag[0];
    features.rpFrag[1] = rpSpecial.numFrag[1];

    if (features.anchorName == NULL)
    {
        features.anchorName = anchorNames[rpSpecial.refID];
        features.spRefName = (*pSpecialRefs)[rpSpecial.familyID];
        
        features.pos = rpSpecial.pos;
        features.len = rpSpecial.length;

        features.pos5[0] = rpSpecial.pos5[0];
        features.pos5[1] = rpSpecial.pos5[1];

        features.pos3[0] = rpSpecial.pos3[0];
        features.pos3[1] = rpSpecial.pos3[1];
    }
}

void Printer::PrintGenotype(const Genotype& genotype, bool hasGenotype)
{
    formatted.clear();
    formatted.str("");

    if (hasGenotype)
    {
        unsigned int size = genotype.genotypes.size();
        for (unsigned int i = 0; i != size; ++i)
        {
            switch(genotype.genotypes[i])
            {
                case -1:
                    formatted << "\t.:";
                    break;
                case 0:
                    formatted << "\t0/0:";
                    break;
                case 1:
                    formatted << "\t0/1:";
                    break;
                case 2:
                    formatted << "\t1/1:";
                    break;
                default:
                    break;
            }

            unsigned int idx = 3 * i;
            formatted << setiosflags(ios::fixed) << setprecision(2) << genotype.likelihoods[idx] << ",";
            formatted << setiosflags(ios::fixed) << setprecision(2) << genotype.likelihoods[idx + 1] << ",";
            formatted << setiosflags(ios::fixed) << setprecision(2) << genotype.likelihoods[idx + 2];

            formatted << ":" << genotype.sampleCount[i].nonSupport
                      << ":" << genotype.sampleCount[i].rp3 
		      << ":" << genotype.sampleCount[i].rp5
                      << ":" << genotype.sampleCount[i].sr3 
		      << ":" << genotype.sampleCount[i].sr5;
        }
    }
    else
    {
        const Array<char*>& sampleNames = libTable.GetSampleNames();

        unsigned int numSamples = sampleNames.Size();
        for (unsigned int i = 0; i != numSamples; ++i)
        {
            formatted << "\t.:.";

	    formatted << "0.00,0.00,0.00";

            formatted << ":" << genotype.sampleCount[i].nonSupport
                      << ":" << genotype.sampleCount[i].rp3 
		      << ":" << genotype.sampleCount[i].rp5
                      << ":" << genotype.sampleCount[i].sr3 
		      << ":" << genotype.sampleCount[i].sr5;
        }
    }

    if (outputGrp.fpSpecial == NULL)
    {
        printf("\tGT:GL:RO:R5:R3:S5:S3");
        printf("%s\n", formatted.str().c_str());
    }
    else
    {
        fprintf(outputGrp.fpSpecial, "\tGT:GL:RO:R5:R3:S5:S3");
        fprintf(outputGrp.fpSpecial, "%s\n", formatted.str().c_str());
    }
}

void Printer::PrintSampleInfo(const Genotype& genotype)
{
    formatted.clear();
    formatted.str("");

    unsigned int size = genotype.sampleCount.Size();
    for (unsigned int i = 0; i != size; ++i)
    {
        if (genotype.sampleCount[i].support > 0) {
            formatted << "\t./1";
	    
            formatted << ":" << genotype.sampleCount[i].nonSupport
                      << ":" << genotype.sampleCount[i].rp3 
		      << ":" << genotype.sampleCount[i].rp5
                      << ":" << genotype.sampleCount[i].sr3 
		      << ":" << genotype.sampleCount[i].sr5;
        } else {
            formatted << "\t0/0";

            formatted << ":" << genotype.sampleCount[i].nonSupport
                      << ":" << genotype.sampleCount[i].rp3 
		      << ":" << genotype.sampleCount[i].rp5
                      << ":" << genotype.sampleCount[i].sr3 
		      << ":" << genotype.sampleCount[i].sr5;
	}
    }

    if (outputGrp.fpSpecial == NULL)
    {
        printf("\tGT:GL:RO:R5:R3:S5:S3");
        printf("%s\n", formatted.str().c_str());
    }
    else
    {
        fprintf(outputGrp.fpSpecial, "\tGT:GL:RO:R5:R3:S5:S3");
        fprintf(outputGrp.fpSpecial, "%s\n", formatted.str().c_str());
    }
}
