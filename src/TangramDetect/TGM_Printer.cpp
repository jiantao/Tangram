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
#include "TGM_Printer.h"

using namespace std;
using namespace Tangram;


Printer::Printer(const Detector* pDetector, const DetectPars& detectPars, const Aligner* pAligner, const Reference* pRef, const LibTable& libTable, const BamPairTable& bamPairTable)
                : pDetector(pDetector), detectPars(detectPars), pAligner(pAligner), pRef(pRef), libTable(libTable), bamPairTable(bamPairTable)
{
    InitFeatures();
}

Printer::~Printer()
{

}

void Printer::Print(void)
{
    for (unsigned int i = 0; i != NUM_SV_TYPES; ++i)
    {
        switch (i)
        {
            case SV_SPECIAL:
                PrintSpecial();
                break;
            default:
                break;
        }
    }
}

void Printer::PrintSpecial(void)
{
    FILE* fpOutput = NULL;
    if (detectPars.outputPrefix != NULL)
    {
        string outputFile(detectPars.outputPrefix);
        outputFile += ".mei.vcf";
        fpOutput = fopen(outputFile.c_str(), "w");
        if (fpOutput == NULL)
            TGM_ErrQuit("Error: Cannot open the MEI VCF file: %s\n", outputFile.c_str());
    }

    unsigned int numSp = libTable.GetNumSpecialRef();
    unsigned int numSamples = libTable.GetNumSamples();

    PrintSpecialHeader(fpOutput);
    sampleMap.Init(numSamples);
    sampleMap.SetSize(numSamples);

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

        InitPrintIdx(numEvents, splitLen);

        while (GetNextSpecial(pRpSpecials, pSplitSpecials))
        {
            sampleMap.MemSet(0);
            InitFeatures();

            if (printIdx.pRpSpecial != NULL)
                SetSampleInfoRpSpecial(*(printIdx.pRpSpecial));

            if (printIdx.pSplitEvent != NULL)
                SetSampleInfoSplit(*(printIdx.pSplitEvent));

            // SetSampleString();
            SetSpecialFeatures(i);
            PrintSpecialBody(fpOutput);

            /*  
            printf("chr%s\t%d\t%d\t%c\t%d\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%s\n", features.anchorName, features.pos, features.pos + features.len + 1, features.strand, 
                    features.rpFrag[0], features.rpFrag[1], features.splitFrag[0], features.splitFrag[1], features.spRefName, features.pos5[0], features.pos5[1], 
                    features.pos3[0], features.pos3[1], formatted.str().c_str());
            */
        }
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

                for (unsigned int j = 0; j != len; ++j)
                {
                    sampleMap.MemSet(0);
                    InitFeatures();

                    SetSampleInfoSplit(pSplitSpecials[j]);

                    // SetSampleString();
                    SetSpecialFeaturesFromSplit(pSplitSpecials[j]);

                    PrintSpecialBody(fpOutput);

                    /*
                    printf("chr%s\t%d\t%d\t%c\t%d\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%s\n", features.anchorName, features.pos, features.pos + features.len + 1, features.strand, 
                            features.rpFrag[0], features.rpFrag[1], features.splitFrag[0], features.splitFrag[1], features.spRefName, features.pos5[0], features.pos5[1], 
                            features.pos3[0], features.pos3[1], formatted.str().c_str());
                    */
                }
            }
        }
    }
}

void Printer::PrintSpecialHeader(FILE* fpOutput)
{
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%Y%m%d", &tstruct);

    if (fpOutput == NULL)
    {
        printf("##fileformat=VCFv4.1\n"
               "##fileDate=%s\n"
               "##source=Tangram\n"
               "##ALT=<ID=INS:ME:AL,Description=\"Insertion of ALU element\">\n"
               "##ALT=<ID=INS:ME:L1,Description=\"Insertion of L1 element\">\n"
               "##ALT=<ID=INS:ME:SV,Description=\"Insertion of SVA element\">\n"
               "##ALT=<ID=INS:ME:HE,Description=\"Insertion of HERV element\">\n"
               "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">\n"
               "##INFO=<ID=STRAND,Number=1,Type=String,Description=\"Orientation of the inserted mobile elements. '/' means the strand information is not available\">\n"
               "##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants. Only presents if the 'IMPRECISE' flag is set\">\n"
               "##INFO=<ID=MEILEN,Number=1,Type=Integer,Description=\"Inserted length of MEI. -1 means the inserted length is not available.\">\n"
               "##INFO=<ID=FRAG,Number=4,Type=Integer,Description=\"Detailed information of supporting fragments: 5' read-pair fragments, 3' read-pair fragments,"
               "5' split fragments and 3' split fragments\">\n"
               "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
               "##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Allele Depth, how many reads support this allele\">\n",
               buf);

        printf("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    }
    else
    {
        fprintf(fpOutput, 
               "##fileformat=VCFv4.1\n"
               "##fileDate=%s\n"
               "##source=Tangram\n"
               "##ALT=<ID=INS:ME:AL,Description=\"Insertion of ALU element\">\n"
               "##ALT=<ID=INS:ME:L1,Description=\"Insertion of L1 element\">\n"
               "##ALT=<ID=INS:ME:SV,Description=\"Insertion of SVA element\">\n"
               "##ALT=<ID=INS:ME:HE,Description=\"Insertion of HERV element\">\n"
               "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">\n"
               "##INFO=<ID=STRAND,Number=1,Type=String,Description=\"Orientation of the inserted mobile elements. '/' means the strand information is not available\">\n"
               "##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants. Only presents if the 'IMPRECISE' flag is set\">\n"
               "##INFO=<ID=MEILEN,Number=1,Type=Integer,Description=\"Inserted length of MEI. -1 means the inserted length is not available.\">\n"
               "##INFO=<ID=FRAG,Number=4,Type=Integer,Description=\"Detailed information of supporting fragments: 5' read-pair fragments, 3' read-pair fragments,"
               "5' split fragments and 3' split fragments\">\n"
               "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
               "##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Allele Depth, how many reads support this allele\">\n",
               buf);

        fprintf(fpOutput, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    }

    const Array<char*>& sampleNames = libTable.GetSampleNames();

    unsigned int numSamples = sampleNames.Size();
    for (unsigned int i = 0; i != numSamples; ++i)
    {
        if (fpOutput == NULL)
            printf("\t%s", sampleNames[i]);
        else
            fprintf(fpOutput, "\t%s", sampleNames[i]);

    }

    if (fpOutput == NULL)
        printf("\n");
    else
        fprintf(fpOutput, "\n");
}

void Printer::PrintSpecialBody(FILE* fpOutput)
{
    int insertedLen = -1;
    if (printIdx.pSplitEvent != NULL)
    {
        if (printIdx.pSplitEvent->pSpecialData->end >= 0 && printIdx.pSplitEvent->pSpecialData->pos >=0)
            insertedLen = printIdx.pSplitEvent->pSpecialData->end - printIdx.pSplitEvent->pSpecialData->pos + 1;
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

    if (fpOutput == NULL)
    {
        printf("chr%s\t"
               "%d\t"
               ".\t"
               "%c\t"
               "<INS:ME:%s>\t"
               ".\t"
               ".\t"
               "%s%s",
               features.anchorName,
               features.pos,
               refChar,
               features.spRefName,
               splitFrag > 0 ? "" : "IMPRECISE;",
               formatted.str().c_str()
              );
    }
    else
    {
        fprintf(fpOutput, "chr%s\t"
               "%d\t"
               ".\t"
               "%c\t"
               "<INS:ME:%s>\t"
               ".\t"
               ".\t"
               "%s%s",
               features.anchorName,
               features.pos,
               refChar,
               features.spRefName,
               splitFrag > 0 ? "" : "IMPRECISE;",
               formatted.str().c_str()
              );
    }

    formatted.clear();
    formatted.str("");

    formatted << "FRAG=" << features.rpFrag[0] << "," << features.rpFrag[1] << "," << features.splitFrag[0] << "," << features.splitFrag[1] << ";";

    if (fpOutput == NULL)
    {
        printf("%sSTRAND=%c;MEILEN=%d\t"
               "GT:AD",
               formatted.str().c_str(),
               features.strand,
               insertedLen
              );
    }
    else
    {
        fprintf(fpOutput, "%sSTRAND=%c;MEILEN=%d\t"
               "GT:AD",
               formatted.str().c_str(),
               features.strand,
               insertedLen
              );
    }

    unsigned int numSamples = libTable.GetNumSamples();
    char buff[4];
    buff[1] = '/';
    buff[3] = '\0';

    for (unsigned int i = 0; i != numSamples; ++i)
    {
        if (sampleMap[i] > 0)
        {
            buff[0] = '1';
            buff[2] = '.';
        }
        else
        {
            buff[0] = '0';
            buff[2] = '0';
        }

        if (fpOutput == NULL)
            printf("\t%s:%d", buff, sampleMap[i]);
        else
            fprintf(fpOutput, "\t%s:%d", buff, sampleMap[i]);
    }

    if (fpOutput == NULL)
        printf("\n");
    else
        fprintf(fpOutput, "\n");
}

bool Printer::GetNextSpecial(const SpecialEvent* pRpSpecials, const SplitEvent* pSplitSpecials)
{
    if (printIdx.rpIdx == printIdx.rpSize)
        printIdx.pRpSpecial = NULL;
    else
        printIdx.pRpSpecial = pRpSpecials + printIdx.rpIdx;

    if (printIdx.pRpSpecial != NULL)
    {
        printIdx.pRpSpecial = NULL;
        while (printIdx.rpIdx < printIdx.rpSize)
        {
            if (pRpSpecials[printIdx.rpIdx].splitIdx >= 0)
                ++(printIdx.rpIdx);
            else
            {
                printIdx.pRpSpecial = pRpSpecials + printIdx.rpIdx;
                break;
            }
        }
    }

    if (printIdx.splitIdx == printIdx.splitSize)
        printIdx.pSplitEvent = NULL;
    else
        printIdx.pSplitEvent = pSplitSpecials + printIdx.splitIdx;

    if (printIdx.pRpSpecial == NULL && printIdx.pSplitEvent == NULL)
        return false;
    if (printIdx.pRpSpecial != NULL && printIdx.pSplitEvent != NULL)
    {
        if (printIdx.pRpSpecial->pos <= (unsigned int) printIdx.pSplitEvent->pos)
            printIdx.pSplitEvent = NULL;
        else
            printIdx.pRpSpecial = NULL;
    }

    if (printIdx.pRpSpecial != NULL)
        ++(printIdx.rpIdx);

    if (printIdx.pSplitEvent != NULL)
    {
        int rpIdx = printIdx.pSplitEvent->rpIdx;
        if (rpIdx >= 0)
        {
            printIdx.pRpSpecial = pRpSpecials + rpIdx;
        }

        ++(printIdx.splitIdx);
    }

    return true;
}

void Printer::SetSampleInfoRpSpecial(const SpecialEvent& rpSpecial)
{
    for (unsigned int k = 0; k != rpSpecial.numFrag[0]; ++k)
    {
        const SpecialPair& specialPair = bamPairTable.specialPairs[rpSpecial.origIndex[0][k]];
        unsigned int readGrpID = specialPair.readGrpID;
        unsigned int sampleID = 0;

        if (libTable.GetSampleID(sampleID, readGrpID))
            ++sampleMap[sampleID];
    }

    for (unsigned int k = 0; k != rpSpecial.numFrag[1]; ++k)
    {
        const SpecialPair& specialPair = bamPairTable.specialPairs[rpSpecial.origIndex[1][k]];
        unsigned int readGrpID = specialPair.readGrpID;
        unsigned int sampleID = 0;

        if (libTable.GetSampleID(sampleID, readGrpID))
            ++sampleMap[sampleID];
    }
}

void Printer::SetSampleInfoSplit(const SplitEvent& splitEvent)
{
    for (unsigned int i = 0; i != splitEvent.size3; ++i)
    {
        if (splitEvent.first3[i].isMajor)
        {
            unsigned int readGrpID = 0;
            unsigned int sampleID = 0;
            unsigned int origIdx = splitEvent.first3[i].origIdx;

            if (splitEvent.first3[i].isSoft)
                readGrpID = bamPairTable.softPairs[origIdx].readGrpID;
            else
                readGrpID = bamPairTable.orphanPairs[origIdx].readGrpID;

            if (libTable.GetSampleID(sampleID, readGrpID))
                ++sampleMap[sampleID];

            ++(features.splitFrag[0]);
        }
    }

    for (unsigned int i = 0; i != splitEvent.size5; ++i)
    {
        if (splitEvent.first5[i].isMajor)
        {
            unsigned int readGrpID = 0;
            unsigned int sampleID = 0;
            unsigned int origIdx = splitEvent.first5[i].origIdx;

            if (splitEvent.first5[i].isSoft)
                readGrpID = bamPairTable.softPairs[origIdx].readGrpID;
            else
                readGrpID = bamPairTable.orphanPairs[origIdx].readGrpID;

            if (libTable.GetSampleID(sampleID, readGrpID))
                ++sampleMap[sampleID];

            ++(features.splitFrag[1]);
        }
    }
}

/*
void Printer::SetSampleString(void)
{
    const Array<char*>& sampleNames = libTable.GetSampleNames();

    formatted.clear();
    formatted.str("");

    for (std::map<unsigned int, unsigned int>::const_iterator itor = sampleMap.begin(); itor != sampleMap.end(); ++itor)
    {
        unsigned int sampleID = itor->first;
        const char* name = sampleNames[sampleID];
        formatted << name << "_" << itor->second << "_";
    }
}
*/

void Printer::SetSpecialFeatures(unsigned int zaSpRefID)
{
    if (printIdx.pSplitEvent != NULL)
    {
        const SplitEvent& splitEvent = *(printIdx.pSplitEvent);
        SetSpecialFeaturesFromSplit(splitEvent);
    }

    if (printIdx.pRpSpecial != NULL)
    {
        const SpecialEvent& rpSpecial = *(printIdx.pRpSpecial);
        SetSpecialFeaturesFromRp(rpSpecial, zaSpRefID);
    }
}

void Printer::SetSpecialFeaturesFromSplit(const SplitEvent& splitEvent)
{
    const Array<char*>* pSpecialRefs = libTable.GetSpecialRefNames();

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

void Printer::SetSpecialFeaturesFromRp(const SpecialEvent& rpSpecial, unsigned int zaSpRefID)
{
    const Array<char*>& anchorNames = libTable.GetAnchorNames();
    const Array<char*>* pSpecialRefs = libTable.GetSpecialRefNames();

    features.rpFrag[0] = rpSpecial.numFrag[0];
    features.rpFrag[1] = rpSpecial.numFrag[1];

    if (features.anchorName == NULL)
    {
        features.anchorName = anchorNames[rpSpecial.refID];
        features.spRefName = (*pSpecialRefs)[zaSpRefID];
        
        features.pos = rpSpecial.pos;
        features.len = rpSpecial.length;

        features.pos5[0] = rpSpecial.pos5[0];
        features.pos5[1] = rpSpecial.pos5[1];

        features.pos3[0] = rpSpecial.pos3[0];
        features.pos3[1] = rpSpecial.pos3[1];
    }
}
