/*
 * =====================================================================================
 *
 *       Filename:  TGM_Tangram.cpp
 *
 *    Description:  
 *
 *        Created:  08/22/2012 08:38:37 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *   Inistitution:  Boston College
 *
 * =====================================================================================
 */

#include <cstdio>
#include <cstdlib>

#include "TGM_Array.h"
#include "TGM_Parameters.h"
#include "TGM_FragLenTable.h"
#include "TGM_LibTable.h"
#include "TGM_BamPair.h"
#include "TGM_Detector.h"
#include "TGM_Reference.h"
#include "TGM_Aligner.h"
#include "TGM_Printer.h"
#include "TGM_Genotype.h"

#include "api/BamMultiReader.h"

using namespace std;
using namespace Tangram;
using namespace BamTools;

int main(int argc, char *argv[])
{
    // load the command line arguments
    DetectPars detectPars;
    AlignerPars alignerPars;
    GenotypePars genotypePars;
    Parameters parameters(detectPars, alignerPars, genotypePars);
    parameters.Set((const char**) argv, argc);

    // load the input bam file names
    vector<string> filenames;
    parameters.SetBamFilenames(filenames);

    // open the bam files with bam multi reader
    BamMultiReader bamMultiReader;
    bamMultiReader.Open(filenames);

    // make sure we have the index files at with the input bam files
    if (!bamMultiReader.LocateIndexes())
    {
        // fail to locate the index files
        // we need to create them if necessary
        if (!bamMultiReader.CreateIndexes()) 
            TGM_ErrQuit("ERROR: Cannot locate or create the index files for the input bam files.\n");
    }

    // read the library information table
    LibTable libTable;
    libTable.Read(detectPars.fpLibInput, detectPars.maxFragDiff);

    // read the fragment length distribution table
    FragLenTable fragLenTable;
    fragLenTable.Read(detectPars.fpHistInput);

    // where should we start to call the SV events
    parameters.ParseRangeStr(bamMultiReader);
    parameters.SetRange(bamMultiReader, libTable.GetFragLenMax());

    BamPairTable bamPairTable(detectPars, libTable, fragLenTable);

    // iterate through the bam files and fill the bam pair table
    BamAlignment alignment;
    while(bamMultiReader.GetNextAlignment(alignment))
    {
        bamPairTable.Update(alignment);
    }

    // clean the fragnment length table
    fragLenTable.Destory();

    // call the SV events with read-pair signal
    Detector detector(detectPars, libTable, bamPairTable);
    detector.Init();
    detector.CallEvents();

    const Reference* pRef = NULL;
    const Aligner* pAligner = NULL;

    Reference reference;
    Aligner aligner(detector, bamPairTable,  alignerPars, reference, libTable);

    // call the SV events with split-read signal
    if (alignerPars.fpRefInput != NULL)
    {
        // always load the whole chromosome reference
        reference.Read(alignerPars.fpRefInput, detectPars.refID, -1, -1);
        reference.CreateFamilyToZA(*(libTable.GetSpecialRefNames()));

        aligner.Map();

        pRef = &reference;
        pAligner = &aligner;
    }

    // print out the events (vcf format)
    Genotype genotype(bamMultiReader, genotypePars, libTable, bamPairTable);
    genotype.Init();

    Printer printer(&detector, detectPars, pAligner, pRef, libTable, bamPairTable, genotypePars, genotype);
    printer.Init();
    printer.Print();

    // clean up and quit
    bamMultiReader.Close();
    return EXIT_SUCCESS;
}

