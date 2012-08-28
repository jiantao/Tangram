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

#include "api/BamMultiReader.h"

using namespace std;
using namespace Tangram;
using namespace BamTools;

int main(int argc, char *argv[])
{
    DetectPars detectPars;
    AlignerPars alignerPars;
    Parameters parameters(detectPars, alignerPars);
    parameters.Set((const char**) argv, argc);

    vector<string> filenames;
    parameters.SetBamFilenames(filenames);

    BamMultiReader bamMultiReader;
    bamMultiReader.Open(filenames);
    bamMultiReader.LocateIndexes();

    LibTable libTable;
    libTable.Read(detectPars.fpLibInput, detectPars.maxFragDiff);

    FragLenTable fragLenTable;
    fragLenTable.Read(detectPars.fpHistInput);

    parameters.ParseRangeStr(bamMultiReader);
    parameters.SetRange(bamMultiReader, libTable.GetFragLenMax());

    BamPairTable bamPairTable(libTable, fragLenTable, detectPars.minSoftSize, detectPars.minMQ, detectPars.spMinMQ);

    BamAlignment alignment;
    while(bamMultiReader.GetNextAlignment(alignment))
    {
        bamPairTable.Update(alignment);
    }

    fragLenTable.Destory();

    Detector detector;
    detector.Init(&detectPars, &libTable, &bamPairTable);

    detector.CallEvents();

    const Reference* pRef = NULL;
    const Aligner* pAligner = NULL;

    Reference reference;
    Aligner aligner(detector, bamPairTable,  alignerPars, reference, libTable);

    if (alignerPars.fpRefInput != NULL)
    {
        reference.Read(alignerPars.fpRefInput, detectPars.refID, detectPars.range[0], detectPars.range[1]);
        reference.CreateFamilyToZA(*(libTable.GetSpecialRefNames()));

        aligner.Map();

        pRef = &reference;
        pAligner = &aligner;
    }

    Printer printer(&detector, pAligner, pRef, libTable, bamPairTable);
    printer.Print();

    bamMultiReader.Close();
    return EXIT_SUCCESS;
}

