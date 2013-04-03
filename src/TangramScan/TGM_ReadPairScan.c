/*
 * =====================================================================================
 *
 *       Filename:  TGM_ReadPairScan.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  04/11/2012 02:44:36 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include "khash.h"
#include "TGM_Error.h"
#include "TGM_Utilities.h"
#include "TGM_ReadPairScan.h"
#include "TGM_BamPairAux.h"

static const char* TGM_LibTableFileName = "lib_table.dat";

static const char* TGM_HistFileName = "hist.dat";

KHASH_MAP_INIT_STR(name, uint32_t);

void TGM_ReadPairScan(const TGM_ReadPairScanPars* pScanPars)
{
    // some default capacity of the containers
    unsigned int capAnchor = 150;
    unsigned int capSample = 20;
    unsigned int capReadGrp = 20;
    unsigned int capHist = 10;


    // initialize the library information table
    TGM_LibInfoTable* pLibTable = TGM_LibInfoTableAlloc(capAnchor, capSample, capReadGrp, pScanPars->specialPrefix);
    TGM_LibInfoTableSetCutoff(pLibTable, pScanPars->cutoff);
    TGM_LibInfoTableSetTrimRate(pLibTable, pScanPars->trimRate);

    TGM_Status status = TGM_CheckWorkingDir(pScanPars->workingDir);
    if (status != TGM_OK)
        TGM_ErrQuit("ERROR: Error found during creating the working directory.\n");

    // get the library table output file name
    char* libTableOutputFile = TGM_CreateFileName(pScanPars->workingDir, TGM_LibTableFileName);

    // write the library information table into the file
    FILE* libTableOutput = fopen(libTableOutputFile, "wb");
    if (libTableOutput == NULL)
        TGM_ErrQuit("ERROR: Cannot open the library output file: %s\n", libTableOutputFile);

    free(libTableOutputFile);

    // structure initialization
    TGM_BamInStreamLite* pBamInStreamLite = TGM_BamInStreamLiteAlloc();

    TGM_SpecialID* pSpecialID = TGM_SpecialIDAlloc(10);

    // buffer used to hold the bam file name
    char bamFileName[TGM_MAX_LINE];

    TGM_FragLenHistArray* pHistArray = TGM_FragLenHistArrayAlloc(capHist);
    TGM_BamHeader* pBamHeader = NULL;

    // open the fragment length histogram output file
    char* histOutputFile = TGM_CreateFileName(pScanPars->workingDir, TGM_HistFileName);
    FILE* histOutput = fopen(histOutputFile, "w");
    if (histOutput == NULL)
        TGM_ErrQuit("ERROR: Cannot open fragment length histogram file: %s\n", histOutputFile);

    free(histOutputFile);

    uint32_t readGrpCount = 0;
    TGM_FragLenHistArrayWriteHeader(readGrpCount, histOutput);

    while (TGM_GetNextLine(bamFileName, TGM_MAX_LINE, pScanPars->fileListInput) == TGM_OK)
    {
        // open the bam file
        TGM_BamInStreamLiteOpen(pBamInStreamLite, bamFileName);

        // load the bam header before read any alignments
        pBamHeader = TGM_BamInStreamLiteLoadHeader(pBamInStreamLite);

        // process the header information
        unsigned int oldSize = 0;
        if (TGM_LibInfoTableSetRG(pLibTable, &oldSize, pBamHeader) != TGM_OK)
            TGM_ErrQuit("ERROR: Found an error when loading the bam file.\n");

        // initialize the fragment length histogram array with the number of newly added libraries in the bam file
        TGM_FragLenHistArrayInit(pHistArray, pLibTable->size - oldSize);

        // do not load cross pairs when build the fragment length distribution
        TGM_Bool loadCross = FALSE;

        // filter data for the no-za sort mode
        TGM_FilterDataNoZA filterData = {pLibTable, TRUE};

        // get the sorting order from the bam header
        TGM_SortMode sortMode = TGM_BamHeaderGetSortMode(pBamHeader);
        if (sortMode == TGM_SORTED_COORDINATE_NO_ZA || sortMode == TGM_SORTED_NAME || sortMode == TGM_SORTED_SPLIT)
        {
            // test if the bam file has za tag or not
            TGM_Status status = TGM_BamInStreamLiteTestZA(pBamInStreamLite);

            if (status == TGM_OK)
                sortMode = TGM_SORTED_COORDINATE_ZA;
            else if (status == TGM_ERR)
                continue;

            TGM_BamInStreamLiteSetSortMode(pBamInStreamLite, sortMode);
        }
        else
            TGM_ErrQuit("ERROR: Invalid sorting order.\n");

        // set the sort order for the bam instream
        if (sortMode != TGM_SORTED_COORDINATE_NO_ZA)
        {
            TGM_BamInStreamLiteSetFilter(pBamInStreamLite, TGM_ReadPairFilter);
            TGM_BamInStreamLiteSetFilterData(pBamInStreamLite, &loadCross);
        }
        else
        {
            TGM_BamInStreamLiteSetFilter(pBamInStreamLite, TGM_ReadPairNoZAFilter);
            TGM_BamInStreamLiteSetFilterData(pBamInStreamLite, &filterData);
        }

        int retNum = 0;
        const bam1_t* pAlgns[3] = {NULL, NULL, NULL};

        TGM_Status bamStatus = TGM_OK;

        // read the primary bam the first time to build fragment length distribution
        do
        {
            int64_t index = -1;
            bamStatus = TGM_BamInStreamLiteRead(pAlgns, &retNum, &index, pBamInStreamLite);

            if (retNum > 0)
            {
                // check if the incoming read pair is normal (unique-unique pair)
                // if yes, then update the corresponding fragment length histogram
                TGM_PairStats pairStats;
                TGM_ZAtag zaTag;

                const TGM_MateInfo* pMateInfo = TGM_BamInStreamLiteGetMateInfo(pBamInStreamLite, index);

                unsigned int backHistIndex = 0;
                TGM_Status zaStatus = TGM_ERR;
                if (TGM_IsNormalPair(&pairStats, &zaTag, &zaStatus, pMateInfo, &backHistIndex, pAlgns, retNum, pLibTable, pScanPars->minMQ))
                    TGM_FragLenHistArrayUpdate(pHistArray, backHistIndex, pairStats.fragLen);

                if (zaStatus == TGM_OK)
                    TGM_SpecialIDUpdate(pSpecialID, &zaTag);
            }

        }while(bamStatus == TGM_OK);

        // finish the process of the histogram and update the library information table
        TGM_FragLenHistArrayFinalize(pHistArray);
        TGM_LibInfoTableUpdate(pLibTable, pHistArray, oldSize, pScanPars->checkLib);

        // write the fragment length histogram into the file
        TGM_FragLenHistArrayWrite(pHistArray, histOutput);

        // close the bam file
        TGM_BamInStreamLiteClose(pBamInStreamLite);
        TGM_BamHeaderFree(pBamHeader);
    }

    // write the library information into file
    TGM_FragLenHistArrayWriteHeader(pLibTable->size, histOutput);
    TGM_LibInfoTableWrite(pLibTable, TRUE, libTableOutput);
    TGM_SpecialIDWrite(pSpecialID, libTableOutput);

    // clean up
    fclose(libTableOutput);
    fclose(histOutput);
    TGM_SpecialIDFree(pSpecialID);
    TGM_LibInfoTableFree(pLibTable);
    TGM_FragLenHistArrayFree(pHistArray);
    TGM_BamInStreamLiteFree(pBamInStreamLite);
}

void TGM_SpecialIDUpdate(TGM_SpecialID* pSpecialID, const TGM_ZAtag* pZAtag)
{
    if (pSpecialID->size == pSpecialID->capacity)
    {
        pSpecialID->capacity *= 2;
        pSpecialID->names = (char (*)[3]) realloc(pSpecialID->names, sizeof(char) * 3 * pSpecialID->capacity);
        if (pSpecialID->names == NULL)
            TGM_ErrQuit("ERROR: Not enough memory for the special ID names.\n");
    }

    unsigned int i = pSpecialID->size;
    if (pZAtag->spRef[0][0] != ' ' && pZAtag->spRef[1][0] == ' ')
    {
        pSpecialID->names[i][0] = pZAtag->spRef[0][0];
        pSpecialID->names[i][1] = pZAtag->spRef[0][1];
    }
    else if (pZAtag->spRef[0][0] == ' ' && pZAtag->spRef[1][0] != ' ')
    {
        pSpecialID->names[i][0] = pZAtag->spRef[1][0];
        pSpecialID->names[i][1] = pZAtag->spRef[1][1];
    }
    else
        return;

    pSpecialID->names[i][2] = '\0';
    
    khiter_t khIter = 0;
    int ret = 0;

    khIter = kh_put(name, pSpecialID->pHash, pSpecialID->names[i], &ret);
    if (ret != 0)
    {
        kh_value((khash_t(name)*) pSpecialID->pHash, khIter) = i;
        ++(pSpecialID->size);
    }
}
