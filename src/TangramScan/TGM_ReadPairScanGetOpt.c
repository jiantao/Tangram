/*
 * =====================================================================================
 *
 *       Filename:  TGM_ReadPairScanGetOpt.c
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  04/17/2012 12:44:34 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (),
 *        Company:
 *
 * =====================================================================================
 */

#include <stdlib.h>
#include <string.h>

#include "TGM_Error.h"
#include "TGM_GetOpt.h"
#include "TGM_ReadPairScanGetOpt.h"

// total number of arguments we should expect for the split-read build program
#define OPT_SCAN_TOTAL_NUM 7

// total number of required arguments we should expect for the split-read build program
#define OPT_SCAN_REQUIRED_NUM 2

// the index of show help in the option object array
#define OPT_HELP           0

// the index of the fasta input file in the option object array
#define OPT_WORKING_DIR    1

// the index of the reference output file in the option object array
#define OPT_INPUT_FILE     2

// the index of the output hash table file in the option object array
#define OPT_CUTOFF         3

// the index of the hash size in the option object array
#define OPT_TRIM_RATE      4

#define OPT_MIN_MQ         5

#define OPT_SPECIAL_PREFIX 6

#define DEFAULT_SCAN_CUTOFF 0.01

#define DEFAULT_SCAN_TRIM_RATE 0.002

#define DEFAULT_SCAN_MIN_MQ 20

// set the parameters for the split-read build program from the pScanParsed command line arguments
void TGM_ReadPairScanSetPars(TGM_ReadPairScanPars* pScanPars, int argc, char* argv[])
{
    TGM_Option opts[] =
    {
        {"help", NULL, FALSE},
        {"dir",   NULL, FALSE},
        {"in",   NULL, FALSE},
        {"cf",   NULL, FALSE},
        {"tr",   NULL, FALSE},
        {"mq",  NULL, FALSE},
        {"sp",  NULL, FALSE},
        {NULL,   NULL, FALSE}
    };

    int optNum = TGM_GetOpt(opts, argc, argv);
    if (optNum < OPT_SCAN_REQUIRED_NUM)
        TGM_ReadPairScanHelp();

    for (unsigned int i = 0; i != OPT_SCAN_TOTAL_NUM; ++i)
    {
        switch (i)
        {
            case OPT_HELP:
                if (opts[i].isFound)
                    TGM_ReadPairScanHelp();
                break;
            case OPT_WORKING_DIR:
                if (opts[i].value == NULL)
                    TGM_ErrQuit("ERROR: The working directory is not specified.\n");

                unsigned int dirLen = strlen(opts[i].value);
                if (opts[i].value[dirLen - 1] == '/')
                    pScanPars->workingDir = strdup(opts[i].value);
                else
                {
                    pScanPars->workingDir = (char*) malloc(dirLen + 2);
                    if (pScanPars->workingDir == NULL)
                        TGM_ErrQuit("ERROR: Not enough memory for the working dir name.\n");

                    sprintf(pScanPars->workingDir, "%s/", opts[i].value);
                }

                break;
            case OPT_INPUT_FILE:
                if (opts[i].value == NULL)
                    TGM_ErrQuit("ERROR: The input bam list file is not specified.\n");

                pScanPars->fileListInput = fopen(opts[i].value, "r");
                if (pScanPars->fileListInput == NULL)
                    TGM_ErrQuit("ERROR: Cannot open file \"%s\" for read.\n", opts[i].value);

                break;
            case OPT_CUTOFF:
                if (opts[i].value == NULL)
                {
                    pScanPars->cutoff = DEFAULT_SCAN_CUTOFF;
                }
                else
                {
                    pScanPars->cutoff = atof(opts[i].value);
                    if (pScanPars->cutoff <= 0.0 || pScanPars->cutoff >= 1.0)
                        TGM_ErrQuit("ERROR: %s is an invalid cutoff value.\n", opts[i].value);
                }

                break;
            case OPT_TRIM_RATE:
                if (opts[i].value == NULL)
                {
                    pScanPars->trimRate = DEFAULT_SCAN_TRIM_RATE;
                }
                else
                {
                    pScanPars->trimRate = atof(opts[i].value);
                    if (pScanPars->trimRate < 0.0 || pScanPars->trimRate >= 1.0)
                        TGM_ErrQuit("ERROR: %s is an invalid trim rate.\n", opts[i].value);
                }

                break;
            case OPT_MIN_MQ:
                if (opts[i].value == NULL)
                {
                    pScanPars->minMQ = DEFAULT_SCAN_MIN_MQ;
                }
                else
                {
                    pScanPars->minMQ = atoi(opts[i].value);
                    if (pScanPars->minMQ < 0 || pScanPars->minMQ > 255)
                        TGM_ErrQuit("ERROR: %s is an invalid minimum mapping quality.\n", opts[i].value);
                }

                break;
            case OPT_SPECIAL_PREFIX:
                if (opts[i].isFound)
                {
                    if (opts[i].value != NULL)
                        pScanPars->specialPrefix = strdup(opts[i].value);
                    else
                        TGM_ErrQuit("ERROR: Special reference prefix is not specified.\n");

                    pScanPars->prefixLen = strlen(pScanPars->specialPrefix);
                }
                else
                {
                    pScanPars->specialPrefix = NULL;
                    pScanPars->prefixLen = 0;
                }

                break;
            default:
                TGM_ErrQuit("ERROR: Unrecognized argument.\n");
                break;
        }
    }
}

void TGM_ReadPairScanHelp(void)
{
    printf("Usage: tangram_scan [options] -in <input_file_list> -dir <output_dir>\n\n");

    printf("Mandatory arguments: -in   FILE   the list of input bam files\n");
    printf("                     -dir  STRING the path to the output dir (must be empty)\n\n");

    printf("Options:             -cf   FLOAT  threashold for normal read pair in the fragment length distribution[0.01 total for both side]\n");
    printf("                     -tr   FLOAT  trim rate for the fragment length distribution[0.02 total for both side]\n");
    printf("                     -mq   INT    minimum mapping quality for a normal read pair\n");
    printf("                     -help        print this help message\n");
    exit(0);
}

// clean up the resouses used in the split-read build program
void TGM_ReadPairScanClean(TGM_ReadPairScanPars* pScanPars)
{
    fclose(pScanPars->fileListInput);
    free(pScanPars->workingDir);
    free(pScanPars->specialPrefix);
}
