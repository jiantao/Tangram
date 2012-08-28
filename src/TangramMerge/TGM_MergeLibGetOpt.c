/*
 * =====================================================================================
 *
 *       Filename:  TGM_MergeLibGetOpt.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  04/18/2012 12:44:36 AM
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
#include "TGM_MergeLibGetOpt.h"

// total number of arguments we should expect for the split-read build program
#define OPT_MERGE_TOTAL_NUM 2

// total number of required arguments we should expect for the split-read build program
#define OPT_MERGE_REQUIRED_NUM 1

// the index of show help in the option object array
#define OPT_HELP           0

// the index of the fasta input file in the option object array
#define OPT_WORKING_DIR    1

void TGM_MergeLibSetPars(TGM_MergeLibPars* pMergePars, int argc, char* argv[])
{
    TGM_Option opts[] = 
    {
        {"help", NULL, FALSE},
        {"dir",   NULL, FALSE},
        {NULL,   NULL, FALSE}
    };

    int optNum = TGM_GetOpt(opts, argc, argv);
    if (optNum < OPT_MERGE_REQUIRED_NUM)
        TGM_MergeLibHelp();

    for (unsigned int i = 0; i != OPT_MERGE_TOTAL_NUM; ++i)
    {
        switch (i)
        {
            case OPT_HELP:
                if (opts[i].isFound)
                    TGM_MergeLibHelp();
                break;
            case OPT_WORKING_DIR:
                if (opts[i].value == NULL)
                    TGM_ErrQuit("ERROR: The working directory is not specified.\n");

                unsigned int dirLen = strlen(opts[i].value);
                if (opts[i].value[dirLen - 1] == '/')
                    pMergePars->workingDir = strdup(opts[i].value);
                else
                {
                    pMergePars->workingDir = (char*) malloc(dirLen + 2);
                    if (pMergePars->workingDir == NULL)
                        TGM_ErrQuit("ERROR: Not enough memory for the working dir name.\n");

                    sprintf(pMergePars->workingDir, "%s/", opts[i].value);
                }

                break;
            default:
                TGM_ErrQuit("ERROR: Unrecognized argument.\n");
                break;
        }
    }
}

void TGM_MergeLibHelp(void)
{
    printf("Usage: tangram_merge -dir <input_dir>\n\n");

    printf("Mandatory arguments: -dir  STRING the path to the dir contains all the fragment length distribution files\n");
    printf("                     -help        print this help message\n");

    exit(0);
}

void TGM_MergeLibClean(TGM_MergeLibPars* pMergePars)
{
    free(pMergePars->workingDir);
}
