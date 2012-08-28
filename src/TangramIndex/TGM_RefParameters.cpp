/*
 * =====================================================================================
 *
 *       Filename:  TGM_RefParameters.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  08/27/2012 08:28:51 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include <stdlib.h>

#include "TGM_GetOpt.h"
#include "TGM_Error.h"
#include "TGM_RefParameters.h"

using namespace Tangram;

// total number of arguments we should expect for the split-read build program
#define OPT_TOTAL_ARGS       4

// total number of required arguments we should expect for the split-read build program
#define OPT_REQUIRED_ARGS    3

// the index of show help in the option object array
#define OPT_HELP                 0

// the index of the fasta input file in the option object array

#define OPT_REF_INPUT            1

#define OPT_SP_REF_INPUT         2

#define OPT_OUTPUT               3


RefPars::RefPars()
{
    fpRefInput = NULL;
    fpSpRefInput = NULL;
    fpOutput = NULL;
}

RefPars::~RefPars()
{
    gzclose(fpRefInput);
    gzclose(fpSpRefInput);
    fclose(fpOutput);
}

void RefPars::Set(const char** argv, int argc)
{
    TGM_Option opts[] = 
    {
        {"help", NULL, FALSE},
        {"ref",  NULL, FALSE},
        {"sp",  NULL, FALSE},
        {"out",  NULL, FALSE},
        {NULL,   NULL, FALSE}
    };

    int optNum = TGM_GetOpt(opts, argc, argv);

    if (optNum < OPT_REQUIRED_ARGS)
        ShowHelp();

    for (unsigned int i = 0; i != OPT_TOTAL_ARGS; ++i)
    {
        switch (i)
        {
            case OPT_HELP:
                if (opts[i].isFound)
                    ShowHelp();
                break;
            case OPT_REF_INPUT:
                if (opts[i].value == NULL)
                    TGM_ErrQuit("ERROR: The reference file is not specified.\n");

                fpRefInput = gzopen(opts[i].value, "r");
                if (fpRefInput == NULL)
                    TGM_ErrQuit("ERROR: Cannot open reference file\n");

                break;
            case OPT_SP_REF_INPUT:
                if (opts[i].value == NULL)
                    TGM_ErrQuit("ERROR: The special reference file is not specified.\n");

                fpSpRefInput = gzopen(opts[i].value, "r");
                if (fpSpRefInput == NULL)
                    TGM_ErrQuit("ERROR: Cannot open special reference file\n");

                break;
            case OPT_OUTPUT:
                if (opts[i].value == NULL)
                    TGM_ErrQuit("ERROR: The output file is not specified.\n");

                fpOutput = fopen(opts[i].value, "wb");
                if (fpOutput == NULL)
                    TGM_ErrQuit("ERROR: Cannot open output file\n");

                break;
            default:
                TGM_ErrQuit("ERROR: Unrecognized argument.\n");
                break;
        }
    }
}

void RefPars::ShowHelp(void) const
{
    printf("Usage: tangram_index -ref <input_ref> -sp <input_special_ref> -out <ouput_file>\n\n");

    printf("Mandatory arguments: -ref  FILE  input of reference file\n");
    printf("                     -sp   FILE  input of special reference file\n");
    printf("                     -out  FILE  output file of indexed reference\n");
    printf("                     -help       print this help message\n");

    exit(EXIT_SUCCESS);
}
