/*
 * =====================================================================================
 *
 *       Filename:  TGM_GetOpt.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  04/17/2012 12:25:24 AM
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

// get the options from command line arguemnts
int TGM_GetOpt(TGM_Option opts[], int argc, const char* argv[])
{
    int optNum = 0;
    TGM_Bool hasParsed = FALSE;

    for (int i = 1; i < argc; ++i)
    {
        if (hasParsed)
        {
            hasParsed = FALSE;
            continue;
        }

        if (argv[i][0] != '-')
            TGM_ErrQuit("ERROR: Invalid argument %s.\n", argv[i]);

        const char* currOpt = argv[i] + 1;

        for (unsigned int j = 0; ; ++j)
        {
            if (opts[j].name == NULL)
                TGM_ErrQuit("ERROR: Ivalid argument \"%s\".\n", argv[i]);

            if (strcmp(currOpt, opts[j].name) == 0)
            {
                opts[j].isFound = TRUE;
                if (i + 1 != argc && argv[i+1][0] != '-')
                {
                    hasParsed = TRUE;
                    opts[j].value = argv[i + 1];
                }
                else
                    opts[j].value = NULL;

                ++optNum;
                break;
            }
        }
    }

    return optNum;
}
