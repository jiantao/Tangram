/*
 * =====================================================================================
 *
 *       Filename:  TGM_MergeLibMain.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  04/18/2012 12:20:03 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include "TGM_LibInfo.h"
#include "TGM_MergeLibGetOpt.h"

int main(int argc, char *argv[])
{
    TGM_MergeLibPars mergeLibPars;

    TGM_MergeLibSetPars(&mergeLibPars, argc, argv);

    TGM_LibInfoTableMerge(mergeLibPars.workingDir);

    TGM_MergeLibClean(&mergeLibPars);

    return EXIT_SUCCESS;
}
