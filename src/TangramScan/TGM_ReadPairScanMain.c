/*
 * =====================================================================================
 *
 *       Filename:  TGM_ReadPairScanMain.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  04/16/2012 11:39:18 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include "TGM_ReadPairScan.h"
#include "TGM_ReadPairScanGetOpt.h"

int main(int argc, char *argv[])
{
    TGM_ReadPairScanPars scanPars;
    TGM_ReadPairScanSetPars(&scanPars, argc, argv);

    TGM_ReadPairScan(&scanPars);

    TGM_ReadPairScanClean(&scanPars);

    return EXIT_SUCCESS;
}
