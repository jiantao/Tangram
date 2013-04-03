/*
 * =====================================================================================
 *
 *       Filename:  TGM_ReadPairScanGetOpt.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  04/17/2012 12:44:31 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef  TGM_READPAIRSCANGETOPT_H
#define  TGM_READPAIRSCANGETOPT_H

#include <stdio.h>
#include <stdint.h>
#include "TGM_Types.h"

// parameters used for read pair build
typedef struct
{
    double cutoff;                 // fragment length cutoff p-value

    double trimRate;               // trim rate from fragment length distribution

    unsigned char minMQ;           // minimum mapping quality for a read pair

    FILE* fileListInput;           // input stream of a file list containing all the bam file names

    char* workingDir;              // working directory for the detector

    char* specialPrefix;           // prefix of the special reference

    uint32_t prefixLen;            // length of the prefix of the special reference

    uint32_t minFrags;             // minimum number of normal fragments in a library

}TGM_ReadPairScanPars;

// set the parameters for the split-read build program from the parsed command line arguments 
void TGM_ReadPairScanSetPars(TGM_ReadPairScanPars* pScanPars, int argc, char* argv[]);

void TGM_ReadPairScanHelp(void);

// clean up the resouses used in the split-read build program
void TGM_ReadPairScanClean(TGM_ReadPairScanPars* pScanPars);

#endif  /*TGM_READPAIRSCANGETOPT_H*/
