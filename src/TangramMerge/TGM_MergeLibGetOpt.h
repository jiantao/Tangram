/*
 * =====================================================================================
 *
 *       Filename:  TGM_MergeLibGetOpt.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  04/18/2012 12:42:52 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef  TGM_MERGELIBGETOPT_H
#define  TGM_MERGELIBGETOPT_H

// parameters used for read pair build
typedef struct
{
    char* workingDir;              // working directory for the detector

}TGM_MergeLibPars;

// set the parameters for the split-read build program from the parsed command line arguments 
void TGM_MergeLibSetPars(TGM_MergeLibPars* pMergePars, int argc, char* argv[]);

void TGM_MergeLibHelp(void);

void TGM_MergeLibClean(TGM_MergeLibPars* pMergePars);

#endif  /*TGM_MERGELIBGETOPT_H*/
