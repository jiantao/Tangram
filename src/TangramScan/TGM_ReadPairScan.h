/*
 * =====================================================================================
 *
 *       Filename:  TGM_ReadPairScan.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  04/11/2012 02:44:37 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef  TGM_READPAIRSCAN_H
#define  TGM_READPAIRSCAN_H

#include "TGM_LibInfo.h"
#include "TGM_ReadPairScanGetOpt.h"

void TGM_ReadPairScan(const TGM_ReadPairScanPars* pScanPars);

void TGM_SpecialIDUpdate(TGM_SpecialID* pSpecialID, const TGM_ZAtag* pZAtag);

#endif  /*TGM_READPAIRSCAN_H*/
