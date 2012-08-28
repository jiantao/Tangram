/*
 * =====================================================================================
 *
 *       Filename:  TGM_Types.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/08/2011 21:43:41
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef  TGM_TYPES_H
#define  TGM_TYPES_H

#define NUM_SV_TYPES 6

typedef enum
{
    FALSE = 0,
    TRUE  = 1

}TGM_Bool;

typedef enum
{
    TGM_OK  =           0,
    TGM_EOF =          -1,
    TGM_ERR =          -2,
    TGM_OVER_FLOW =    -97,
    TGM_FULL      =    -98,
    TGM_NOT_FOUND =    -99,
    TGM_OUT_OF_RANGE = -100

}TGM_Status;

typedef enum
{
    SV_NORMAL             = 0,

    SV_DELETION           = 1,

    SV_TANDEM_DUP         = 2,

    SV_INVERSION          = 3,

    SV_SPECIAL            = 4,

    SV_INTER_CHR_TRNSLCTN = 5

}SV_EventType;


#endif  /*TGM_TYPES_H*/
