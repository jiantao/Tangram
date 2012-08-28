/*
 * =====================================================================================
 *
 *       Filename:  TGM_GetOpt.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  04/17/2012 12:24:18 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */
#ifndef  TGM_GETOPT_H
#define  TGM_GETOPT_H

#include <stdio.h>
#include "TGM_Types.h"

// an object hold the options parsed from command line arguments
typedef struct TGM_Option
{
    char* name;          // name of this option, following a '-' character in the command line

    const char* value;   // a pointer points to the value of the option in the argv array

    TGM_Bool isFound;    // boolean variable to indicate that if we fount this option or not

}TGM_Option;

// get the options from command line arguemnts
int TGM_GetOpt(TGM_Option opts[], int argc, char* argv[]);

#endif  /*TGM_GETOPT_H*/

