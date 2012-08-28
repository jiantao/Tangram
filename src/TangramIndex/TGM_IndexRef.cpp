/*
 * =====================================================================================
 *
 *       Filename:  TGM_IndexRef.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  08/27/2012 08:13:07 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */


#include "TGM_RefParameters.h"
#include "TGM_Reference.h"

using namespace Tangram;

int main(int argc, char *argv[])
{
    RefPars refPars;
    refPars.Set((const char**) argv, argc);

    Reference reference;
    reference.Create(refPars.fpRefInput, refPars.fpSpRefInput, refPars.fpOutput);

    return EXIT_SUCCESS;
}
