/*
 * =====================================================================================
 *
 *       Filename:  TGM_RefParameters.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  08/27/2012 08:24:29 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef  TGM_REFPARAMETERS_H
#define  TGM_REFPARAMETERS_H

#include <stdio.h>
#include "zlib.h"

namespace Tangram
{
    class RefPars
    {
        public:
            RefPars();
            ~RefPars();

            void Set(const char** argv, int argc);

        private:

            void ShowHelp(void) const;

        public:

            gzFile fpRefInput;
            gzFile fpSpRefInput;
            FILE*  fpOutput;
    };
};

#endif  /*TGM_REFPARAMETERS_H*/
