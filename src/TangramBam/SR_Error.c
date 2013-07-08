#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>		// for definition of errno 
#include <stdarg.h>		// ANSI C header file 
#include "SR_Error.h"

static void SR_ErrDoit(int, const char *, va_list);
char *pname = NULL;		// caller can set this from argv[0] 

// Nonfatal error related to a system call.
// Print a message and return. 
void SR_ErrRet(const char *fmt, ...)
{
    va_list		ap;

    va_start(ap, fmt);
    SR_ErrDoit(1, fmt, ap);
    va_end(ap);
    return;
}

// Fatal error related to a system call.
// Print a message and terminate. 
void SR_ErrSys(const char *fmt, ...)
{
    va_list		ap;

    va_start(ap, fmt);
    SR_ErrDoit(1, fmt, ap);
    va_end(ap);
    exit(1);
}

// Fatal error related to a system call.
// Print a message, dump core, and terminate. 
void SR_ErrDump(const char *fmt, ...)
{
    va_list		ap;

    va_start(ap, fmt);
    SR_ErrDoit(1, fmt, ap);
    va_end(ap);
    abort();		// dump core and terminate 
    exit(1);		// shouldn't get here 
}

// Nonfatal error unrelated to a system call.
// Print a message and return. 
void SR_ErrMsg(const char *fmt, ...)
{
    va_list		ap;

    va_start(ap, fmt);
    SR_ErrDoit(0, fmt, ap);
    va_end(ap);
    return;
}

// Fatal error unrelated to a system call.
// Print a message and terminate. 
void SR_ErrQuit(const char *fmt, ...)
{
    va_list		ap;

    va_start(ap, fmt);
    SR_ErrDoit(0, fmt, ap);
    va_end(ap);
    exit(1);
}

// Print a message and return to caller.
// Caller specifies "errnoflag". 
static void SR_ErrDoit(int errnoflag, const char *fmt, va_list ap)
{
    int		errno_save;
    char	buf[MAXLINE_ERR];

    errno_save = errno;		// value caller might want printed 
    vsprintf(buf, fmt, ap);
    if (errnoflag)
        sprintf(buf+strlen(buf), ": %s", strerror(errno_save));
    strcat(buf, "\n");
    fflush(stdout);		// in case stdout and stderr are the same 
    fputs(buf, stderr);
    fflush(NULL);		// flushes all stdio output streams 
    return;
}
