#ifndef  SR_ERROR_H
#define  SR_ERROR_H

// max length for error message
#define	MAXLINE_ERR 4096

// Fatal error related to a system call.
// Print a message, dump core, and terminate. 
void SR_ErrDump(const char *, ...);	

// Nonfatal error unrelated to a system call.
// Print a message and return. 
void SR_ErrMsg(const char *, ...);

// Fatal error unrelated to a system call.
// Print a message and terminate. 
void SR_ErrQuit(const char *, ...);

// Nonfatal error related to a system call.
// Print a message and return. 
void SR_ErrRet(const char *, ...);

// Fatal error related to a system call.
// Print a message and terminate. 
void SR_ErrSys(const char *, ...);

#endif  //SR_ERROR_H
