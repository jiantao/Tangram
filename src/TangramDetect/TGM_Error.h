#ifndef  TGM_ERROR_H
#define  TGM_ERROR_H

// max length for error message
#define	MAXLINE_ERR 4096

#ifdef __cplusplus
extern "C" {
#endif	// __cplusplus

// Fatal error related to a system call.
// Print a message, dump core, and terminate. 
void TGM_ErrDump(const char *, ...);	

// Nonfatal error unrelated to a system call.
// Print a message and return. 
void TGM_ErrMsg(const char *, ...);

// Fatal error unrelated to a system call.
// Print a message and terminate. 
void TGM_ErrQuit(const char *, ...);

// Nonfatal error related to a system call.
// Print a message and return. 
void TGM_ErrRet(const char *, ...);

// Fatal error related to a system call.
// Print a message and terminate. 
void TGM_ErrSys(const char *, ...);

#ifdef __cplusplus
}
#endif	// __cplusplus

#endif  //TGM_ERROR_H
