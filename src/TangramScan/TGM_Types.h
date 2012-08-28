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
    TGM_A = 1,
    TGM_C = 2,
    TGM_G = 4,
    TGM_T = 8,
    TGM_N = 15

}TGM_Base;

// search direction for query region relative to the positon of anchor mate
typedef enum
{
    TGM_UPSTREAM,     // search the query region uptream to the anchor mate (smaller coordinate)

    TGM_DOWNSTREAM    // search the query region downstream to the anchor mate (larger coordinate)

}TGM_Direction;

// action applied on a certain DNA sequence
typedef enum
{
    TGM_FORWARD,

    TGM_REVERSE_COMP,    // change the sequence to its reverse complement format

    TGM_INVERSE,         // change the sequence to its inverse format

    TGM_COMP

}TGM_SeqAction;

typedef enum
{
    TGM_1F = 0,
    TGM_1R = 1,
    TGM_2F = 2,
    TGM_2R = 3

}TGM_SingleMode;

typedef enum
{
    TGM_BAD_PAIR_MODE = -1,
    TGM_1F2F = 0,
    TGM_1F2R = 1,
    TGM_1R2F = 2,
    TGM_1R2R = 3,
    TGM_2F1F = 4,
    TGM_2F1R = 5,
    TGM_2R1F = 6,
    TGM_2R1R = 7

}TGM_PairMode;

typedef enum
{
    STREAM_PASS   = -1,

    STREAM_RETRUN = 0,

    STREAM_KEEP = 1

}TGM_StreamCode;

typedef enum
{
    ORDINARY_ONLY = 0,

    ALL_SV = 1

}TGM_DetectMode;

typedef TGM_SeqAction TGM_Strand;

#define MD5_CHECKSUM_LEN 16

#define MD5_STR_LEN 32

#define MAX_HASH_SIZE 12

#define TGM_MAX_LINE 1024

#define TGM_EMPTY 0

#define NUM_TOTAL_PAIR_MODE 8

#define NUM_ALLOWED_PAIR_MODE 2

#define NUM_ALLOWED_HIST 2

// mismatch code for cigar
#define BAM_CMISMATCH 8

#define TGM_BLOCK_SIZE 5120

#endif  /*TGM_TYPES_H*/
