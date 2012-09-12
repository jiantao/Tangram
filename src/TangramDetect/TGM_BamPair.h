/*
 * =====================================================================================
 *
 *       Filename:  TGM_BamPair.h
 *
 *    Description:  Table used to store bam pairs
 *
 *        Created:  05/09/2012 06:52:05 PM
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Jiantao Wu (), 
 *   Inistitution:  Boston College
 *
 * =====================================================================================
 */

#ifndef  TGM_BAMPAIR_H
#define  TGM_BAMPAIR_H

#include <string>

#include "api/BamAlignment.h"

//#include "TGM_ObjPool.h"
#include "TGM_Parameters.h"
#include "TGM_LibTable.h"
#include "TGM_Sequence.h"
#include "TGM_FragLenTable.h"

#define BAM_CIGAR_SHIFT 4
#define BAM_CIGAR_MASK  ((1 << BAM_CIGAR_SHIFT) - 1)

/*
  CIGAR operations.
 */
/*! @abstract CIGAR: match */
#define BAM_CMATCH      0
/*! @abstract CIGAR: insertion to the reference */
#define BAM_CINS        1
/*! @abstract CIGAR: deletion from the reference */
#define BAM_CDEL        2
/*! @abstract CIGAR: skip on the reference (e.g. spliced alignment) */
#define BAM_CREF_SKIP   3
/*! @abstract CIGAR: clip on the read with clipped sequence present in qseq */
#define BAM_CSOFT_CLIP  4
/*! @abstract CIGAR: clip on the read with clipped sequence trimmed off */
#define BAM_CHARD_CLIP  5
/*! @abstract CIGAR: padding */
#define BAM_CPAD        6
/*! @abstract CIGAR: sequence match */
#define BAM_CSEQ_MATCH  7        
/*! @abstract CIGAR: sequence mismatch */
#define BAM_CMISMATCH   8

namespace Tangram
{
    // a map used to map the pair mode into its corresponding number
    // negative value means invalid mode
    static const int8_t TGM_PairModeMap[16] = 
    { 
        -1, -1, 0, 1,

        -1, -1, 2, 3,

        4, 5, -1, -1,

        6, 7, -1, -1
    };  

    static const int8_t TGM_SeqTechMap[4][8] = 
    {
        {0, 1, 2, 3, 4, 5, 6, 7},     // ILLUMINA

        {2, 3, 0, 1, 5, 4, 7, 6},     // 454
 
        {1, 0, 3, 2, 6, 7, 4, 5},     // SOLiD

        {3, 2, 1, 0, 7, 6, 5, 4}      // ILLUMINA LONG
    };

    typedef enum
    {
        TGM_1F = 0,   // first mate forward
        TGM_1R = 1,   // first mate reversed
        TGM_2F = 2,   // second mate forward
        TGM_2R = 3    // second mate reversed

    }SingleOrient;

    typedef enum
    {
        TGM_UP = 0,   // up stream

        TGM_DOWN = 1  // down stream

    }Direction;

    static const int8_t TGM_SingleSeqTechMap[4][4] = 
    {
        {TGM_2R, TGM_2F, TGM_1R, TGM_1F},

        {TGM_2F, TGM_2R, TGM_1F, TGM_1R},

        {TGM_2F, TGM_2R, TGM_1F, TGM_1R},

        {TGM_2R, TGM_2F, TGM_1R, TGM_1F},
    };

    static const int8_t TGM_SingleDirectionhMap[4][4] = 
    {
        {TGM_DOWN, TGM_UP, TGM_DOWN, TGM_UP},

        {TGM_UP, TGM_DOWN, TGM_DOWN, TGM_UP},

        {TGM_DOWN, TGM_UP, TGM_UP, TGM_DOWN},

        {TGM_UP, TGM_DOWN, TGM_UP, TGM_DOWN},
    };

    typedef enum
    {
        TGM_BAD_PAIR_MODE = -1, // bad pair mode
        TGM_1F2F = 0,           // first mate forward, second mate forward (the position of first and second mate here corresponding to their genomic position)
        TGM_1F2R = 1,           // first mate forward, second mate reversed
        TGM_1R2F = 2,           // first mate reversed, second mate forward
        TGM_1R2R = 3,           //  .
        TGM_2F1F = 4,           //  .
        TGM_2F1R = 5,           //  .
        TGM_2R1F = 6,           //  .
        TGM_2R1R = 7            //  .

    }PairOrient;

    typedef enum
    {
        PT_UNKNOWN    = -1,
        PT_NORMAL     = 0,
        PT_LONG       = 1,
        PT_SHORT      = 2,
        PT_REVERSED   = 3,
        PT_INVERTED3  = 4,
        PT_INVERTED5  = 5,
        PT_SPECIAL3   = 6,
        PT_SPECIAL5   = 7,
        PT_CROSS      = 8,
        PT_SOFT3      = 9,
        PT_SOFT5      = 10

    }PairType;

    static const int8_t PairTypeMap[2][8] = 
    {
        {PT_INVERTED3, PT_NORMAL, PT_UNKNOWN, PT_INVERTED5, PT_UNKNOWN, PT_UNKNOWN, PT_REVERSED, PT_UNKNOWN},

        {PT_UNKNOWN, PT_UNKNOWN, PT_REVERSED, PT_UNKNOWN, PT_INVERTED3, PT_NORMAL, PT_UNKNOWN, PT_INVERTED5}
    };

    // status report of a bam pair
    struct BamPairStat
    {
        uint8_t bestMQ[2];         // mapping quality

        uint8_t secMQ[2];          // depreciated features

        int16_t numMM[2];          // number of mismatches

        uint16_t numMappings[2];   // number of mappings in the whole genome

        int32_t end[2];            // ending position of each mate

        uint32_t readGrpID;        // name of the read group

        int32_t fragLen;           // fragment length of the pair

        char spRef[2][3];          // special reference name

        int8_t readPairType;       // read pair type (for sv type)

        int softSize[2][2];        // soft clipping size

        PairOrient orient;         // orientation mode of the pair
    };

    struct LocalPair
    {
        int32_t readGrpID;                  // read group ID

        int32_t refID:16, fragLenQual:16;   // reference ID, fragment length quality (two-sided p value)

        int32_t pos[2];                     // mapping position for each mate (0-based)

        int32_t end[2];                     // mapping end for each mate (0-based)

        int16_t numMM[2];                   // number of mismatches of the up mate, number of mismatches of the down
                                                                                                                                             
        uint8_t bestMQ[2];                  // best mapping quality of the up and down mate

        int8_t readPairType;                // read pair type (for sv type)

        int8_t orient;                      // fragment length quality(-10log(p-value)), read pairt type, pair mode
    };

    // special pair structure (one mate is unique the other mate hits the special reference)
    struct SpecialPair
    {
        int32_t readGrpID;                                           // read group ID
                                                                                                                                             
        int16_t refID[2];                                            // up reference ID, down reference ID
                                                                                                                                             
        int32_t pos[2];                                              // alignment position of the up and down mate
                                                                     
        int32_t end[2];                                              // alignment end of the up and down mate

        int16_t numMM[2];                                            // number of mismatches of the up mate, number of mismatches of the down
                                                                                                                                             
        uint8_t bestMQ[2];                                           // best mapping quality of the up and down mate

        uint16_t numSpeicalHits;                                     // number of hits on the special reference

        uint32_t specialID:16, pairType:8, orient:6, aStrand:1, sStrand:1;   // special reference ID, pair type, pair orientation, anchor mate strand, mutiple mate strand
    };

    // orphan pair structure (one mate is unique the other mate is unmapped)
    struct OrphanPair
    {
        int32_t readGrpID;                        // read group ID

        int32_t anchorPos;                        // mapping position of the anchor mate

        int32_t anchorEnd;                        // mapping end of the anchor mate

        int32_t refID:20, aOrient:4, bestMQ:8;    // reference ID, anchor mate orientation, mapping quality

        TGM_Sequence read;                        // read sequence of the orphan read
    };

    struct SoftPair
    {
        int32_t readGrpID;                                                    // read group ID

        uint32_t refID:16, readPairType:8, orient:6, aStrand:1, sStrand:1;    // reference ID, read pair type, pair orientation, anchor mate strand, soft mate strand

        int32_t pos;                                                          // mapping position of anchor mate

        int32_t end;                                                          // mapping end of the anchor mate

        int32_t matePos;                                                      // mapping position of the soft mate

        uint32_t* cigar;                                                      // cigar information of the soft mate

        int32_t cigarLen;                                                     // cigar length

        TGM_Sequence read;                                                    // read sequence of the soft mate
    };

    typedef struct
    {
        uint32_t poolIdx;

        uint32_t idx;

    }ObjIndex;

    class BamPairTable
    {
        public:
            BamPairTable(const DetectPars& detectPars, const LibTable& libTable, const FragLenTable& inFragLenTable);
            ~BamPairTable();

            // update the bam pair table with the incoming alignment
            void Update(const BamTools::BamAlignment& alignment);

        private:

            // basic filter of bam alignment
            bool BamPairFilter(void) const;

            // check if the bam pair is orphan pair
            inline bool IsOrphanPair(void) const
            {
                if (!pAlignment->IsMapped() || !pAlignment->IsMateMapped())
                    return true;

                return false;
            }

            // set the pair status with the information from a bam alignment
            bool SetPairStat(void);

            // parse the za tag in the bam alignment
            void ParseZAstr(bool isUpMate);

            // calculate the number of mismatches from the za tag
            int GetNumMismatchFromZA(int32_t* pLen, int softSize[2], const char* cigarStr, unsigned int cigarLen, const char* mdStr, unsigned int mdLen);

            // calculate the number of mismatches from the bam alignment
            int GetNumMismatchFromBam(int softSize[2]);
            
            // get the pair type
            PairType GetPairType(bool isUpMate);

            // get the pair orientation
            PairOrient GetPairOrient(bool isUpMate);

            // update the orphan pair array
            void UpdateOrphanPair(void);

            // update the local pair array
            void UpdateLocalPair(Array<LocalPair>& localPairs);

            // update the inverted pair array
            void UpdateInvertedPair(void);

            // update the special pair array
            void UpdateSpecialPair(void);

            // update the cross pair array
            void UpdateCrossPair(void);

            // update the soft pair array
            void UpdateSoftPair(void);

            // transfer the cigar format
            uint32_t* TransferCigar(const std::vector<BamTools::CigarOp>& cigarData);

        public:
            // long pair array (deletions)
            Array<LocalPair> longPairs;

            // short pair array (insertions or duplications)
            Array<LocalPair> shortPairs;

            // reversed pair array (duplications)
            Array<LocalPair> reversedPairs;

            // inverted pair array (inversions)
            Array<LocalPair> invertedPairs;

            // special pair array (MEI insertions)
            Array<SpecialPair> specialPairs;

            // orphan pair array (candidates for split read algorithm)
            Array<OrphanPair> orphanPairs;

            // soft pair array (candidates for split read algorithm)
            Array<SoftPair> softPairs;

            int numInverted3;

        private:

            // void* readNameHash;

            // ObjPool<std::string> readNames;

            // ObjIndex objIndex;
            
            // detector parameters
            const DetectPars& detectPars;

            // reference to a library information table
            const LibTable& libTable;

            // reference to a fragment length distribution table
            const FragLenTable& fragLenTable;

            // pointer to the incoming bam alignment
            const BamTools::BamAlignment* pAlignment;

            // bam pair status
            BamPairStat pairStat;

            // read group name
            std::string readGrpName;

            // ZA string
            std::string zaStr;
    };
};

#endif  /*TGM_BAMPAIR_H*/
