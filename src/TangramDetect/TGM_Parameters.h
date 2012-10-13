/*
 * =====================================================================================
 *
 *       Filename:  TGM_Parameters.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  07/28/2012 12:45:57 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#ifndef  TGM_PARAMETERS_H
#define  TGM_PARAMETERS_H

#include <cstdio>
#include <stdint.h>
#include <vector>
#include <string>

#include "api/BamMultiReader.h"

namespace Tangram
{
    static const int8_t default_score_matrix[] = { 2, -2, -2, -2, -2, // A
                                                  -2,  2, -2, -2, -2, // C
                                                  -2, -2,  2, -2, -2, // G
                                                  -2, -2, -2,  2, -2, // T
                                                  -2, -2, -2, -2, -2};// N

    static const int8_t second_score_matrix[] = {  2, -4, -4, -4, -4, // A
                                                  -4,  2, -4, -4, -4, // C
                                                  -4, -4,  2, -4, -4, // G
                                                  -4, -4, -4,  2, -4, // T
                                                  -4, -4, -4, -4, -4};// N

    static const unsigned char DEFAULT_MIN_MQ = 20;

    static const uint8_t DEFAULT_GAP_OPEN = 3;

    static const uint8_t DEFAULT_GAP_EXT = 1;

    static const uint8_t DEFAULT_SEC_GAP_OPEN = 6;

    static const uint8_t DEFAULT_SEC_GAP_EXT = 4;

    // static const int DEFAULT_MIN_PARTIAL_SIZE = 10;

    static const float DEFAULT_MAX_MISMATCH_RATE = 0.1;

    static const int DEFAULT_MIN_ALIGNED_LEN = 15;

    static const int DEFAULT_MIN_TRIGGER_LEN = 30;

    static const int EDGE_TOLERANCE = 5;

    static const double DEFAULT_MIN_SCORE_RATE = 0.8;

    static const double DEFAULT_MIN_ENTROPY = 1.2;

    static const double DEFAULT_MIN_COVER_RATE = 0.85;

    static const double MIN_POLYA_RATE = 0.85;

    static double DEFAULT_MIN_AGREE_RATE = 0.5;

    static const unsigned int DEFAULT_PARTIAL_SIZE = 20;

    static const unsigned int DEFAULT_SV_TYPE_NUM = 5;

    static const int MIN_POLYA_SIZE = 5;
    
    static const int DEFAULT_THREAD_NUM = 1;

    static const int32_t DEFAULT_MIN_JUMP_LEN = -1;

    static const int DEFAULT_MIN_RP_FRAG = 2;

    static const int DEFAULT_MIN_SR_FRAG = 5;

    class DetectPars
    {
        public:

            DetectPars();
            ~DetectPars();

            void Clean(void);

        public:

            FILE* fpLibInput;

            FILE* fpHistInput;

            FILE* fpBamListInput;

            const char* outputPrefix;

            const char* pRangeStr;

            uint32_t detectSet;

            int refID;

            int32_t range[2];

            int numThread;

            int minSoftSize;

            int minClusterSize;

            int minEventLength;

            int maxFragDiff;

            bool useSplitRead;

            unsigned char minMQ;

            unsigned char spMinMQ;
    };

    struct AlignerPars
    {
        const int8_t* mat;

        const int8_t* secMat;

        FILE* fpRefInput;

        uint8_t gapOpen;

        uint8_t gapExt;

        uint8_t secGapOpen;

        uint8_t secGapExt;

        int8_t flag;

        int16_t scoreFilter;

        int16_t distFilter;

        uint32_t detectSet;

        int minAlignedLen;

        int minTriggerLen;

        int numThread;

        double minScoreRate;

        double minEntropy;

        double maxMisMatchRate;

        double minCoverRate;

        double minAgreeRate;

        AlignerPars()
        {
            mat = default_score_matrix;
            secMat = second_score_matrix;
            fpRefInput = NULL;
            gapOpen = DEFAULT_GAP_OPEN;
            gapExt = DEFAULT_GAP_EXT;
            secGapOpen = DEFAULT_SEC_GAP_OPEN;
            secGapExt = DEFAULT_SEC_GAP_EXT;
            flag = 0;
            flag |= 0x08;
            flag |= 0x0f;
            scoreFilter = 0;
            distFilter = 32767;
            detectSet = 0xffffffff;
            minAlignedLen = DEFAULT_MIN_ALIGNED_LEN;
            minTriggerLen = DEFAULT_MIN_TRIGGER_LEN;
            numThread = DEFAULT_THREAD_NUM;
            minScoreRate = DEFAULT_MIN_SCORE_RATE;
            minEntropy = DEFAULT_MIN_ENTROPY;
            maxMisMatchRate = DEFAULT_MAX_MISMATCH_RATE;
            minCoverRate = DEFAULT_MIN_COVER_RATE;
            minAgreeRate = DEFAULT_MIN_AGREE_RATE;
        }

        ~AlignerPars()
        {
            if (fpRefInput != NULL)
                fclose(fpRefInput);
        }
    };

    struct GenotypePars
    {
        bool doGenotype;

        unsigned char minMQ;

        int minCrossLen;

        unsigned int minRpFrag;

        unsigned int minSrFrag;

        int32_t minJumpLen;

        double p[3];  // parameters for calculating the binomial pdf

        GenotypePars()
        {
            doGenotype = false;

            minMQ = DEFAULT_MIN_MQ;

            minRpFrag = DEFAULT_MIN_RP_FRAG;

            minSrFrag = DEFAULT_MIN_SR_FRAG;

            minCrossLen = DEFAULT_MIN_ALIGNED_LEN;

            minJumpLen = DEFAULT_MIN_JUMP_LEN;

            p[0] = 0.99;
            p[1] = 0.5;
            p[2] = 0.01;
        }
    };

    class Parameters
    {
        public:

            Parameters(DetectPars& detectPars, AlignerPars& aligerPars, GenotypePars& genotypePars);

            ~Parameters();

            void Set(const char** argv, int argc);

            void ShowHelp(void) const;

            void ParseRangeStr(const BamTools::BamMultiReader& multiReader);

            void SetRange(BamTools::BamMultiReader& multiReader, const int32_t& maxFragLen) const;

            void SetBamFilenames(std::vector<std::string>& filenames);

            void Clean(void);

        public:

            DetectPars& detectPars;

            AlignerPars& alignerPars;

            GenotypePars& genotypePars;
    };
};

#endif  /*TGM_PARAMETERS_H*/
