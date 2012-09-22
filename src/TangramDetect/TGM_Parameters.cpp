/*
 * =====================================================================================
 *
 *       Filename:  TGM_Parameters.cpp
 *
 *    Description:  Set the parameters with arguments from the command line
 *
 *        Created:  07/28/2012 01:55:55 PM
 *       Revision:  none
 *       Compiler:  g++
 *
 *         Author:  Jiantao Wu (), 
 *   Inistitution:  Boston College
 *
 * =====================================================================================
 */

#include <cstdlib>

#include "TGM_Version.h"
#include "TGM_Parameters.h"
#include "TGM_Error.h"
#include "TGM_GetOpt.h"
#include "TGM_Utilities.h"

using namespace std;
using namespace Tangram;
using namespace BamTools;

// total number of arguments we should expect for the split-read build program
#define OPT_TOTAL_ARGS       15

// total number of required arguments we should expect for the split-read build program
#define OPT_REQUIRED_ARGS    3

// the index of show help in the option object array
#define OPT_HELP                 0

// the index of the fasta input file in the option object array

#define OPT_LIB_INPUT            1

#define OPT_HIST_INPUT           2

#define OPT_INPUT_FILE           3

#define OPT_RANGE_STR            4

#define OPT_REF_INPUT            5

#define OPT_CHECK_LIB            6

#define OPT_MIN_CLUSTER_SIZE     7

#define OPT_MIN_EVENT_LEN        8

#define OPT_MIN_MQ               9

#define OPT_SPECIAL_MIN_MQ       10

#define OPT_DETECT_SET           11

#define OPT_MIN_SCORE_RATE       12

#define OPT_THREAD_NUM           13

#define OPT_OUTPUT               14


#define DEFAULT_MIN_CLUSTER_SIZE 2

#define DEFAULT_MIN_EVENT_LEN 100

#define DEFAULT_MIN_SOFT_SIZE 10

#define DEFAULT_MIN_MQ 20

#define DEFAULT_SPECIAL_MIN_MQ 20

#define DEFAULT_DETECT_SET 0xffffffff

#define DETECT_SET_MASK 0x0000001f

#define MAX_REF_NAME_LEN 512

#define MAX_LINE 1024

DetectPars::DetectPars()
{

    fpLibInput = NULL;

    fpHistInput = NULL;

    fpBamListInput = NULL;

    pRangeStr = NULL;

    detectSet = DEFAULT_DETECT_SET;

    refID = -1;

    range[0] = -1;

    range[1] = -1;

    numThread = DEFAULT_THREAD_NUM;

    minSoftSize = DEFAULT_MIN_SOFT_SIZE;

    minClusterSize = DEFAULT_MIN_CLUSTER_SIZE;

    minEventLength = DEFAULT_MIN_EVENT_LEN;

    maxFragDiff = 0;

    useSplitRead = false;

    minMQ = DEFAULT_MIN_MQ;

    spMinMQ = DEFAULT_SPECIAL_MIN_MQ;
}

DetectPars::~DetectPars()
{
    Clean();
}

void DetectPars::Clean(void)
{
    fclose(fpLibInput);
    fclose(fpHistInput);
    fclose(fpBamListInput);
}



Parameters::Parameters(DetectPars& dPars, AlignerPars& aPars)
                : detectPars(dPars), alignerPars(aPars)
{

}

Parameters::~Parameters()
{

}

void Parameters::Set(const char** argv, int argc)
{
    TGM_Option opts[] = 
    {
        {"help", NULL, FALSE},
        {"lb",   NULL, FALSE},
        {"ht",   NULL, FALSE},
        {"in",   NULL, FALSE},
        {"rg",   NULL, FALSE},
        {"ref",   NULL, FALSE},
        {"cl",   NULL, FALSE},
        {"mcs",   NULL, FALSE},
        {"mel",   NULL, FALSE},
        {"mq",   NULL, FALSE},
        {"smq",   NULL, FALSE},
        {"dt",  NULL, FALSE},
        {"msr",  NULL, FALSE},
        {"p",  NULL, FALSE},
        {"out",  NULL, FALSE},
        {NULL,   NULL, FALSE}
    };

    int optNum = TGM_GetOpt(opts, argc, argv);
    if (optNum < OPT_REQUIRED_ARGS)
        ShowHelp();

    for (unsigned int i = 0; i != OPT_TOTAL_ARGS; ++i)
    {
        switch (i)
        {
            case OPT_HELP:
                if (opts[i].isFound)
                    ShowHelp();
                break;
            case OPT_LIB_INPUT:
                if (opts[i].value == NULL)
                    TGM_ErrQuit("ERROR: The library information file is not specified.\n");

                detectPars.fpLibInput = fopen(opts[i].value, "rb");
                if (detectPars.fpLibInput == NULL)
                    TGM_ErrQuit("ERROR: Cannot open the libraray information file \"%s\" for reading.\n", opts[i].value);

                break;
            case OPT_HIST_INPUT:
                if (opts[i].value == NULL)
                    TGM_ErrQuit("ERROR: The histogram file is not specified.\n");

                detectPars.fpHistInput = fopen(opts[i].value, "rb");
                if (detectPars.fpHistInput == NULL)
                    TGM_ErrQuit("ERROR: Cannot open the histogram file \"%s\" for read.\n", opts[i].value);

                break;
            case OPT_INPUT_FILE:
                if (opts[i].value != NULL)
                {
                    detectPars.fpBamListInput = fopen(opts[i].value, "r");
                    if (detectPars.fpBamListInput == NULL)
                        TGM_ErrQuit("ERROR: Cannot open file \"%s\" for read.\n", opts[i].value);
                }
                else
                    TGM_ErrQuit("ERROR: Bam input is not specified.\n");

                break;
            case OPT_RANGE_STR:
                if (opts[i].value != NULL)
                    detectPars.pRangeStr = opts[i].value;

                break;
            case OPT_REF_INPUT:
                if (opts[i].value != NULL)
                {
                    alignerPars.fpRefInput = fopen(opts[i].value, "rb");
                    if (alignerPars.fpRefInput == NULL)
                        TGM_ErrQuit("ERROR: Cannot open file \"%s\" for read.\n", opts[i].value);

                    detectPars.useSplitRead = true;
                }

                break;

            case OPT_CHECK_LIB:
                if (opts[i].value != NULL)
                    detectPars.maxFragDiff = atoi(opts[i].value);
                else
                    detectPars.maxFragDiff = 0;

                break;
            case OPT_MIN_CLUSTER_SIZE:
                if (opts[i].value != NULL)
                {
                    detectPars.minClusterSize = atoi(opts[i].value);
                    if (detectPars.minClusterSize <= 0)
                        TGM_ErrQuit("ERROR: %s is not a valid argument for minimum cluster size.\n");
                }

                break;

            case OPT_MIN_EVENT_LEN:
                if (opts[i].value != NULL)
                {
                    detectPars.minEventLength = atoi(opts[i].value);
                    if (detectPars.minEventLength <= 0)
                        TGM_ErrQuit("ERROR: %s is not a valid argument for minimum event length.\n", opts[i].value);
                }

                break;

            case OPT_MIN_MQ:
                if (opts[i].value != NULL)
                {
                    int tempMQ = atoi(opts[i].value);
                    if (tempMQ < 0 || tempMQ > 255)
                        TGM_ErrQuit("ERROR: %s is an invalid minimum mapping quality.\n", opts[i].value);

                    detectPars.minMQ = tempMQ;
                }

                break;
            case OPT_SPECIAL_MIN_MQ:
                if (opts[i].value != NULL)
                {
                    int tempMQ = atoi(opts[i].value);
                    if (tempMQ < 0 || tempMQ > 255)
                        TGM_ErrQuit("ERROR: %s is an invalid minimum mapping quality.\n", opts[i].value);

                    detectPars.spMinMQ = tempMQ;
                }

                break;
            case OPT_DETECT_SET:
                if (opts[i].value != NULL)
                {
                    if (opts[i].value[0] == '0')
                        sscanf(opts[i].value, "%x", &(detectPars.detectSet));
                    else
                        sscanf(opts[i].value, "%d", &(detectPars.detectSet));
                }

                // the detection set can not be 0, otherwise no events will be called
                if ((detectPars.detectSet & DETECT_SET_MASK) == 0)
                    TGM_ErrQuit("ERROR: Please specify at least one type of SV to detect.\n");

                alignerPars.detectSet = detectPars.detectSet;

                break;
            case OPT_MIN_SCORE_RATE:
                if (opts[i].value != NULL)
                    alignerPars.minScoreRate = atof(opts[i].value);

                if (alignerPars.minScoreRate <= 0.0 || alignerPars.minScoreRate > 1.0)
                    TGM_ErrQuit("ERROR: %s is an invalid minimum score rate. (0.0 - 1.0]\n", opts[i].value);

                break;
            case OPT_THREAD_NUM:
                if (opts[i].value != NULL)
                {
                    int numThread = atoi(opts[i].value);
                    if (numThread <= 0)
                        TGM_ErrQuit("ERROR: Invalid number of threads.\n");

                    detectPars.numThread = numThread;
                    alignerPars.numThread = numThread;
                }

                break;
            case OPT_OUTPUT:
                detectPars.outputPrefix = opts[i].value;
                break;
            default:
                TGM_ErrQuit("ERROR: Unrecognized argument.\n");
                break;
        }
    }
}

void Parameters::ParseRangeStr(const BamMultiReader& multiReader)
{
    const char* pRangeStr = detectPars.pRangeStr;
    if (pRangeStr == NULL)
        return;

    char buff[MAX_REF_NAME_LEN];
    const char* rangePos = NULL;
    for (unsigned int i = 0; i != MAX_REF_NAME_LEN; ++i)
    {
        if (pRangeStr[i] != '\0' && pRangeStr[i] != ':')
            buff[i] = pRangeStr[i];
        else
        {
            buff[i] = '\0';
            rangePos = pRangeStr + i;
            break;
        }
    }

    detectPars.refID = 0;
    detectPars.refID = multiReader.GetReferenceID(buff);
    if (detectPars.refID < 0)
        TGM_ErrQuit("ERROR: %s is not a valid reference name.\n", buff);

    if (*rangePos == ':')
    {
        ++rangePos;
        detectPars.range[0] = atoi(rangePos);
        if (detectPars.range[0] < 0)
            TGM_ErrQuit("ERROR: Invalid range %s\n", detectPars.pRangeStr);

        const char* nextRangePos = strpbrk(rangePos, "-");
        if (nextRangePos != NULL)
        {
            ++nextRangePos;
            detectPars.range[1] = atoi(nextRangePos);
            if (detectPars.range[0] < 0)
                TGM_ErrQuit("ERROR: Invalid range %s\n", pRangeStr);

            if (detectPars.range[0] > detectPars.range[1])
                TGM_ErrQuit("ERROR: Invalid range %s\n", pRangeStr);
        }
    }

    detectPars.range[0] -= 1;
}

void Parameters::SetRange(BamMultiReader& multiReader, const int32_t& maxFragLen) const
{
    if (detectPars.refID < 0)
        return;

    int32_t start = detectPars.range[0];
    int32_t end = detectPars.range[1];

    if (start < 0)
        start = 0;

    const RefVector& refVector = multiReader.GetReferenceData();
    if (end < 0)
    {
        end = refVector[detectPars.refID].RefLength;
    }
    else
    {
        if (end + 2 * maxFragLen - 1 < refVector[detectPars.refID].RefLength)
            end += 2 * maxFragLen - 1;
        else
            end = refVector[detectPars.refID].RefLength;
    }

    if (!multiReader.SetRegion(detectPars.refID, start, detectPars.refID, end))
        TGM_ErrQuit("ERROR: Cannot set the detection region.\n");
}

void Parameters::SetBamFilenames(vector<string>& filenames)
{
    string bamFilename;
    char line[MAX_LINE];
    filenames.reserve(100);

    while (TGM_GetNextLine(line, MAX_LINE, detectPars.fpBamListInput) == TGM_OK)
    {
        bamFilename.assign(line);
        filenames.push_back(bamFilename);
    }
}

void Parameters::ShowHelp(void) const
{
    printf("\nProgram: tangram (Toolbox for structural variation detection)\n");
    printf("Version: %s\n\n", TGM_VERSION);

    printf("Usage: tangram_detect [options] -lb <library_info> -ht <fragment_length_histogram> -in <bam_list> -rg <regions>\n\n");

    printf("Mandatory arguments: -lb   FILE   library information file\n");
    printf("                     -ht   FILE   fragment length histogram file\n");
    printf("                     -in   FILE   list of all input bam files\n");
    printf("                     -rg   STRING chromosome region (all the bam files must be sorted by chromosome positions and indexed)\n\n");

    printf("Options:             -ref  FILE   transfered reference sequence, required for split alignment[null]\n");
    printf("                     -out  STR    prefix to the output files, including the path.[stdout]\n");
    printf("                     -cl   INT    check for invalid libraries\n");
    printf("                     -mcs  INT    minimum cluster size[2]\n");
    printf("                     -mel  INT    minimum event lenth[100]\n");
    printf("                     -mq   INT    minimum mapping quality for pairs other than special pairs[20]\n");
    printf("                     -smq  INT    minimum mapping quality for special pairs[20]\n");
    printf("                     -dt   INT    detection set [0xffffffff: report all types of SV]\n");
    printf("                     -msr  FLOAT  minimum score rate for split alignments[0.8]\n");
    printf("                     -p    INT    number of processors (threads)[1]\n");
    printf("                     -help        print this help message\n");

    printf("Notes:\n\n");

    printf("  1. A region should be presented in one of the following formats:\n\
    `1', `2:1000' and `X:1000-2000' (1-based). When a region is specified,\n\
    the input alignment file must be an indexed BAM file.\n\n");

    printf("  2. Detection set is a bit set to indicate which types of SV will be detected.\n");
    printf("     Each bit in this bit set corresponding to a type of SV event:\n\n");
    printf("         DELETION        0x1 (hexadecimal) or 1 (decimal)\n");
    printf("         TANDUM DUP      0x2 (hexadecimal) or 2 (decimal)\n");
    printf("         INVERSION       0x4 (hexadecimal) or 4 (decimal)\n");
    printf("         MEI             0x8 (hexadecimal) or 8 (decimal)\n\n");
    printf("     To detect multiple types of SV events, just add all the corresponding numbers up.\n");
    printf("     For example, if you only want to call DELETIONS and MEI, then the input value for\n");
    printf("     `-dt' option will be 0x1(1) + 0x8(8) = 0x9(9). The `-dt' option can take either decimal\n");
    printf("     or hexadecimal (start with `0x') number for input.\n\n");

    exit(EXIT_SUCCESS);
}
