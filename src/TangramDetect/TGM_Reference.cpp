/*
 * =====================================================================================
 *
 *       Filename:  TGM_Reference.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  05/17/2012 04:05:26 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jiantao Wu (), 
 *        Company:  
 *
 * =====================================================================================
 */

#include <cstring>
#include <map>

#include "md5.h"
#include "kseq.h"
#include "khash.h"
#include "TGM_Error.h"
#include "TGM_Reference.h"

KSEQ_INIT(gzFile, gzread)

using namespace Tangram;
using namespace std;

// maximum number of character will be load in a line from the fasta file
#define MAX_REF_LINE 1024

#define MD5_STR_LEN 32

#define MD5_CHECKSUM_LEN 16

namespace {
static int CompareInt(const void* a, const void* b)
{
    const int* first = (const int*) a;
    const int* second = (const int*) b;

    return (*first - *second);
}
} // end namespace

namespace Tangram {
static void SeqTransfer(char* seq, unsigned int len)
{
    for (unsigned int i = 0; i != len; ++i)
        seq[i] = nt_table[ (int) seq[i]];
}
} // end namespace Tangram

Reference::Reference()
{

}

Reference::~Reference()
{
    for (unsigned int i = 0; i != refHeader.names.Size(); ++i)
        free(refHeader.names[i]);

    for (unsigned int i = 0; i != spRefHeader.names.Size(); ++i)
        free(spRefHeader.names[i]);
}

void Reference::Create(gzFile fpRefFastaInput, gzFile fpSpRefFastaInput, FILE* fpRefOutput)
{
    int64_t headerPos = 0;

    unsigned int writeSize = fwrite(&headerPos, sizeof(int64_t), 1, fpRefOutput);
    if (writeSize != 1)
        TGM_ErrQuit("ERROR: Cannot write the reference header position.\n");


    bool hasPadding = true;
    if (fpSpRefFastaInput != NULL)
        CreateRef(spRefHeader, fpSpRefFastaInput, fpRefOutput, hasPadding);

    hasPadding = false;
    CreateRef(refHeader, fpRefFastaInput, fpRefOutput, hasPadding);

    WriteHeader(fpRefOutput);
}

void Reference::CreateRef(RefHeader& header, gzFile fpRefFastaInput, FILE* fpRefOutput, bool hasPadding)
{
    InitRefHeader(header);

    int8_t padding[DEFAULT_PADDING_LEN];
    memset(padding, 4, DEFAULT_PADDING_LEN);

    uint64_t seqLen = 0;
    int ret = 0;
    kseq_t* seq = kseq_init(fpRefFastaInput);

    bool isFirst = true;

    while ((ret = kseq_read(seq)) >= 0)
    {
        if (seq->name.l <= 0 || seq->seq.l <= 0)
            continue;

        unsigned int writeSize = 0;
        if (hasPadding && !isFirst)
        {
            writeSize = fwrite(padding, sizeof(int8_t), DEFAULT_PADDING_LEN, fpRefOutput);
            if (writeSize != DEFAULT_PADDING_LEN)
                TGM_ErrQuit("ERROR: Cannot write the reference sequence into the file.\n");

            seqLen += DEFAULT_PADDING_LEN;
        }

        isFirst = false;

        unsigned int numRef = header.names.Size();
        if (header.names.IsFull())
        {
            header.names.Resize(numRef * 2);
            header.endPos.Resize(numRef * 2);
            header.md5.Resize(numRef * MD5_STR_LEN * 2 + 1);
        }

        SetName(header, seq->name.s, seq->name.l);

        seqLen += seq->seq.l;
        header.endPos[numRef] = seqLen - 1;
        header.endPos.Increment();

        SetMd5(header, seq->seq.s, seq->seq.l);

        SeqTransfer(seq->seq.s, seq->seq.l);

        writeSize = fwrite(seq->seq.s, sizeof(int8_t), seq->seq.l, fpRefOutput);
        if (writeSize != seq->seq.l)
            TGM_ErrQuit("ERROR: Cannot write the reference sequence into the file.\n");
    }

    kseq_destroy(seq);
}

void Reference::InitRefHeader(RefHeader& header)
{
    header.names.Init(20);
    header.endPos.Init(20);
    header.md5.Init(20 * MD5_STR_LEN + 1);
}

void Reference::SetName(RefHeader& header, const char* buff, int len)
{
    header.names.End() = (char*) malloc((len + 1) * sizeof(char));
    if (header.names.End() == NULL)
        TGM_ErrQuit("ERROR: Not enough memory for the reference name.\n");

    strcpy(header.names.End(), buff);
    header.names.Increment();
}

void Reference::SetMd5(RefHeader& header, char* sequence, uint32_t seqLen)
{
    unsigned char MD5[MD5_CHECKSUM_LEN];
    memset(MD5, 0, MD5_CHECKSUM_LEN);

    MD5_CTX context;
    MD5Init(&context);
    MD5Update(&context, (unsigned char*) sequence, seqLen);
    MD5Final(MD5, &context);

    unsigned int md5Len = header.md5.Size();
    char* md5String = header.md5.GetPointer(md5Len);
    for (unsigned int i = 0; i != MD5_CHECKSUM_LEN; ++i)
    {
        sprintf(md5String, "%02X", MD5[i]);
        md5String += 2;
    }

    header.md5.SetSize(md5Len + MD5_STR_LEN);
}

void Reference::WriteHeader(FILE* fpRefOutput)
{
    int64_t filePos = ftello(fpRefOutput);

    uint32_t numRef = spRefHeader.names.Size();
    unsigned int writeSize = fwrite(&numRef, sizeof(uint32_t), 1, fpRefOutput);
    if (writeSize != 1)
        TGM_ErrQuit("ERROR: Cannot write the number of special references into the file.\n");

    for (unsigned int i = 0; i != numRef; ++i)
    {
        uint32_t nameLen = strlen(spRefHeader.names[i]);
        writeSize = fwrite(&nameLen, sizeof(uint32_t), 1, fpRefOutput);
        if (writeSize != 1)
            TGM_ErrQuit("ERROR: Cannot write the length special reference name into the file.\n");

        writeSize = fwrite(spRefHeader.names[i], sizeof(char), nameLen, fpRefOutput);
        if (writeSize != nameLen)
            TGM_ErrQuit("ERROR: Cannot write the special reference names into the file.\n");
    }

    writeSize = fwrite(spRefHeader.endPos.GetPointer(0), sizeof(int64_t), numRef, fpRefOutput);
    if (writeSize != numRef)
        TGM_ErrQuit("ERROR: Cannot write the end position of special reference into the file.\n");

    writeSize = fwrite(spRefHeader.md5.GetPointer(0), sizeof(char), numRef * MD5_STR_LEN, fpRefOutput);
    if (writeSize != numRef * MD5_STR_LEN)
        TGM_ErrQuit("ERROR: Cannot write the MD5 of special reference into the file.\n");

    numRef = refHeader.names.Size();
    writeSize = fwrite(&numRef, sizeof(uint32_t), 1, fpRefOutput);
    if (writeSize != 1)
        TGM_ErrQuit("ERROR: Cannot write the number of special references into the file.\n");

    for (unsigned int i = 0; i != numRef; ++i)
    {
        uint32_t nameLen = strlen(refHeader.names[i]);
        writeSize = fwrite(&nameLen, sizeof(uint32_t), 1, fpRefOutput);
        if (writeSize != 1)
            TGM_ErrQuit("ERROR: Cannot write the length special reference name into the file.\n");

        writeSize = fwrite(refHeader.names[i], sizeof(char), nameLen, fpRefOutput);
        if (writeSize != nameLen)
            TGM_ErrQuit("ERROR: Cannot write the special reference names into the file.\n");
    }

    writeSize = fwrite(refHeader.endPos.GetPointer(0), sizeof(int64_t), numRef, fpRefOutput);
    if (writeSize != numRef)
        TGM_ErrQuit("ERROR: Cannot write the end position of special reference into the file.\n");

    writeSize = fwrite(refHeader.md5.GetPointer(0), sizeof(char), numRef * MD5_STR_LEN, fpRefOutput);
    if (writeSize != numRef * MD5_STR_LEN)
        TGM_ErrQuit("ERROR: Cannot write the MD5 of special reference into the file.\n");

    fseeko(fpRefOutput, 0, SEEK_SET);

    writeSize = fwrite(&filePos, sizeof(int64_t), 1, fpRefOutput);
    if (writeSize != 1)
        TGM_ErrQuit("ERROR: Cannot write position of reference header into the file.\n");
}

void Reference::Read(FILE* fpRefInput, const int32_t& refID, const int32_t& start, const int32_t& end)
{
    int64_t headerPos = 0;
    unsigned int readSize = fread(&headerPos, sizeof(int64_t), 1, fpRefInput);
    if (readSize != 1)
        TGM_ErrQuit("ERROR: Cannot read the header position.\n");

    fseeko(fpRefInput, headerPos, SEEK_SET);

    ReadRefHeader(fpRefInput);

    fseeko(fpRefInput, sizeof(int64_t), SEEK_SET);

    if (spRefHeader.names.Size() > 0)
        ReadSpecialRef(fpRefInput);

    ReadRef(fpRefInput, refID, start, end);
}

void Reference::ReadRefHeader(FILE* fpRefInput)
{
    uint32_t numRef = 0;;
    unsigned int readSize = fread(&numRef, sizeof(uint32_t), 1, fpRefInput);
    if (readSize != 1)
        TGM_ErrQuit("ERROR: Cannot read the number of special references from the file.\n");

    if (numRef > 0)
    {
        spRefHeader.names.Init(numRef);
        spRefHeader.endPos.Init(numRef);
        spRefHeader.md5.Init(numRef * MD5_STR_LEN);

        spRefHeader.names.SetSize(numRef);
        spRefHeader.endPos.SetSize(numRef);
        spRefHeader.md5.SetSize(numRef * MD5_STR_LEN);

        for (unsigned int i = 0; i != numRef; ++i)
        {
            uint32_t nameLen = 0;
            readSize = fread(&nameLen, sizeof(uint32_t), 1, fpRefInput);
            if (readSize != 1)
                TGM_ErrQuit("ERROR: Cannot read the length special reference name from the file.\n");

            spRefHeader.names[i] = (char*) calloc(sizeof(char), nameLen + 1);
            if (spRefHeader.names[i] == NULL)
                TGM_ErrQuit("ERROR: Not enough memory for the refernece name.\n");

            readSize = fread(spRefHeader.names[i], sizeof(char), nameLen, fpRefInput);
            if (readSize != nameLen)
                TGM_ErrQuit("ERROR: Cannot read the special reference names from the file.\n");
        }

        readSize = fread(spRefHeader.endPos.GetPointer(0), sizeof(int64_t), numRef, fpRefInput);
        if (readSize != numRef)
            TGM_ErrQuit("ERROR: Cannot read the end position of special reference from the file.\n");

        readSize = fread(spRefHeader.md5.GetPointer(0), sizeof(char), numRef * MD5_STR_LEN, fpRefInput);
        if (readSize != numRef * MD5_STR_LEN)
            TGM_ErrQuit("ERROR: Cannot read the MD5 of special reference from the file.\n");
    }

    numRef = refHeader.names.Size();
    readSize = fread(&numRef, sizeof(uint32_t), 1, fpRefInput);
    if (readSize != 1)
        TGM_ErrQuit("ERROR: Cannot read the number of special references from the file.\n");

    refHeader.names.Init(numRef);
    refHeader.endPos.Init(numRef);
    refHeader.md5.Init(numRef * MD5_STR_LEN);

    refHeader.names.SetSize(numRef);
    refHeader.endPos.SetSize(numRef);
    refHeader.md5.SetSize(numRef * MD5_STR_LEN);

    for (unsigned int i = 0; i != numRef; ++i)
    {
        uint32_t nameLen = 0;
        readSize = fread(&nameLen, sizeof(uint32_t), 1, fpRefInput);
        if (readSize != 1)
            TGM_ErrQuit("ERROR: Cannot read the length special reference name from the file.\n");

        refHeader.names[i] = (char*) calloc(sizeof(char), nameLen + 1);
        if (refHeader.names[i] == NULL)
            TGM_ErrQuit("ERROR: Not enough memory for the refernece name.\n");

        readSize = fread(refHeader.names[i], sizeof(char), nameLen, fpRefInput);
        if (readSize != nameLen)
            TGM_ErrQuit("ERROR: Cannot read the special reference names from the file.\n");
    }

    readSize = fread(refHeader.endPos.GetPointer(0), sizeof(int64_t), numRef, fpRefInput);
    if (readSize != numRef)
        TGM_ErrQuit("ERROR: Cannot read the end position of special reference from the file.\n");

    readSize = fread(refHeader.md5.GetPointer(0), sizeof(char), numRef * MD5_STR_LEN, fpRefInput);
    if (readSize != numRef * MD5_STR_LEN)
        TGM_ErrQuit("ERROR: Cannot read the MD5 of special reference from the file.\n");
}

void Reference::ReadSpecialRef(FILE* fpRefInput)
{
    int64_t spRefLen = spRefHeader.endPos.Last() + 1;

    spRefSeq.Init(spRefLen);
    spRefSeq.SetSize(spRefLen);

    unsigned int readSize = fread(spRefSeq.GetPointer(0), sizeof(int8_t), spRefLen, fpRefInput);
    if (readSize != spRefLen)
        TGM_ErrQuit("ERROR: Cannot read the special reference from the file.\n");

    CreatFamily();
}

void Reference::ReadRef(FILE* fpRefInput, const int32_t& refID, int32_t start, int32_t end)
{
    int64_t refBegin = GetRefBeginPos(refID);
    int64_t refEnd = refHeader.endPos[refID];
    int64_t refLen = refEnd - refBegin + 1;

    if (start < 0)
        start = 0;

    if (end < 0)
        end = refLen - 1;

    if (start >= refLen || end >= refLen)
        TGM_ErrQuit("ERROR: Invalid reference region.\n");

    int64_t regionBegin = start + refBegin;
    int64_t regionEnd = end + refBegin;

    uint64_t regionLen = regionEnd - regionBegin + 1;

    refSeq.Init(regionLen);
    refSeq.SetSize(regionLen);

    int ret = fseeko(fpRefInput, regionBegin, SEEK_CUR);
    if (ret < 0)
        TGM_ErrQuit("ERROR: Cannot jump int the reference file.\n");

    unsigned int readSize = fread(refSeq.GetPointer(0), sizeof(int8_t), regionLen, fpRefInput);
    if (readSize != regionLen)
        TGM_ErrQuit("ERROR: Cannot read the reference from the file.\n");

    this->refID = refID;
    pos = start;
}

void Reference::CreatFamily(void)
{
    unsigned int spNameSize = spRefHeader.names.Size();
    map<string, int> nameHash;

    familyMap.Init(spNameSize);
    familyMap.SetSize(spNameSize);
    string name;
    int nameIndex = 0;

    for (unsigned int i = 0; i != spNameSize; ++i)
    {
        const char* pName = strpbrk(spRefHeader.names[i], "_");

        if (pName != NULL)
        {
            ++pName;
            const char* end = pName;
            int len = 0;

            while (*end != '\0')
            {
                if (*end == '.')
                    break;

                ++end;
            }

            len = end - pName;

            if (len > 0)
            {
                name.assign(pName, len);
                map<string, int>::iterator itor = nameHash.find(name);
                if (itor != nameHash.end())
                    familyMap[i] = itor->second;
                else
                {
                    familyName.push_back(name);
                    nameHash[name] = nameIndex;
                    familyMap[i] = nameIndex;

                    ++nameIndex;
                }
            }
            else
                familyMap[i] = -1;
        }
        else
            familyMap[i] = -1;
    }
}

const char* Reference::GetSpRefName(int& spRefID, int32_t spRefPos) const
{
    spRefID = spRefHeader.endPos.UpperBound(spRefPos, CompareInt);
    if (spRefID >= 0)
        return spRefHeader.names[spRefID];
    else
        return NULL;
}

void Reference::CreateFamilyToZA(const Array<char*>& spZANames)
{
    int spZANum = spZANames.Size();
    int numFamily = familyName.size();

    if (numFamily == 0 || spZANum == 0)
        return;

    familyToZA.Init(numFamily);
    familyToZA.SetSize(numFamily);

    ZAToFamily.Init(spZANum);
    ZAToFamily.SetSize(spZANum);

    for (int i = 0; i != spZANum; ++i)
        ZAToFamily[i] = -1;

    for (int i = 0; i != numFamily; ++i)
    {
        familyToZA[i] = -1;
        for (int j = 0; j != spZANum; ++j)
        {
            const char* name = familyName[i].c_str();
            const char* zaName = spZANames[j];

            if (toupper(zaName[0]) == toupper(name[0]) && toupper(zaName[1]) == toupper(name[1]))
            {
                familyToZA[i] = j;
                ZAToFamily[j] = i;
            }
        }
    }
}
