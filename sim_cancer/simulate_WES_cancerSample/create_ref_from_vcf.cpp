// ==========================================================================
//                             create_ref_from_vcf
// ==========================================================================
// Copyright (c) 20014-2015, Jochen Singer, ETH Zürich
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Jochen Singer <jochen.singer@bsse.ethz.ch>
// ==========================================================================

#include <map>
#include <iostream>
#include <fstream>
#include <tuple>

#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/stream.h>
#include <seqan/basic.h>
#include <seqan/arg_parse.h>
#include <seqan/vcf_io.h>

using namespace seqan;

struct AppOptions
{
    String<char> refFileIn;
    String<char> refFileOut;
    String<char> vcfFileIn;
    String<char> mapFileOut;
};

ArgumentParser::ParseResult
parseCommandLine(AppOptions & options, int argc, char const ** argv)
{
    // Setup ArgumentParser.
    ArgumentParser parser("get_ref_from_vcf");
    // Set short description, version, and date.
    setShortDescription(parser, "Generating a fasta reference file from a vcf file.");
    setVersion(parser, "0.1");
    setDate(parser, "2015");

    ArgParseArgument refInArgument(ArgParseArgument::INPUT_FILE, "IN");
    setValidValues(refInArgument, "fa fasta fna");
    addArgument(parser, refInArgument);

    ArgParseArgument vcfInArgument(ArgParseArgument::INPUT_FILE, "IN");
    setValidValues(vcfInArgument, "vcf");
    addArgument(parser, vcfInArgument);

    ArgParseArgument refOutArgument(ArgParseArgument::OUTPUT_FILE, "OUT");
    setValidValues(refOutArgument, "fa fasta fna");
    addArgument(parser, refOutArgument);

    ArgParseArgument vcfOutArgument(ArgParseArgument::OUTPUT_FILE, "OUT");
    setValidValues(vcfOutArgument, "csv");
    addArgument(parser, vcfOutArgument);

    // Parse command line.
    ArgumentParser::ParseResult res = parse(parser, argc, argv);

    // Only extract  options if the program will continue after parseCommandLine()
    if (res != ArgumentParser::PARSE_OK)
        return res;

    // Extract option values.
    getArgumentValue(options.refFileIn, parser, 0);
    getArgumentValue(options.vcfFileIn, parser, 1);
    getArgumentValue(options.refFileOut, parser, 2);
    getArgumentValue(options.mapFileOut, parser, 3);

    return ArgumentParser::PARSE_OK;
}

int main(int argc, const char* argv[]){

    typedef String<Iupac>    TRefSeq;

    //---------------------------------------------------------------
    // parse the command line.
    //---------------------------------------------------------------
    ArgumentParser parser;
    AppOptions options;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If there was an error parsing or built-in argument parser functionality
    // was triggered then we exit the program.  The return code is 1 if there
    // were errors and 0 if there were none.
    if (res != ArgumentParser::PARSE_OK)
        return res == ArgumentParser::PARSE_ERROR;

    std::cout << "create_ref_from_vcf\n"
              << "===================\n\n";

    StringSet<CharString> refIds;
    StringSet<TRefSeq> refSeqs;
    SeqFileIn seqFileIn(toCString(options.refFileIn));
    // the map helps to find the correct reference in refIds and refSeqs
    std::map<CharString, int> refMap;
    std::cout << "Loading references ... \n";
    for (unsigned i = 0; !atEnd(seqFileIn); ++i)
    {
        resize(refIds, length(refIds) + 1);
        resize(refSeqs, length(refSeqs) + 1);
        readRecord(back(refIds), back(refSeqs), seqFileIn);
        refMap.insert(std::pair<CharString, int>(back(refIds), i));
    }
    std::cout << "done!\n";

    shrinkToFit(refIds);
    shrinkToFit(refSeqs);

    StringSet<TRefSeq> refSeqsOut;
    resize(refSeqsOut, length(refSeqs), Exact());

    StringSet<String<std::tuple<unsigned, unsigned> > > mappingPos;
    resize(mappingPos, length(refSeqs));
    std::ofstream mappingFile(toCString(options.mapFileOut));
    if (!mappingFile.is_open())
        std::cout << "Unable to open file\n";

    seqan::VcfFileIn vcfStream(toCString(options.vcfFileIn));
    seqan::VcfHeader header;
    readHeader(header, vcfStream);

    unsigned tempRefId = -1;
    unsigned currentPos = -1;
    String<char> currentChrom = "";
    VcfRecord record;
    VcfRecord longestRecord;

    SeqFileOut seqFileOut(toCString(options.refFileOut));

    while (!atEnd(vcfStream))
    {
        readRecord(record, vcfStream);
        if (currentChrom != contigNames(context(vcfStream))[record.rID])
        {
            if (currentChrom != "")
            {
                append(refSeqsOut[tempRefId], suffix(refSeqs[tempRefId], currentPos));
            }
            currentChrom = contigNames(context(vcfStream))[record.rID];
            tempRefId=refMap[contigNames(context(vcfStream))[record.rID]];
            if (infix(refSeqs[tempRefId], record.beginPos, record.beginPos + length(record.ref)) != record.ref)
                std::cerr << "Sequences do not overlap! " << infix(refSeqs[tempRefId], record.beginPos, record.beginPos + length(record.ref)) << " " << record.ref << std::endl;
            append(refSeqsOut[tempRefId], prefix(refSeqs[tempRefId], record.beginPos));
            appendValue(mappingPos[tempRefId], std::tuple<unsigned, unsigned>(record.beginPos, length(refSeqsOut[tempRefId])));
            append(refSeqsOut[tempRefId], record.alt);
            longestRecord = record;
            currentPos = record.beginPos + length(record.ref);
        }
        else
        {
            // this if is necessary because vcf entries may overlap
            // no overlap
            if (record.beginPos >= longestRecord.beginPos + length(longestRecord.ref))
            {
                if (infix(refSeqs[tempRefId], record.beginPos, record.beginPos + length(record.ref)) != record.ref)
                    std::cerr << "Sequences do not overlap!" << std::endl;
                append(refSeqsOut[tempRefId], infix(refSeqs[tempRefId], currentPos, record.beginPos));
                appendValue(mappingPos[tempRefId], std::tuple<unsigned, unsigned>(record.beginPos, length(refSeqsOut[tempRefId])));
                append(refSeqsOut[tempRefId], record.alt);
                longestRecord = record;
                currentPos = record.beginPos + length(record.ref);
            }
            else
            {
                std::cerr << "This should not have happened! - Two variants overlap!" << std::endl;
                unsigned overlapStart = record.beginPos - longestRecord.beginPos;
                if (length(longestRecord.alt) < overlapStart + length(record.alt))
                {
                    int overlapEnd = length(longestRecord.alt) - overlapStart;
                    if (overlapEnd < 0)
                    {
                        std::cerr << "overlapEnd is negative" << std::endl;
                    }

                    if (suffix(longestRecord.alt, overlapStart) != prefix(record.alt, overlapEnd))
                    {
                        std::cerr << "Problem! -> Overlap" << std::endl;
                    }
                    append(refSeqsOut[tempRefId], suffix(record.alt, overlapEnd));
                    appendValue(mappingPos[tempRefId], std::tuple<unsigned, unsigned>(record.beginPos, length(refSeqsOut[tempRefId]) - length(record.alt)));
                    longestRecord = record;
                    currentPos = record.beginPos + length(record.ref);
                }
                // new record is contained
                else
                {
                    unsigned overlapEnd = length(record.alt);
                    appendValue(mappingPos[tempRefId], std::tuple<unsigned, unsigned>(record.beginPos, length(refSeqsOut[tempRefId]) - length(longestRecord.alt) + overlapStart));
                    if (infix(longestRecord.alt, overlapStart, overlapStart + overlapEnd) != record.alt)
                    {
                        std::cout << "Problem! -> Contained" << std::endl;
                    }
                }
            }
        }
    }
    append(refSeqsOut[tempRefId], suffix(refSeqs[tempRefId], currentPos));

    for (unsigned i = 0; i < length(refSeqsOut); ++i)
    {

        if (!empty(refSeqsOut[i]))
        {
            writeRecord(seqFileOut, refIds[i], refSeqsOut[i]);
            for (unsigned j = 0; j < length(mappingPos[i]); ++j)
            {
                mappingFile << refIds[i] << "\t";
                mappingFile << std::get<0>(mappingPos[i][j]) << "\t";
                mappingFile <<std::get<1>(mappingPos[i][j]) << "\n";
            }

        }
        else
        {
            writeRecord(seqFileOut, refIds[i], refSeqs[i]);
        }
    }

    return 0;
}



