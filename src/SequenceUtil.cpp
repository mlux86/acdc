#include <algorithm>

#include "SequenceUtil.h"
#include "Logger.h"

#include <seqan/alignment_free.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

void SequenceUtil::allPermsRepetition(std::vector<std::string> & perms, const std::vector<char> & alphabet, std::string elem, const unsigned length)
{
    if (length == 0)
    {
        perms.push_back(elem);
        return;
    }

    for (auto c : alphabet)
    {
        SequenceUtil::allPermsRepetition(perms, alphabet, elem + c, length - 1);
    }
}

std::vector<std::string> SequenceUtil::allKmers(const unsigned kmerLength)
{
    std::vector<char> alphabet = {'A', 'C', 'G', 'T'};
    std::vector<std::string> perms;
    SequenceUtil::allPermsRepetition(perms, alphabet, "", kmerLength);
    return perms;
}

std::string SequenceUtil::reverseComplement(const std::string & seq)
{
    std::string revCompl(seq);
    std::reverse(revCompl.begin(), revCompl.end());

    for (auto & c : revCompl)
    {
        switch(c)
        {
            case 'A':
                c = 'T';
                break;
            case 'C':
                c = 'G';
                break;
            case 'G':
                c = 'C';
                break;
            case 'T':
                c = 'A';
                break;
        }
    }

    return revCompl;
}

std::vector<std::string> SequenceUtil::filterFasta(const std::string & fasta, const std::set<std::string> contigs)
{

    std::stringstream ss;

    seqan::SeqFileIn seqFileIn(fasta.c_str());
    seqan::StringSet<seqan::CharString> ids;
    seqan::StringSet<seqan::String<seqan::Iupac> > seqs;

    seqan::readRecords(ids, seqs, seqFileIn);

    std::vector<std::string> sequences;

    unsigned n = seqan::length(ids);

    for (unsigned i = 0; i < n; i++)
    {
        std::string id;
        move(id, ids[i]);
        if (std::find(contigs.begin(), contigs.end(), id) != contigs.end())
        {
            std::string seq;
            move(seq, seqs[i]);
            sequences.push_back(seq);
        }
    }

    return sequences;
}

void SequenceUtil::exportFilteredFasta(const std::string & fasta, const std::set<std::string> contigs, const std::string & exportFilename)
{

    std::stringstream ss;

    seqan::SeqFileIn seqFileIn(fasta.c_str());
    seqan::StringSet<seqan::CharString> ids;
    seqan::StringSet<seqan::String<seqan::Iupac> > seqs;

    seqan::readRecords(ids, seqs, seqFileIn);

    unsigned n = seqan::length(ids);

    for (unsigned i = 0; i < n; i++)
    {
        std::string id;
        move(id, ids[i]);
        if (std::find(contigs.begin(), contigs.end(), id) != contigs.end())
        {
            std::string seq;
            move(seq, seqs[i]);
            ss << '>' << id << '\n' << seq << '\n';
        }
    }

    std::ofstream ofs(exportFilename, std::ofstream::out);
    ofs << ss.str();
    ofs.close();
}

SequenceStats SequenceUtil::calculateStats(const std::string & fasta, unsigned minContigLength)
{
    SequenceStats stats;
    stats.numBasepairs = 0;
    unsigned long gc = 0;

    seqan::SeqFileIn seqFileIn(fasta.c_str());
    seqan::StringSet<seqan::CharString> ids;
    seqan::StringSet<seqan::String<seqan::Iupac> > seqs;

    seqan::readRecords(ids, seqs, seqFileIn);

    unsigned n = seqan::length(ids);

    for (unsigned i = 0; i < n; i++)
    {
        std::string id;
        move(id, ids[i]);
        std::string seq;
        move(seq, seqs[i]);

        stats.contigs.push_back(id);

        unsigned len = seqan::length(seq);
        if (len < minContigLength)
        {
            continue;
        }

        stats.numBasepairs += seq.size();
        stats.contigLength[id] = seq.size();

        unsigned long contigGc = 0;

        for (unsigned long j = 0; j < seq.size(); j++)
        {
            if (seq[j] == 'C' || seq[j] == 'G')
            {
                contigGc++;
            }
        }

        stats.contigGcContent[id] = ((double)contigGc) / ((double)seq.size());

        gc += contigGc;
    }

    stats.gcContent = ((double)gc) / ((double)stats.numBasepairs);

    return stats;
}
