#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <cstdio>
#include <string>
#include <sstream>
#include <set>

#include "RnammerAdapter.h"
#include "SequenceUtil.h"
#include "SequenceVectorizer.h"
#include "Logger.h"
#include "IOUtil.h"

RnammerAdapter::RnammerAdapter()
{
}

RnammerAdapter::~RnammerAdapter()
{
}

bool RnammerAdapter::rnammerExists()
{
    std::string rnammerCommand = "rnammer -v > /dev/null 2>&1";
    FILE * f = popen(rnammerCommand.c_str(), "r");
    if (!f)
    {
        throw std::runtime_error("fork() or pipe() failed!");
    }
    int r = pclose(f);
    if (WEXITSTATUS(r) != 0)
    {
        return false;
    }
    return true;
}

std::map<unsigned, std::string> RnammerAdapter::find16S(const std::string & fasta, const SequenceVectorizationResult & svr)
{
    // Temporary files
    boost::filesystem::path temp = boost::filesystem::unique_path();
    const std::string fname = temp.native();

    // Run RNAmmer

	std::string rnammerCommand = "rnammer -s bac,arch -m lsu,ssu,tsu -gff '" + fname + "' '" + fasta + "' > /dev/null 2>&1";
    DLOG << "Executing: " << rnammerCommand << "\n";
 	FILE * f = popen(rnammerCommand.c_str(), "r");
 	if (!f)
 	{
 		throw std::runtime_error("fork() or pipe() failed!");
 	}
    int r = pclose(f);
    if (WEXITSTATUS(r) != 0)
    {
        throw std::runtime_error("Rnammer finished abnormally. Use -vv or -vvv switch to show more information.");    
    }

    // Parse output

    auto rnammerOut = IOUtil::fileLinesToVec(fname);
    std::remove(fname.c_str());

    std::vector<RnammerResult> result;

    for (const auto & line : rnammerOut)
    {
        std::vector<std::string> parts;
        boost::split(parts, line, boost::is_any_of("\t"));

        if (parts.size() < 10)
        {
            continue;
        }

        std::string contig = parts[0];
        unsigned startPos = boost::lexical_cast<unsigned>(parts[3]);
        unsigned endPos = boost::lexical_cast<unsigned>(parts[4]); 
        std::string attr = parts[8];

        if (attr == "16s_rRNA") 
        {
            VLOG << "[rnammer] Found 16S sequence in contig " << contig << std::endl;
            RnammerResult res;
            res.contig = contig;
            res.startPos = startPos;
            res.endPos = endPos;
            result.push_back(res);
        }
    }

    unsigned n = svr.contigs.size();
    std::map<unsigned, std::string> _16s;
    std::stringstream ss;

    // Process hits and assemble result vector with 16S sequences

    std::set<unsigned> processedResults;

    for (unsigned i = 0; i < n; i++)
    {
        const std::string & contig = svr.contigs.at(i);
        auto win = svr.windows.at(i);

        for (unsigned j = 0; j < result.size(); j++)            
        {
            const auto & rr = result.at(j);
            if ((rr.contig == contig) && (processedResults.find(j) == processedResults.end()) && ( 
                                        (win.from >= rr.startPos && win.from <= rr.endPos) || //overlap left
                                        (win.to >= rr.startPos && win.to <= rr.endPos) || //overlap right
                                        (rr.startPos >= win.from && rr.endPos <= win.to))) //contained
            {
               _16s[i] = get16S(fasta, rr);
               processedResults.insert(j);
            }
        }
    }

    return _16s;
}

std::string RnammerAdapter::get16S(const std::string & fasta, const RnammerResult & rr)
{
    std::set<std::string> contigs = { rr.contig };
    auto fRes = SequenceUtil::filterFasta(fasta, contigs);
    if (fRes.size() != 1)
    {
        throw std::runtime_error("[export16S] Couldn't find unique contig.");
    }
    std::string seq = fRes[0];
    return seq.substr(rr.startPos, rr.endPos-rr.startPos);
}