#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <cstdio>
#include <string>

#include "RnammerAdapter.h"
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

std::vector<std::string> RnammerAdapter::find16SContigs(const std::string & fasta)
{
    boost::filesystem::path temp = boost::filesystem::unique_path();
    const std::string fname = temp.native();

	std::string rnammerCommand = "rnammer -s bacterial -m lsu,ssu,tsu -gff '" + fname + "' '" + fasta + "' > /dev/null 2>&1";

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

    auto rnammerOut = IOUtil::fileLinesToVec(fname);
    std::remove(fname.c_str());

    std::vector<std::string> result;

    for (const auto & line : rnammerOut)
    {
        std::vector<std::string> parts;
        boost::split(parts, line, boost::is_any_of("\t"));

        if (parts.size() < 10)
        {
            continue;
        }

        std::string contig = parts[0];
        std::string attr = parts[8];

        if (attr == "16s_rRNA")
        {
            result.push_back(contig);
            VLOG << "[rnammer] Found 16S sequence in contig " << contig << std::endl;
        }
    }

    return result;
}