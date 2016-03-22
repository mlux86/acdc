#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <cstdio>
#include <string>

#include "KrakenAdapter.h"
#include "Logger.h"
#include "IOUtil.h"
#include "Opts.h"

KrakenAdapter::KrakenAdapter()
{
}

KrakenAdapter::~KrakenAdapter()
{
}

bool KrakenAdapter::krakenExists()
{
    if (Opts::krakenDb().empty())
    {
        return false;
    }

    // Try to call kraken and kraken-translate executables and check their return values

    std::string krakenCommand = "kraken -h > /dev/null 2>&1";
    FILE * f1 = popen(krakenCommand.c_str(), "r");
    if (!f1)
    {
        throw std::runtime_error("fork() or pipe() failed!");
    }
    int r1 = pclose(f1);
    if (WEXITSTATUS(r1) != 0)
    {
        return false;
    }


    std::string krakenTranslateCommand = "kraken-translate -h > /dev/null 2>&1";
    FILE * f2 = popen(krakenTranslateCommand.c_str(), "r");
    if (!f2)
    {
        throw std::runtime_error("fork() or pipe() failed!");
    }
    int r2 = pclose(f2);
    if (WEXITSTATUS(r2) != 0)
    {
        return false;
    }    

    // No error, Kraken seems to work
    return true;
}

KrakenResult KrakenAdapter::runKraken(const std::string & fasta)
{
    // use temporary files for Kraken results
    boost::filesystem::path temp = boost::filesystem::unique_path();
    const std::string fname = temp.native(); 
    const std::string fnameT = fname + ".t";

    // generate shell commands
	std::string krakenCommand = "kraken --db '" + Opts::krakenDb() + "' --output '" + fname + "' '" + fasta + "'";
    if (!(Logger::getInstance().getLevel() == Verbose || Logger::getInstance().getLevel() == Debug))
    {
        krakenCommand = krakenCommand + " > /dev/null 2>&1";
    }
    std::string krakenTranslateCommand = "kraken-translate --db '" + Opts::krakenDb() + "' '" + fname + "' > '" + fnameT + "'";

    // execute kraken
    DLOG << "Executing: " << krakenCommand << "\n";
 	FILE * f1 = popen(krakenCommand.c_str(), "r");
 	if (!f1)
 	{
 		throw std::runtime_error("fork() or pipe() failed!");
 	}
    int r1 = pclose(f1);
    if (WEXITSTATUS(r1) != 0)
    {
        throw std::runtime_error("Kraken finished abnormally. Use -vv or -vvv switch to show more information.");    
    }

    // execute kraken-translate
    DLOG << "Executing: " << krakenTranslateCommand << "\n";
    FILE * f2 = popen(krakenTranslateCommand.c_str(), "r");
    if (!f2)
    {
        throw std::runtime_error("fork() or pipe() failed!");
    }
    int r2 = pclose(f2);
    if (WEXITSTATUS(r2) != 0)
    {
        throw std::runtime_error("Kraken-translate finished abnormally. Use -vv or -vvv switch to show more information.");    
    }

    // parse kraken-translate output
    auto krakenOut = IOUtil::fileLinesToVec(fnameT);
    std::remove(fname.c_str());
    std::remove(fnameT.c_str());

    KrakenResult res;

    for (const auto & line : krakenOut) // each line is one classification of a contig
    {
        std::vector<std::string> parts;
        boost::split(parts, line, boost::is_any_of("\t"));

        std::string contig = parts[0];

        boost::split(parts, parts[1], boost::is_any_of(";"));

        // split phylogenentic classification to find species
        std::string species;
        if (parts.size() >= 9)
        {
            species = parts[8];
        } else
        {
            species = "unknown";
        }
        res.classification[contig] = species;

        // split phylogenentic classification to find domain (for bacterial background)
        std::string domain;
        if (parts.size() >= 3)
        {
            domain = parts[2];
        } else
        {
            domain = "unknown";
        }
        if (domain == "Bacteria")
        {
            res.bacterialBackground++;
        }

    }

    res.bacterialBackground /= krakenOut.size();

    return res;
}