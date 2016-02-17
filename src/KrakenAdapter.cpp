#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <cstdio>
#include <string>

#include "KrakenAdapter.h"
#include "Logger.h"
#include "IOUtil.h"

KrakenAdapter::KrakenAdapter()
{
}

KrakenAdapter::~KrakenAdapter()
{
}

bool KrakenAdapter::krakenExists()
{
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

    return true;
}

KrakenResult KrakenAdapter::runKraken(const std::string & fasta, const Opts & opts)
{
    boost::filesystem::path temp = boost::filesystem::unique_path();
    const std::string fname = temp.native(); 
    const std::string fnameT = fname + ".t";

	std::string krakenCommand = "kraken --db '" + opts.krakenDb() + "' --output '" + fname + "' '" + fasta + "'";
    if (!(Logger::getInstance().getLevel() == Verbose || Logger::getInstance().getLevel() == Debug))
    {
        krakenCommand = krakenCommand + " > /dev/null 2>&1";
    }
    std::string krakenTranslateCommand = "kraken-translate --db '" + opts.krakenDb() + "' '" + fname + "' > '" + fnameT + "'";

    DLOG << "Executing: " << krakenCommand << "\n";
 	FILE * f1 = popen(krakenCommand.c_str(), "r");
 	if (!f1)
 	{
 		throw std::runtime_error("fork() or pipe() failed!");
 	}
    int r1 = pclose(f1);
    if (WEXITSTATUS(r1) != 0)
    {
        throw std::runtime_error("kraken finished abnormally.");    
    }

    DLOG << "Executing: " << krakenTranslateCommand << "\n";
    FILE * f2 = popen(krakenTranslateCommand.c_str(), "r");
    if (!f2)
    {
        throw std::runtime_error("fork() or pipe() failed!");
    }
    int r2 = pclose(f2);
    if (WEXITSTATUS(r2) != 0)
    {
        throw std::runtime_error("kraken-translate finished abnormally.");    
    }

    auto krakenOut = IOUtil::fileLinesToVec(fnameT);

    KrakenResult res;

    for (const auto & line : krakenOut)
    {
    	std::vector<std::string> parts;
    	boost::split(parts, line, boost::is_any_of("\t"));

    	std::string contig = parts[0];

		boost::split(parts, parts[1], boost::is_any_of(";"));

		std::string species;
		if (parts.size() >= 9)
		{
			species = parts[8];
		} else
		{
			species = "unknown";
		}

		res.classification[contig] = species;
    }

    std::remove(fname.c_str());
    std::remove(fnameT.c_str());

    return res;
}