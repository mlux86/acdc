#include "KrakenAdapter.h"
#include "Logger.h"
#include <boost/algorithm/string.hpp>

KrakenAdapter::KrakenAdapter()
{
}

KrakenAdapter::~KrakenAdapter()
{
}

KrakenResult KrakenAdapter::runKraken(const std::string & fasta, const Opts & opts)
{
	std::string krakenCommand = opts.krakenScript() + " " + fasta;

 	FILE * file = popen(krakenCommand.c_str(), "r");

 	if (!file)
 	{
 		throw std::runtime_error("Cannot run Kraken!");
 	}

    char buffer[4096];
    fgets(buffer, sizeof(buffer), file);
    pclose(file);

    std::string krakenResultFile(buffer, strlen(buffer)-1);

    auto krakenOut = Util::fileLinesToVec(krakenResultFile);

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

    return res;
}