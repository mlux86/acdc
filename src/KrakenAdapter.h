#pragma once

#include <string>
#include <map>
#include "Opts.h"
#include "Clustering.h"

struct KrakenResult
{
	std::map<std::string, std::string> classification;
};

class KrakenAdapter
{

private:
	const Opts & opts;
	
public:
	KrakenAdapter(const Opts & opts);
	~KrakenAdapter();
	
	bool krakenExists();
	KrakenResult runKraken(const std::string & fasta);
	
};
