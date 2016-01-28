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
	KrakenAdapter();
	~KrakenAdapter();

public:
	static KrakenResult runKraken(const std::string & fasta, const Opts & opts);
	
};
