#pragma once

#include <string>
#include <map>
#include "Clustering.h"

struct KrakenResult
{
	double bacterialBackground = 0.0;
	std::map<std::string, std::string> classification;
};

class KrakenAdapter
{

private:
	
public:
	KrakenAdapter();
	~KrakenAdapter();
	
	bool krakenExists();
	KrakenResult runKraken(const std::string & fasta);
	
};
