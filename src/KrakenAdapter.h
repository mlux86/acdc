#pragma once

#include <string>
#include <map>

// Result container for Kraken
struct KrakenResult
{
	// Bacterial background as reported by Kraken
	double bacterialBackground = 0.0;

	// Species as reported by Kraken in the form of a map [Contig => Species]
	std::map<std::string, std::string> classification;
};

// Run Kraken via this Adapter
class KrakenAdapter
{

private:
	
public:
	KrakenAdapter();
	~KrakenAdapter();
	
	// Returns true if kraken can be called on this system
	bool krakenExists();

	// Runs kraken on the specified fasta file
	KrakenResult runKraken(const std::string & fasta);
	
};
