#pragma once

#include <string>
#include <map>

// forward declaration
struct SequenceVectorizationResult; 

// Representation of an RNAmmer hit
struct RnammerResult
{
	// contig where 16S was found on
	std::string contig;

	// start position of 16S sequence on the contig
	unsigned startPos;

	// end position of 16S sequence on the contig
	unsigned endPos;
};

// Calls RNAmmer and reads its results
class RnammerAdapter
{
private:
	RnammerAdapter();
	~RnammerAdapter();
	
public:
	// Returns true if RNAmmer exists on the system and is callable
	static bool rnammerExists();

	// Calls RNAmmer and returns a map of 16S sequences
	// Key is the index of a data point, value is the corresponding 16S sequence found at that data point
	static std::map<unsigned, std::string> find16S(const std::string & fasta, const SequenceVectorizationResult & svr);

	// Extracts a 16S sequence from the given fasta and RNAmmer hit
	static std::string get16S(const std::string & fasta, const RnammerResult & rr);
	
};
