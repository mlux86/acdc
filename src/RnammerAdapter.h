#pragma once

#include <string>
#include <vector>

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

	// Calls RNAmmer and returns a vector vec of 16S sequences
	// vec[i].empty() => no 16S gene found in oneshot data point i
	// otherwise       => 16S sequence vec[i] found in oneshot data point i	
	static std::vector<std::string> find16S(const std::string & fasta, const SequenceVectorizationResult & svr);

	// Extracts a 16S sequence from the given fasta and RNAmmer hit
	static std::string get16S(const std::string & fasta, const RnammerResult & rr);
	
};
