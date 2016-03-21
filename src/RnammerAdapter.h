#pragma once

#include <string>
#include <vector>

#include "Opts.h"
#include "SequenceVectorizer.h"

struct RnammerResult
{
	std::string contig;
	unsigned startPos;
	unsigned endPos;
};

class RnammerAdapter
{
private:
	RnammerAdapter();
	~RnammerAdapter();
	
public:
	static bool rnammerExists();
	static std::vector<std::string> find16S(const std::string & fasta, const SequenceVectorizationResult & svr);
	static std::string get16S(const std::string & fasta, const RnammerResult & rr);
	
};
