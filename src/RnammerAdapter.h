#pragma once

#include <string>
#include <vector>

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
	static std::vector<bool> mark16S(const std::string & fasta, const SequenceVectorizationResult & svr);
	
};
