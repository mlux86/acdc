#pragma once

#include <seqan/sequence.h> 
#include <Eigen/Dense>
#include <set>
#include <map>
#include <string>
#include <vector>

struct WindowRange
{
	unsigned from;
	unsigned to;
};

struct SequenceVectorizationResult
{
	Eigen::MatrixXd data;
	std::vector<std::string> contigs;
	std::vector<WindowRange> windows;
	std::map<std::string, unsigned> contigSizes;
};

class SequenceVectorizer
{
protected:
	std::string inputFasta = "";

	unsigned minContigLength = 0;
	unsigned kmerLength = 0;
	unsigned windowWidth = 0;
	unsigned windowStep = 0;
	unsigned targetNumPoints = 0;

	std::map<std::string, unsigned> kmerIndexes;
	std::set<std::string> uniqueKmers;

	bool normalize = true;

	std::vector<std::string> ids;
	std::vector<seqan::Dna5String> sequences;	

	void buildParams();
	void loadFasta();
	void estimateWindowParams();
	void buildFeatureKmers();

public:
	SequenceVectorizer(const unsigned kmerLength, const unsigned windowWidth, const unsigned windowStep);
	SequenceVectorizer(const std::string & inputFasta);
	~SequenceVectorizer();

	unsigned getDim() const;
	std::set<std::string> getFeatures() const;
	unsigned getKmerLength() const;
	unsigned getWindowWidth() const;
	unsigned getWindowStep() const;

	void setNormalize(const bool normalize_);
	bool getNormalize() const;

	std::pair<Eigen::MatrixXd, std::vector<WindowRange>> vectorize(seqan::Dna5String & sequence) const;
	SequenceVectorizationResult vectorize() const;
	
};
