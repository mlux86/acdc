#pragma once

#include <seqan/sequence.h> 
#include <eigen3/Eigen/Dense>
#include <set>
#include <map>
#include <string>
#include <vector>
#include "Opts.h"

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

	void buildParams(const Opts & opts);
	void loadFasta();
	void estimateWindowParams();
	void buildFeatureKmers();

public:
	SequenceVectorizer(const unsigned kmerLength, const unsigned windowWidth, const unsigned windowStep);
	SequenceVectorizer(const std::string & inputFasta, const Opts & opts);
	~SequenceVectorizer();

	unsigned getDim() const;
	std::set<std::string> getFeatures() const;
	unsigned getKmerLength() const;
	unsigned getWindowWidth() const;
	unsigned getWindowStep() const;

	void setNormalize(const bool normalize_);
	bool getNormalize() const;

	Eigen::MatrixXd vectorize(seqan::Dna5String & sequence) const;
	std::pair< Eigen::MatrixXd, std::vector<std::string> > vectorize() const;
	
};
