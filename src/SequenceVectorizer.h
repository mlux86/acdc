#ifndef __SequenceVectorizer__
#define __SequenceVectorizer__

#include <seqan/sequence.h> 
#include <eigen3/Eigen/Dense>
#include <set>
#include <map>
#include <string>

class SequenceVectorizer
{
protected:
	unsigned kmerLength;
	unsigned windowSize;
	unsigned windowStep;

	std::map<std::string, unsigned> kmerIndexes;
	std::set<std::string> uniqueKmers;

	bool normalize;

public:
	SequenceVectorizer(const unsigned kmerLength_, const unsigned windowSize_, const unsigned windowStep_);
	~SequenceVectorizer();

	unsigned getDim() const;
	std::set<std::string> getFeatures() const;
	unsigned getKmerLength() const;
	unsigned getWindowSize() const;
	unsigned getWindowStep() const;

	void setNormalize(const bool normalize_);
	bool getNormalize() const;

	Eigen::MatrixXd vectorize(seqan::Dna5String & sequence) const;
	std::pair< Eigen::MatrixXd, std::vector<std::string> > vectorize(const std::string & fasta) const;
	
};

#endif