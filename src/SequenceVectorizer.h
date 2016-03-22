#pragma once

#include <seqan/sequence.h> 
#include <Eigen/Dense>
#include <set>
#include <map>
#include <string>
#include <vector>

// Represents a vectorization window
struct WindowRange
{
	unsigned from;
	unsigned to;
};

// Represents a sequence vectorization
struct SequenceVectorizationResult
{
	// Vectorized data
	Eigen::MatrixXd data;

	// Contig labels
	std::vector<std::string> contigs;

	// Window for each data point
	std::vector<WindowRange> windows;

	// Contig sizes (in nucleotides)
	std::map<std::string, unsigned> contigSizes;
};

// Transforms DNA sequences into vectorial data
class SequenceVectorizer
{
protected:
	// Input fasta file
	std::string inputFasta = "";

	// Minimal contig length
	unsigned minContigLength = 0;

	// K-mer length
	unsigned kmerLength = 0;

	// Window width
	unsigned windowWidth = 0;

	// Window step
	unsigned windowStep = 0;

	// Target number of points
	unsigned targetNumPoints = 0;

	// Indices of k-mer dimensions
	std::map<std::string, unsigned> kmerIndexes;

	// All unique k-mers
	std::set<std::string> uniqueKmers;

	// Normalize vectorial data?
	bool normalize = true;

	// Contig IDs
	std::vector<std::string> ids;

	// Sequences
	std::vector<seqan::Dna5String> sequences;	

	// Build parameters from given Opts
	void buildParams();

	// Load fasta file
	void loadFasta();

	// Estimate window parameters if incomplete
	void estimateWindowParams();

	// Build k-mers for dimensions once
	void buildFeatureKmers();

public:
	// Use pre-computed window parameters
	SequenceVectorizer(const unsigned kmerLength, const unsigned windowWidth, const unsigned windowStep);

	// Estimate window parameters from Opts
	SequenceVectorizer(const std::string & inputFasta);

	~SequenceVectorizer();

	// Returns dimension of vectorial data
	unsigned getDim() const;

	// Returns all features as k-mers
	std::set<std::string> getFeatures() const;

	unsigned getKmerLength() const;
	unsigned getWindowWidth() const;
	unsigned getWindowStep() const;

	void setNormalize(const bool normalize_);
	bool getNormalize() const;

	// Vectorize the given sequence using all pre-computed parameters
	std::pair<Eigen::MatrixXd, std::vector<WindowRange>> vectorize(seqan::Dna5String & sequence) const;

	// Vectorize the given fasta member
	SequenceVectorizationResult vectorize() const;
	
};
