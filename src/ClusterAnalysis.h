#pragma once 

#include <vector>
#include <map>

#include "Clustering.h"

// Forward declaration
struct SequenceVectorizationResult; 

// Contains the result of a cluster analysis for one given data matrix 
struct ClusterAnalysisResult
{
	// Original data before dimensionality reduction
	Eigen::MatrixXd dataOrig;

	// t-Sne reduced representation
	Eigen::MatrixXd dataSne;

	// PCA reduced representation
	Eigen::MatrixXd dataPca;

	// Bootstrap indices (empty if data is not part of a bootstrap)
	std::vector<unsigned> bootstrapIndexes;

	// True if t-Sne or PCA indicate multi modality
	bool isMultiModal; 

	// Number of clusters found by clustering PCA data. One if not multi-modal.
	unsigned numClustPca; 
	// Clustering results for different number of clusters
	std::vector<ClusteringResult> clustsPca;

	// Number of clusters found by clustering t-Sne data. One if not multi-modal.
	unsigned numClustSne; 
	// Clustering results for different number of clusters
	std::vector<ClusteringResult> clustsSne;

	// Does this result have separated components?
	bool hasSeparatedComponents; 
	// Number of connected components
	unsigned numClustCC; 
	// Connected components clustering
	ClusteringResult clustCC;
};

// Main class to perform cluster analysis on a given sample
class ClusterAnalysis
{

private:
	ClusterAnalysis();
	~ClusterAnalysis();
	
	// Generates a set of randomized but stratified indices for bootstrapping
	// Stratified means that it is guaranteed that each data point appears at least once in some bootstrap
	static std::vector< std::vector<unsigned> > stratifiedSubsamplingIndices(const unsigned n, const unsigned k, const double ratio = 0.8);

	// Performs cluster analysis for one bootstrap with the given indices, possibly called by some thread pool
	static ClusterAnalysisResult bootstrapTask(const SequenceVectorizationResult & svr, const std::vector<unsigned> indices);

public:
	// Performs cluster analysis for one specific data matrix and given contigs (/sizes)
	// Contigs are necessary to resolve inter-cluster contig connections which indicate a wrong clustering
	static ClusterAnalysisResult analyze(const Eigen::MatrixXd & data, const std::vector<std::string> & contigs,  const std::map<std::string, unsigned> & contigSizes);

	// Analyze one sequence vectorization by generating bootstraps and clustering them
	static std::vector<ClusterAnalysisResult> analyzeBootstraps(const SequenceVectorizationResult & svr);
	
};
