#pragma once 

#include <vector>
#include <string>
#include <map>

#include "Clustering.h"

// Forward declaration
struct SequenceVectorizationResult; 

// Contains the contamination detection result for one given data matrix 
struct ContaminationDetectionResult
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

	// Does this result have separated components?
	bool hasSeparatedComponents; 
};

struct ContaminationDetectionSummary
{
	double confidenceCC = 0;
	double confidenceDip = 0;
	std::string state;
};

// Main class to perform cluster analysis on a given sample
class ContaminationDetection
{

private:
	ContaminationDetection();
	~ContaminationDetection();
	
	// Generates a set of randomized but stratified indices for bootstrapping
	// Stratified means that it is guaranteed that each data point appears at least once in some bootstrap
	static std::vector< std::vector<unsigned> > stratifiedSubsamplingIndices(const unsigned n, const unsigned k, const double ratio = 0.8);

	// Performs cluster analysis for one bootstrap with the given indices, possibly called by some thread pool
	static ContaminationDetectionResult bootstrapTask(const SequenceVectorizationResult & svr, const std::vector<unsigned> indices, const unsigned bootstrapId);

public:
	// Performs cluster analysis for one specific data matrix and given contigs (/sizes)
	// Contigs are necessary to resolve inter-cluster contig connections which indicate a wrong clustering
	static ContaminationDetectionResult analyze(const Eigen::MatrixXd & data, const std::vector<std::string> & contigs, const std::map<std::string, unsigned> & contigSizes, const unsigned randomSeed);

	// Analyze one sequence vectorization by generating bootstraps and clustering them
	static std::vector<ContaminationDetectionResult> analyzeBootstraps(const SequenceVectorizationResult & svr);

	static ContaminationDetectionSummary summarizeBootstraps(const std::vector<ContaminationDetectionResult> & bootstraps);
	
};
