#pragma once 

#include <vector>
#include "Clustering.h"
#include "Opts.h"
#include <eigen3/Eigen/Dense>

struct ClusterAnalysisResult
{
	Eigen::MatrixXd dataOrig;
	Eigen::MatrixXd dataSne;
	Eigen::MatrixXd dataPca;

	std::vector<unsigned> bootstrapIndexes;

	bool isMultiModal; // true if sne or pca multi modal
	unsigned numClustPca; // one if not multi-modal
	std::vector<ClusteringResult> clustsPca;
	unsigned numClustSne; // one if not multi-modal
	std::vector<ClusteringResult> clustsSne;

	bool hasSeparatedComponents; // true if numClustCC > 1
	unsigned numClustCC; 
	ClusteringResult clustCC;
};

class ClusterAnalysis
{

private:
	ClusterAnalysis();
	~ClusterAnalysis();
	
	static std::vector< std::vector<unsigned> > stratifiedSubsamplingIndices(const unsigned n, const unsigned k, const double ratio = 0.8);
	static ClusterAnalysisResult bootstrapTask(const Eigen::MatrixXd & dataOrig, const std::vector<std::string> & contigs, const Opts & opts, const std::vector<unsigned> indices);

public:
	static ClusterAnalysisResult analyze(const Eigen::MatrixXd & data, const std::vector<std::string> & contigs, const Opts & opts);
	static std::vector<ClusterAnalysisResult> analyzeBootstraps(const Eigen::MatrixXd & data, const std::vector<std::string> & contigs, const Opts & opts);
	
};
