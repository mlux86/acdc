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

	ClusteringResult resConnComponents;
	ClusteringResult resDipMeans;
};

class ClusterAnalysis
{

private:
	ClusterAnalysis();
	~ClusterAnalysis();
	
	static std::vector< std::vector<unsigned> > stratifiedSubsamplingIndices(const unsigned n, const unsigned k, const double ratio = 0.8);
	static ClusterAnalysisResult bootstrapTask(const Eigen::MatrixXd & dataOrig, const Opts & opts, const std::vector<unsigned> indices);

public:
	static ClusterAnalysisResult analyze(const Eigen::MatrixXd & data, const Opts & opts);
	static std::vector<ClusterAnalysisResult> analyzeBootstraps(const Eigen::MatrixXd & data, const Opts & opts);
	
};
