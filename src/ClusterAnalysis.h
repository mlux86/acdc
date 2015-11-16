#pragma once 

#include <vector>
#include "Clustering.h"
#include "Opts.h"
#include <eigen3/Eigen/Dense>

struct ClusterAnalysisResult
{
	ClusteringResult resConnComponents;
	ClusteringResult resDipMeans;
};

class ClusterAnalysis
{

private:
	ClusterAnalysis();
	~ClusterAnalysis();
	static ClusterAnalysisResult bootstrapTask(const unsigned taskId, const Eigen::MatrixXd & data, const Opts & opts, const std::vector<unsigned> indices);

public:
	static ClusterAnalysisResult analyze(const Eigen::MatrixXd & data, const Opts & opts, const std::string & descPrefix);
	static std::vector<ClusterAnalysisResult> analyzeBootstraps(const Eigen::MatrixXd & data, const Opts & opts);
	
};
