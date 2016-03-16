#pragma once

#include <eigen3/Eigen/Dense>
#include <vector>

struct ClusteringResult
{
	unsigned numClusters = 0;
	std::vector<unsigned> labels;
};

class Clustering
{

private:
	Clustering();
	~Clustering();

public:
	static ClusteringResult connComponents(const Eigen::MatrixXd & data, unsigned knnK, const std::vector<std::string> & contigs);
	static bool isMultiModal(const Eigen::MatrixXd & data, double alpha, double splitThreshold);
	static double daviesBouldin(const Eigen::MatrixXd & data, const std::vector<unsigned> & labels);
	static std::pair<unsigned, std::vector<ClusteringResult>> estimateK(const Eigen::MatrixXd & data, unsigned maxK, const std::vector<std::string> & contigs);
	static void postprocess(ClusteringResult & cr, const std::vector<std::string> & contigs);
};
