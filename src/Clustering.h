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
	// static ClusteringResult spectralClustering(Eigen::MatrixXd & adjacencies);
    static ClusteringResult dipMeans(const Eigen::MatrixXd& data, double alpha = 0.0, double splitThreshold = 0.01, unsigned maxClusters = 0);    

};
