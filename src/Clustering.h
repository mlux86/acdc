#pragma once

#include <Eigen/Dense>
#include <vector>
#include <map>

// Contains result of the clustering on some (unconnected) data matrix
struct ClusteringResult
{
	// (Estimated) number of clusters
	unsigned numClusters = 0; 

	// Labels of this clustering
	std::vector<unsigned> labels;

	// Eventual outlier clusters, contains outlier cluster labels
	std::vector<unsigned> outlierClusters;
};

// Main class to perform clustering 
class Clustering
{

private:
	// Data to cluster on
	const Eigen::MatrixXd & data;

	// Contigs to use for resolving inter-cluster contig connections
	const std::vector<std::string> & contigs;

	// Contigs sizes to use for re-assigning data points with inter-cluster contig connections
	const std::map<std::string, unsigned> & contigSizes;

public:
	// Initializes members
	Clustering(const Eigen::MatrixXd & data_, const std::vector<std::string> & contigs_, const std::map<std::string, unsigned> & contigSizes_);
	~Clustering();

	// Runs connected component clustering on a mutual nearest neighbor graph with knnK neighbors
	ClusteringResult connComponents(unsigned knnK);

	// Estimates the optimal number of clusters and returns a ClusteringResult for each of the tested numbers [1, maxK]
	// The optimal number of clusters is found by looking at the minimum Davies Bouldin index
	std::pair<unsigned, std::vector<ClusteringResult>> estimateK(unsigned maxK);

	// Post-processes a ClusteringResult by resolcing inter-cluster contig connections and re-labeling outlier points
	void postprocess(ClusteringResult & cr);
};
