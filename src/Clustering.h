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



class ClusteringAlgorithm
{

public:
	ClusteringAlgorithm();
	virtual ~ClusteringAlgorithm();

	virtual void initialize(const Eigen::MatrixXd & data) = 0;
	virtual std::vector<unsigned> cluster(unsigned k) = 0;

};



class ClusterPostProcessing
{

private:
	// Contigs to use for resolving inter-cluster contig connections
	const std::vector<std::string> & contigs;

	// Contigs sizes to use for re-assigning data points with inter-cluster contig connections
	const std::map<std::string, unsigned> & contigSizes;

public:
	// Initializes members
	ClusterPostProcessing(const std::vector<std::string> & contigs_, const std::map<std::string, unsigned> & contigSizes_);
	~ClusterPostProcessing();

	// Post-processes a ClusteringResult by resolcing inter-cluster contig connections and re-labeling outlier points
	void run(ClusteringResult & cr);
};



class ClusterEstimator
{

protected:
	ClusterPostProcessing cpp;

public:
	ClusterEstimator(ClusterPostProcessing cpp_);
	~ClusterEstimator();

	virtual std::pair<unsigned, std::vector<ClusteringResult>> estimateK(const Eigen::MatrixXd & data) = 0;

};



class ConnectedComponentsEstimator : public ClusterEstimator
{

private:
	unsigned knnK;

public:
	ConnectedComponentsEstimator(ClusterPostProcessing cpp_, unsigned knnK_);
	~ConnectedComponentsEstimator();

	std::pair<unsigned, std::vector<ClusteringResult>> estimateK(const Eigen::MatrixXd & data);

};


class ExternalValidityEstimator : public ClusterEstimator
{

private:
	ClusteringAlgorithm & clustAlgo;
	unsigned maxK = 9;

public:
	ExternalValidityEstimator(ClusterPostProcessing cpp_, ClusteringAlgorithm & clustAlgo_, unsigned maxK_);
	~ExternalValidityEstimator();
	
	std::pair<unsigned, std::vector<ClusteringResult>> estimateK(const Eigen::MatrixXd & data);

};
