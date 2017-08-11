#pragma once

#include <Eigen/Dense>
#include <vector>
#include <map>
#include <string>
#include <memory>

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

	virtual std::string name() = 0;
	virtual std::map<std::string, std::string> parameters() = 0;

};



class ClusterPostProcessing
{

private:
	// Contigs to use for resolving inter-cluster contig connections
	const std::vector<std::string> contigs;

	// Contigs sizes to use for re-assigning data points with inter-cluster contig connections
	const std::map<std::string, unsigned> contigSizes;

public:
	// Initializes members
	ClusterPostProcessing(const std::vector<std::string> contigs_, const std::map<std::string, unsigned> contigSizes_);
	~ClusterPostProcessing();

	// Post-processes a ClusteringResult by resolcing inter-cluster contig connections and re-labeling outlier points
	void run(ClusteringResult & cr);
};



class ClusterEstimator
{

protected:
	std::unique_ptr<ClusterPostProcessing> cpp;

public:
	ClusterEstimator(std::unique_ptr<ClusterPostProcessing> & cpp_);
	~ClusterEstimator();

	virtual std::pair<unsigned, std::vector<ClusteringResult>> estimateK(const Eigen::MatrixXd & data) = 0;

	virtual std::string name() = 0;
	virtual std::map<std::string, std::string> parameters() = 0;	
};



class ConnectedComponentsEstimator : public ClusterEstimator
{

protected:
	unsigned knnK;

public:
	ConnectedComponentsEstimator(std::unique_ptr<ClusterPostProcessing> & cpp_, unsigned knnK_);
	~ConnectedComponentsEstimator();

	std::pair<unsigned, std::vector<ClusteringResult>> estimateK(const Eigen::MatrixXd & data);

	std::string name();
	std::map<std::string, std::string> parameters();	
};


class ExternalValidityEstimator : public ClusterEstimator
{

protected:
	std::unique_ptr<ClusteringAlgorithm> clustAlgo;
	unsigned maxK = 5;

public:	
	ExternalValidityEstimator(std::unique_ptr<ClusterPostProcessing> & cpp_, std::unique_ptr<ClusteringAlgorithm> & clustAlgo_, unsigned maxK_);
	~ExternalValidityEstimator();

	std::pair<unsigned, std::vector<ClusteringResult>> estimateK(const Eigen::MatrixXd & data);

	std::string name();
	std::map<std::string, std::string> parameters();	
};
