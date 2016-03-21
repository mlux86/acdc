#pragma once

#include <Eigen/Dense>
#include <vector>
#include <map>

#include "Opts.h"

struct ClusteringResult
{
	unsigned numClusters = 0;
	std::vector<unsigned> labels;
	std::vector<unsigned> outlierClusters;
};

class Clustering
{

private:
	const Opts & opts;
	const Eigen::MatrixXd & data;
	const std::vector<std::string> & contigs;
	const std::map<std::string, unsigned> & contigSizes;

public:
	Clustering(const Opts & opts_, const Eigen::MatrixXd & data_, const std::vector<std::string> & contigs_, const std::map<std::string, unsigned> & contigSizes_);
	~Clustering();

	ClusteringResult connComponents(unsigned knnK);
	bool isMultiModal(double alpha, double splitThreshold);
	std::pair<unsigned, std::vector<ClusteringResult>> estimateK(unsigned maxK);

	void postprocess(ClusteringResult & cr);
	double daviesBouldin(const std::vector<unsigned> & labels);
};
