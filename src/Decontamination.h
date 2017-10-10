#pragma once

#include <vector>
#include <Eigen/Dense>

#include "SequenceVectorizer.h"
#include "Clustering.h"

struct DecontaminationResult
{
	std::unique_ptr<ClusterEstimator> estimatorCC = nullptr;
	std::unique_ptr<ClusterEstimator> estimatorClusterValidity = nullptr;

	unsigned numClustCC; 
	std::vector<ClusteringResult> clustsCC;

	unsigned numClustSne; 
	std::vector<ClusteringResult> clustsSne;

	unsigned numClustPca; 
	std::vector<ClusteringResult> clustsPca;

	std::string mostLikelyClusteringName;
	ClusteringResult* mostLikelyClustering;
};

class Decontamination
{

private:
	Decontamination();
	~Decontamination();
	
public:
	static DecontaminationResult findLikelyClusterings(const Eigen::MatrixXd & dataSne, const Eigen::MatrixXd & dataPca, const SequenceVectorizationResult & svr);

};