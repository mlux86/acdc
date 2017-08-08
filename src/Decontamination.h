#pragma once

#include <vector>
#include <Eigen/Dense>

#include "SequenceVectorizer.h"
#include "Clustering.h"

struct DecontaminationResult
{
	unsigned numClustCC; 
	ClusteringResult clustCC;

	unsigned numClustSne; 
	std::vector<ClusteringResult> clustsSne;

	unsigned numClustPca; 
	std::vector<ClusteringResult> clustsPca;	
};

class Decontamination
{

private:
	Decontamination();
	~Decontamination();
	
public:
	static DecontaminationResult findLikelyClusterings(const Eigen::MatrixXd & dataSne, const Eigen::MatrixXd & dataPca, const SequenceVectorizationResult & svr);

};