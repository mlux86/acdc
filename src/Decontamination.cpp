#include "HierarchicalClustering.h"
#include "Decontamination.h"

Decontamination::Decontamination()
{

}

Decontamination::~Decontamination()
{

}

DecontaminationResult Decontamination::findLikelyClusterings(const Eigen::MatrixXd & dataSne, const Eigen::MatrixXd & dataPca, const SequenceVectorizationResult & svr)
{
	DecontaminationResult res;

	ClusterPostProcessing cpp(svr.contigs, svr.contigSizes);
	HierarchicalClustering hc;
	ExternalValidityEstimator eve(cpp, hc, 5);

	auto resPca = eve.estimateK(dataPca);
	res.numClustPca = resPca.first;
	res.clustsPca = resPca.second;

	auto resSne = eve.estimateK(dataSne);
	res.numClustSne = resSne.first;
	res.clustsSne = resSne.second;

	ConnectedComponentsEstimator cce(cpp, 9);
	auto resCC = cce.estimateK(dataSne);
	res.numClustCC = resCC.first;
	res.clustCC = resCC.second.at(0);

	return(res);
}