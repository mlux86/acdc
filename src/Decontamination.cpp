#include <memory>

#include "HierarchicalClustering.h"
#include "Decontamination.h"
#include "Logger.h"

Decontamination::Decontamination()
{

}

Decontamination::~Decontamination()
{

}

DecontaminationResult Decontamination::findLikelyClusterings(const Eigen::MatrixXd & dataSne, const Eigen::MatrixXd & dataPca, const SequenceVectorizationResult & svr)
{
	DecontaminationResult res;

	VLOG << "Initializing clustering estimators and post-processing" << std::endl;

	std::unique_ptr<ClusterPostProcessing> cpp1(new ClusterPostProcessing(svr.contigs, svr.contigSizes));
	std::unique_ptr<ClusteringAlgorithm> hc(new HierarchicalClustering());
	res.estimatorClusterValidity = std::unique_ptr<ClusterEstimator>(new ExternalValidityEstimator(cpp1, hc, 5));

	std::unique_ptr<ClusterPostProcessing> cpp2(new ClusterPostProcessing(svr.contigs, svr.contigSizes));
	res.estimatorCC = std::unique_ptr<ClusterEstimator>(new ConnectedComponentsEstimator(cpp2, 9));

	VLOG << "Running cluster validity estimator on PCA data" << std::endl;
	auto resPca = res.estimatorClusterValidity->estimateK(dataPca);
	res.numClustPca = resPca.first;
	res.clustsPca = resPca.second;

	VLOG << "Running cluster validity estimator on t-SNE data" << std::endl;
	auto resSne = res.estimatorClusterValidity->estimateK(dataSne);
	res.numClustSne = resSne.first;
	res.clustsSne = resSne.second;

	VLOG << "Running connected components estimator on t-SNE data" << std::endl;
	auto resCC = res.estimatorCC->estimateK(dataSne);
	res.numClustCC = resCC.first;
	res.clustsCC = resCC.second;

	return(res);
}