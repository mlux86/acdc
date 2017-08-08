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

	Clustering cSne(dataSne, svr.contigs, svr.contigSizes);
	Clustering cPca(dataPca, svr.contigs, svr.contigSizes);

	auto resPca = cPca.estimateK(5);
	res.numClustPca = resPca.first;
	res.clustsPca = resPca.second;

	auto resSne = cSne.estimateK(5);
	res.numClustSne = resSne.first;
	res.clustsSne = resSne.second;

	res.clustCC = cSne.connComponents(9);
	res.numClustCC = res.clustCC.numClusters;

	return(res);
}