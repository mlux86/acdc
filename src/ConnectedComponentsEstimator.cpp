#include <Eigen/Dense>
#include <string>

#include "TarjansAlgorithm.h"
#include "MLUtil.h"
#include "Clustering.h"

ConnectedComponentsEstimator::ConnectedComponentsEstimator(std::unique_ptr<ClusterPostProcessing> & cpp_, unsigned knnK_) : ClusterEstimator(cpp_), knnK(knnK_)
{

}

ConnectedComponentsEstimator::~ConnectedComponentsEstimator()
{

}

std::pair<unsigned, std::vector<ClusteringResult>> ConnectedComponentsEstimator::estimateK(const Eigen::MatrixXd & data)
{
	TarjansAlgorithm ta;
    Eigen::MatrixXd aff = MLUtil::knnAffinityMatrix(data, knnK, true);
    auto res = ta.run(aff);
    cpp->run(res); 

    std::vector<ClusteringResult> clusterings;
    clusterings.push_back(res);

    return std::make_pair(res.numClusters, clusterings);
}

std::string ConnectedComponentsEstimator::name()
{
	return "ConnectedComponentsEstimator";
}

std::map<std::string, std::string> ConnectedComponentsEstimator::parameters()
{
    std::map<std::string, std::string> params = {{"nn_graph_k", std::to_string(knnK)}};
    return params;
}
