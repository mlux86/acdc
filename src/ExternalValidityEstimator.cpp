#include <Eigen/Dense>
#include <memory>

#include "HierarchicalClustering.h"
#include "ClusteringUtil.h"
#include "Clustering.h"

ExternalValidityEstimator::ExternalValidityEstimator(std::unique_ptr<ClusterPostProcessing> & cpp_, std::unique_ptr<ClusteringAlgorithm> & clustAlgo_, unsigned maxK_) : ClusterEstimator(cpp_), clustAlgo(std::move(clustAlgo_)), maxK(maxK_)
{

}

ExternalValidityEstimator::~ExternalValidityEstimator()
{

}

std::pair<unsigned, std::vector<ClusteringResult>> ExternalValidityEstimator::estimateK(const Eigen::MatrixXd & data)
{
    std::vector<ClusteringResult> results(maxK);

    // label trivial clustering for k == 1
    results[0].numClusters = 1;
    results[0].labels.resize(data.rows());
    std::fill(results[0].labels.begin(),results[0].labels.end(), 1);

    clustAlgo->initialize(data);

    // find minimum Davies Bouldin index
    double minDb = std::numeric_limits<double>::max();
    unsigned optK = 0;

    for (unsigned k = 2; k <= maxK; k++)
    {
        results[k-1].numClusters = k;
        results[k-1].labels = clustAlgo->cluster(k);        
        cpp->run(results[k-1]);

        double db = ClusteringUtil::daviesBouldin(data, results[k-1].labels);

        if (db < minDb)
        {
            minDb = db;
            optK = results[k-1].numClusters;
        }
    }

    return std::make_pair(optK, results);
}

std::string ExternalValidityEstimator::name()
{
    return "ExternalValidityEstimator";
}

std::map<std::string, std::string> ExternalValidityEstimator::parameters()
{
    std::map<std::string, std::string> params = {{"max_k", std::to_string(maxK)}, {"cluster_validity", "davies_bouldin"}};
    params["clustering_algorithm"] = clustAlgo->name();
    for (auto & it : clustAlgo->parameters())
    {
        params["clustering_" + it.first] = it.second;
    }
    return params;
}
