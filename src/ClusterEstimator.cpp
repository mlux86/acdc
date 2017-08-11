#include "Clustering.h"

ClusterEstimator::ClusterEstimator(std::unique_ptr<ClusterPostProcessing> & cpp_) : cpp(std::move(cpp_))
{

}

ClusterEstimator::~ClusterEstimator()
{

}