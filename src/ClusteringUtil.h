#pragma once

#include <Eigen/Dense>

// Matrix handling utilities
class ClusteringUtil
{
private:
	ClusteringUtil();
	~ClusteringUtil();

public:
	static bool isMultiModal(const Eigen::MatrixXd & data, double alpha, double splitThreshold);
	static double daviesBouldin(const Eigen::MatrixXd & data, const std::vector<unsigned> & labels);
};