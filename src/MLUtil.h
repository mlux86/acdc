#pragma once

#include <Eigen/Dense>
#include <vector>

// Machine Learning Utilities
class MLUtil
{
private:
	MLUtil();
	~MLUtil();

public:
	// Returns a k-nearest neighbour affinity matrix
	// Uses a kd-Tree for fast lookup
	static Eigen::MatrixXd knnAffinityMatrix(const Eigen::MatrixXd & data, const unsigned k, bool mutual = false);

	// Computes the ndims first principal components of the given data matrix
	static Eigen::MatrixXd pca(const Eigen::MatrixXd & data, const unsigned ndims);

	// Computes the pairwise distances of data
	// Returns a quadratic matrix with the distances
	static Eigen::MatrixXd pdist(const Eigen::MatrixXd & data);

	// Computes the pairwise distances of data
	// Returns a condensed vector which contains only the upper triangular part of a pairwise quadratic distance matrix
	static Eigen::VectorXd condensedPdist(const Eigen::MatrixXd & data);

	// Performs canonical correlation analysis
	// Returns the same objects as the corresponding MATLAB implementation
	static std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::VectorXd, Eigen::MatrixXd, Eigen::MatrixXd> canonicalCorrelation(const Eigen::MatrixXd & x, const Eigen::MatrixXd & y);

	// Aligns one bootstrap data matrix to a reference using canonical correlation
	static Eigen::MatrixXd alignBootstrap(const Eigen::MatrixXd & reference, const Eigen::MatrixXd & bootstrap, const std::vector<unsigned> & bootstrapIndexes);

	// Aligns bootstrap labels to a reference to have the same labels in each corresponding cluster
	static std::vector<unsigned> alignBootstrapLabels(const std::vector<unsigned> & referenceLabels, const std::vector<unsigned> & bootstrapLabels, const std::vector<unsigned> & bootstrapIndexes);

};
