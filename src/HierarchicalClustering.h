#pragma once

#include <vector>
#include <Eigen/Dense>

#include "Clustering.h"

// Performs hierarchical clustering
class HierarchicalClustering : public ClusteringAlgorithm
{

private:
	// Helper function for recursively assigning the cluster number to leaf nodes in a linkage tree
	static void clusterNum(const Eigen::MatrixXd & linkage, std::vector<unsigned> & labels, unsigned idx, unsigned clsnum);

public:
	Eigen::MatrixXd linkageMat;

	HierarchicalClustering();
	~HierarchicalClustering();

	// Computes a (n-1)-by-dim linkage matrix
	// Same format as SciPys implementation of a linkage matrix
    static Eigen::MatrixXd linkage(const Eigen::MatrixXd & data);

    void initialize(const Eigen::MatrixXd & data);

    // Clusters data by cutting the linkage tree (calculated by linkage()) at level k 
	std::vector<unsigned> cluster(unsigned k);

	std::string name();
	std::map<std::string, std::string> parameters();

};
