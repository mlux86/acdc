#pragma once

#include <vector>
#include <Eigen/Dense>

// Performs hierarchical clustering
class HierarchicalClustering
{

private:
	HierarchicalClustering();
	~HierarchicalClustering();

	// Helper function for recursively assigning the cluster number to leaf nodes in a linkage tree
	static void clusterNum(const Eigen::MatrixXd & linkage, std::vector<unsigned> & labels, unsigned idx, unsigned clsnum);

public:
	// Computes a (n-1)-by-dim linkage matrix
	// Same format as SciPys implementation of a linkage matrix
    static Eigen::MatrixXd linkage(const Eigen::MatrixXd & data);

    // Clusters data by cutting the linkage tree (calculated by linkage()) at level maxK 
	static std::vector<unsigned> cluster(const Eigen::MatrixXd & linkage, unsigned maxK);

};
