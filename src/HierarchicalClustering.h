#pragma once

#include "Clustering.h"

#include <vector>

class HierarchicalClustering
{

private:
	HierarchicalClustering();
	~HierarchicalClustering();
	static void clusterNum(const Eigen::MatrixXd & linkage, std::vector<unsigned> & labels, unsigned idx, unsigned clsnum);

public:
    static Eigen::MatrixXd linkage(const Eigen::MatrixXd & data);
	static std::vector<unsigned> cluster(const Eigen::MatrixXd & linkage, unsigned maxK);

};
