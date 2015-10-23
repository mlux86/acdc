#ifndef __Clustering__
#define __Clustering__

#include <eigen3/Eigen/Dense>

struct ClusteringResult
{
	unsigned numClusters = 0;
	std::vector<unsigned> labels;
};

class Clustering
{

private:
	Clustering();
	~Clustering();
	
public:
	static ClusteringResult spectralClustering(Eigen::MatrixXd & adjacencies);

};


#endif 
