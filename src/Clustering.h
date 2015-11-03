#ifndef __Clustering__
#define __Clustering__

#include <eigen3/Eigen/Dense>

struct ClusteringResult
{
	unsigned numClusters = 0;
	Eigen::VectorXd labels;
};

class Clustering
{

private:
	Clustering();
	~Clustering();

public:
	static ClusteringResult spectralClustering(Eigen::MatrixXd & adjacencies);
    static ClusteringResult dipMeans(const Eigen::MatrixXd& data, double alpha = 0.0, double splitThreshold = 0.01, unsigned maxClusters = 0);    

};


#endif 
