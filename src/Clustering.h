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

	static std::pair<ClusteringResult, double> kMeansIter(const Eigen::MatrixXd & data, unsigned k);
	
public:
	static ClusteringResult spectralClustering(Eigen::MatrixXd & adjacencies);
	static ClusteringResult kMeans(const Eigen::MatrixXd & data, unsigned k, unsigned numBootstraps = 5);

};


#endif 
