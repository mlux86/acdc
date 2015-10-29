#ifndef __Kmeans__
#define __Kmeans__

#include <eigen3/Eigen/Dense>
#include <string>

#include "Clustering.h"

class Kmeans
{

private:
	unsigned k = 1;
	unsigned numBootstraps = 1;
	
	std::string init = "plusplus";
	Eigen::MatrixXd means;

	std::pair<ClusteringResult, double> iteration(const Eigen::MatrixXd & data);

	void sampleInit(const Eigen::MatrixXd & data);
	void plusPlusInit(const Eigen::MatrixXd & data);

public:
	Kmeans(const unsigned k_);
	~Kmeans();

	void initSample();
	void initPlusPlus();
	void initMeans(const Eigen::MatrixXd & initMeans);

	Eigen::MatrixXd getMeans();	

	ClusteringResult run(const Eigen::MatrixXd & data);
	void setNumBootstraps(const unsigned numBootstraps);


};

#endif 
