#pragma once

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


	void sampleInit(const Eigen::MatrixXd & data);
	void plusPlusInit(const Eigen::MatrixXd & data);

public:
	Kmeans(const unsigned k_);
	~Kmeans();

	void setNumBootstraps(const unsigned numBootstraps);
	Eigen::MatrixXd getMeans();	
	
	void initSample();
	void initPlusPlus();
	void initMeans(const Eigen::MatrixXd & initMeans);

	std::pair<ClusteringResult, double> iteration(const Eigen::MatrixXd & data);
	ClusteringResult run(const Eigen::MatrixXd & data);


};
