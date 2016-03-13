#pragma once 

#include <eigen3/Eigen/Dense>

class MLUtil
{
private:
	MLUtil();
	~MLUtil();

public:
	static Eigen::MatrixXd knnAffinityMatrix(const Eigen::MatrixXd & data, const unsigned k, bool mutual = false);

	static Eigen::MatrixXd pca(const Eigen::MatrixXd & data, const unsigned ndims);

	static Eigen::MatrixXd pdist(const Eigen::MatrixXd & data);	

	static Eigen::VectorXd condensedPdist(const Eigen::MatrixXd & data);
	
	static std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::VectorXd, Eigen::MatrixXd, Eigen::MatrixXd> canonicalCorrelation(const Eigen::MatrixXd & x, const Eigen::MatrixXd & y);

	static Eigen::MatrixXd alignBootstrap(const Eigen::MatrixXd & reference, const Eigen::MatrixXd & bootstrap, const std::vector<unsigned> & bootstrapIndexes);

	static std::vector<unsigned> alignBootstrapLabels(const std::vector<unsigned> & referenceLabels, const std::vector<unsigned> & bootstrapLabels, const std::vector<unsigned> & bootstrapIndexes);

};