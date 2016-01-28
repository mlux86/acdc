#pragma once

#include <vector>
#include <string>
#include <memory>

#include <eigen3/Eigen/Dense>
#include <json/json.h>
#include "Clustering.h"

class Util
{

private:
	Util();
	~Util();
	
public:
	static void allPermsRepetition(std::vector<std::string> & perms, const std::vector<char> & alphabet, std::string elem, const unsigned length);
	static std::vector<std::string> allKmers(const unsigned kmerLength);
	static std::string reverseComplement(const std::string & seq);

	static void matrixRemoveRow(Eigen::MatrixXd & matrix, unsigned int rowToRemove);
	static void matrixRemoveColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove);
	static Eigen::MatrixXd loadMatrix(std::string filename, char delimiter = ' ');
	static void saveMatrix(const Eigen::MatrixXd & mat, std::string filename, char delimiter = ' ');

	static unsigned estimateTsnePerplexity(const Eigen::MatrixXd & mat);

	static unsigned long long getFileSizeBytes(const std::string & filename);

	static Eigen::MatrixXd knnAffinityMatrix(const Eigen::MatrixXd & data, const unsigned k, bool mutual = false);

	static Eigen::MatrixXd pca(const Eigen::MatrixXd & data, const unsigned ndims);
	static Eigen::MatrixXd pdist(const Eigen::MatrixXd & data);

	static std::vector< std::vector<unsigned> > stratifiedSubsamplingIndices(const unsigned n, const unsigned k, const double ratio = 0.8);

	static Json::Value clusteringToJson(const Eigen::MatrixXd & mat, const std::vector<unsigned> & labels, const std::vector<std::string> & tooltips);

	template<typename T>
	static std::vector<unsigned> numericLabels(const std::vector<T> & labels);

	static std::vector<std::string> fileLinesToVec(const std::string & filename);

	template<typename T>
	static std::vector<T> mutualLabels(const std::vector<T> & labels1, const std::vector<T> & labels2);

	static Eigen::MatrixXd alignDataset(const Eigen::MatrixXd & reference, const Eigen::MatrixXd & toalign, const std::vector<std::string> & labelsReference, const std::vector<std::string> & labelsToalign);

	static std::unique_ptr<std::string> filterFasta(const std::string & fasta, const std::vector<std::string> contigs);

};

#include "Util.tpp"