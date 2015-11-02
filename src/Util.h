#ifndef __Util__
#define __Util__

#include <vector>
#include <string>

#include <eigen3/Eigen/Dense>

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
	
	template <typename T>
	static std::vector<unsigned> sortIndexes(const std::vector<T> & v) 
	{
	    // initialize original index locations
	    std::vector<unsigned> idx(v.size());
	    for (unsigned i = 0; i != idx.size(); ++i) idx[i] = i;

	    // sort indexes based on comparing values in v
	    std::sort(idx.begin(), idx.end(), [&v](unsigned i1, unsigned i2) {return v[i1] < v[i2];});

	    return idx;
	}

	static std::vector< std::vector<unsigned> > stratifiedSubsamplingIndices(const unsigned n, const unsigned k, const double ratio = 0.8);

};


#endif