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

	static Eigen::MatrixXd loadMatrix(std::string filename, char delimiter = ' ');
	static void saveMatrix(const Eigen::MatrixXd & mat, std::string filename, char delimiter = ' ');

};


#endif