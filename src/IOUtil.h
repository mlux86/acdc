#include <vector>
#include <string>
#include <iomanip>
#include <sstream>
#include <yaml-cpp/yaml.h>
#include <Eigen/Dense>

#define BOOST_NO_CXX11_SCOPED_ENUMS
#include <boost/filesystem.hpp>
#undef BOOST_NO_CXX11_SCOPED_ENUMS

// IO utilities

struct Fixed 
{
	Fixed(double v = 0): value(v) 
	{

	}

	double value;

	std::string toString() const 
	{
		std::stringstream ss;
		ss << std::setprecision(2) << value;
		return ss.str();
	}
};

YAML::Emitter & operator << (YAML::Emitter & out, const Fixed & f);

class IOUtil
{
private:
	IOUtil();
	~IOUtil();

public:	
	// Copies a directories contents from source to destination
	// No recursion!
	static void copyDir(boost::filesystem::path const & source, boost::filesystem::path const & destination, bool overwriteExisting);

	// Reads a files lines into a vector of strings
	static std::vector<std::string> fileLinesToVec(const std::string & filename);

	static std::vector<Fixed> columnToFixed(const Eigen::MatrixXd & m, unsigned colIdx);
};