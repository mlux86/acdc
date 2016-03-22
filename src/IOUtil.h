#include <vector>
#include <string>
#include <boost/filesystem.hpp>

// IO utilities
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
};