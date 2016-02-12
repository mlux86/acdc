#include <vector>
#include <string>
#include <boost/filesystem.hpp>

class IOUtil
{
private:
	IOUtil();
	~IOUtil();

public:	
	static void copyDir(boost::filesystem::path const & source, boost::filesystem::path const & destination, bool overwriteExisting);
	static unsigned long long getFileSizeBytes(const std::string & filename);
	static std::vector<std::string> fileLinesToVec(const std::string & filename);
};