#include <fstream>

#include "IOUtil.h"

unsigned long long IOUtil::getFileSizeBytes(const std::string & filename)
{
    std::ifstream in(filename, std::ifstream::ate | std::ifstream::binary);
    return in.tellg(); 
}

void IOUtil::copyDir(boost::filesystem::path const & source, boost::filesystem::path const & destination)
{
    namespace fs = boost::filesystem;

    // Check whether the function call is valid
    if(!fs::exists(source) || !fs::is_directory(source))
    {
        throw std::runtime_error("Source directory does not exist or is not a directory.");
    }

    // Iterate through the source directory
    for(fs::directory_iterator file(source); file != fs::directory_iterator(); ++file)
    {
        fs::path current(file->path());
        if(!fs::is_directory(current))
        {
            fs::copy_file(current, destination / current.filename());
        }
    }
}

std::vector<std::string> IOUtil::fileLinesToVec(const std::string & filename)
{
    std::ifstream ifs(filename);
    std::string str;
    std::vector<std::string> result;
    while (std::getline(ifs, str))
    {
        if(!str.empty())
        {
            result.push_back(str);
        }
    }     
    ifs.close();
    return result;
}