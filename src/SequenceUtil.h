#pragma once

#include <string>
#include <vector>

class SequenceUtil
{
private:
	SequenceUtil();
	~SequenceUtil();
	
public:	
	static void allPermsRepetition(std::vector<std::string> & perms, const std::vector<char> & alphabet, std::string elem, const unsigned length);
	static std::vector<std::string> allKmers(const unsigned kmerLength);
	static std::string reverseComplement(const std::string & seq);	
};