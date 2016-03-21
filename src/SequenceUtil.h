#pragma once

#include <string>
#include <vector>
#include <set>

class SequenceUtil
{
private:
	SequenceUtil();
	~SequenceUtil();
	
public:	
	static void allPermsRepetition(std::vector<std::string> & perms, const std::vector<char> & alphabet, std::string elem, const unsigned length);
	static std::vector<std::string> allKmers(const unsigned kmerLength);
	static std::string reverseComplement(const std::string & seq);	
	static std::vector<std::string> filterFasta(const std::string & fasta, const std::set<std::string> contigs);
	static void exportFilteredFasta(const std::string & fasta, const std::set<std::string> contigs, const std::string & exportFilename);
};