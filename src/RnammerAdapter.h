#pragma once

#include <string>
#include <vector>

class RnammerAdapter
{
private:
	RnammerAdapter();
	~RnammerAdapter();
	
public:
	static bool rnammerExists();
	static std::vector<std::string> find16SContigs(const std::string & fasta);
	
};
