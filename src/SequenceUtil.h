#pragma once

#include <string>
#include <vector>
#include <set>
#include <map>

#include "TaxonomyAnnotation.h"

struct SequenceStats
{
	// Included contigs (fulfilling all requirements, such as minimum length)
	std::vector<std::string> includedContigs;

	// Discarded contigs 
	std::vector<std::string> discardedContigs;
	
    // Size in bp
    unsigned long numBasepairs;

    // GC% content
    double gcContent;

    // Contig lengths
	std::map<std::string, unsigned long> contigLength;

    // GC% content per contig
	std::map<std::string, double> contigGcContent;

	// taxonomy information per contig
	std::map<std::string, Taxonomy> taxonomies;
};

// Sequence handling utilities
class SequenceUtil
{
private:
	SequenceUtil();
	~SequenceUtil();

public:
	// Generates all permutations with repetition of words from an alphabet with given length
	static void allPermsRepetition(std::vector<std::string> & perms, const std::vector<char> & alphabet, std::string elem, const unsigned length);

	// Generates all possible k-mers with the specified length
	static std::vector<std::string> allKmers(const unsigned kmerLength);

	// Computes the reverse complement of seq
	static std::string reverseComplement(const std::string & seq);

	// Filters a fasta file to return only some given contigs
	static std::vector<std::string> filterFasta(const std::string & fasta, const std::set<std::string> contigs);

	// Exports a filtered fasta file to include only some given contigs
	static void exportFilteredFasta(const std::string & fasta, const std::set<std::string> contigs, const std::string & exportFilename);

    // Calculate statistics on fasta
    static SequenceStats calculateStats(const std::string & fasta, unsigned minContigLength);
};
