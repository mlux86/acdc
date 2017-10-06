#include "TaxonomyAnnotation.h"

#include <fstream>
#include <sstream>
#include <algorithm>
#include <map>
#include <set>

#include "ResultIO.h"

void TaxonomyAnnotation::updateContaminationState(ResultContainer & res)
{
	const auto & contigs = res.stats.includedContigs;
	const auto & optClust = res.clusterings.optimalClustering;
	auto & taxonomies = res.stats.taxonomies;

	// create per cluster label mapping and vice versa
	
	std::map<unsigned, std::set<std::string>> contigsPerCluster; 
	std::map<std::string, unsigned> clusterPerContig; 
    for (unsigned i = 0; i < optClust.labels.size(); i++)
    {
        unsigned lbl = optClust.labels.at(i);
        //TODO initialize value vector??!
    	std::string contig = res.seqVectorization.contigs.at(i);
        contigsPerCluster[lbl].insert(contig);
        clusterPerContig[contig] = lbl;
    }

    // calculate cluster sizes

	std::map<unsigned, unsigned long> clusterSizes;	
    for (const auto & it : contigsPerCluster)
    {
    	unsigned lbl = it.first;
    	unsigned long size = 0;
    	for (const auto & c : it.second)
    	{
    		size += res.stats.contigLength.at(c);
    	}
    	clusterSizes[lbl] = size;
	}    

	// calculate taxonomy frequencies in each cluster

	std::map<std::string, std::map<unsigned, double>> taxFreqs;

	// 1. iterate over contigs and get taxonomy for each
	// 2. find in which cluster contig is located
	// 3. add length of contig to found cluster

	for (auto & c : contigs)
	{
		if(taxonomies.find(c) != taxonomies.end())
		{
			std::string tax = taxonomies.at(c).identifier;
			unsigned lbl = clusterPerContig.at(c);
			taxFreqs[tax][lbl] += res.stats.contigLength.at(c) / clusterSizes.at(lbl);
		}
	}

	// calculate taxonomy confidences based on cluster

	std::map<std::string, std::map<unsigned, double>> taxFreqs;

	// 1. iterate over contigs and get taxonomy for each
	// 2. find in which cluster contig is located
	// 3. add length of contig to found cluster

	for (auto & c : contigs)
	{
		if(taxonomies.find(c) != taxonomies.end())
		{
			std::string tax = taxonomies.at(c).identifier;
			unsigned lbl = clusterPerContig.at(c);
			taxFreqs[tax][lbl] += res.stats.contigLength.at(c) / clusterSizes.at(lbl);
		}
	}


}

void TaxonomyAnnotation::annotateFromFile(ResultContainer & res, const std::string & filename)
{
	// line format: contig id<tab>taxonomy

	std::ifstream infile(filename);
	std::string line;

	const auto & contigs = res.stats.includedContigs;

	while(std::getline(infile, line))
	{
		std::stringstream linestream(line);
		std::string contig;
		std::string taxonomy;
		std::getline(linestream, contig, '\t');
		std::getline(linestream, taxonomy, '\t');
		std::transform(taxonomy.begin(), taxonomy.end(), taxonomy.begin(), ::tolower);

		if(std::find(contigs.begin(), contigs.end(), contig) != contigs.end())
		{
			res.stats.taxonomies[contig].identifier = taxonomy;
		}
	}

	updateContaminationState(res);
}

void TaxonomyAnnotation::annotateContig(ResultContainer & res, const std::string & contigId, const std::string & taxonomy)
{	
	std::string t = taxonomy;
	std::transform(t.begin(), t.end(), t.begin(), ::tolower);
	res.stats.taxonomies[contigId].identifier = t;

	updateContaminationState(res);
}

void TaxonomyAnnotation::annotateUnknown(ResultContainer & res)
{

}