#include "TaxonomyAnnotation.h"

#include <fstream>
#include <sstream>
#include <algorithm>
#include <map>
#include <set>
#include <boost/algorithm/string/predicate.hpp>
#include <queue>
#include <regex>

#include "ResultIO.h"
#include "IOUtil.h"
#include "Opts.h"
#include "Logger.h"

YAML::Emitter & operator << (YAML::Emitter & out, const Taxonomy & t) 
{    
	std::string stateStr;
	switch (t.state)
	{
		case Clean: stateStr = "clean"; break;
		case Contamination: stateStr = "contamination"; break;
		default: stateStr = "NA"; break;
	}

	std::string typeStr;
	switch (t.type)
	{
		case Annotated: typeStr = "annotated"; break;
		default: typeStr = "estimated"; break;
	}	
	
	out << YAML::Flow << YAML::BeginMap;
    out << YAML::Key << "type" << YAML::Value << typeStr;
    out << YAML::Key << "state" << YAML::Value << stateStr;
    out << YAML::Key << "identifier" << YAML::Value << t.identifier;
    out << YAML::Key << "confidence" << YAML::Value << Fixed(t.confidence);
    out << YAML::EndMap;

    return out;
}

void TaxonomyAnnotation::updateContaminationConfidences(ResultContainer & res)
{
	const auto & contigs = res.stats.includedContigs;
	const auto & optClust = *(res.clusterings.mostLikelyClustering);
	auto & taxonomies = res.stats.taxonomies;

	// create per cluster label mapping and vice versa
	
	std::map<unsigned, std::set<std::string>> contigsPerCluster; 
	std::map<std::string, unsigned> clusterPerContig; 
    for (unsigned i = 0; i < optClust.labels.size(); i++)
    {
        unsigned lbl = optClust.labels.at(i);
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
			taxFreqs[tax][lbl] += (double)res.stats.contigLength.at(c) / (double)clusterSizes.at(lbl);
		}
	}

	// calculate taxonomy confidences based on cluster

	for (auto & it : taxFreqs)
	{
		auto & freqs = it.second;

		double freqSum = 0;
		for (const auto & it2 : freqs)
		{
			freqSum += it2.second;
		}

		for (auto & it2 : freqs) 
		{
			it2.second = 0.5 * it2.second + 0.5 * it2.second / freqSum;
		};
	}

	// assign confidences to result

	for (const auto & it : contigsPerCluster)
	{
		unsigned lbl = it.first;
		auto contigs = it.second;
		for (auto & c : contigs)
		{
			auto & t = res.stats.taxonomies[c];
			if (t.identifier != "unknown")
			{
				t.confidence = taxFreqs[t.identifier][lbl];
			}
		}
	}
}

void TaxonomyAnnotation::updateContaminationState(ResultContainer & res)
{	
	std::string target = Opts::targetTaxonomy();
	std::transform(target.begin(), target.end(), target.begin(), ::tolower);
    auto targetRegex = std::regex(target);

	if (target != "")
	{
		for (auto & t : res.stats.taxonomies)
		{
			auto & tax = t.second;
            if (std::regex_match(tax.identifier, targetRegex))
			{
				tax.state = Clean;
			} else if (tax.identifier != "unknown")
			{
				tax.state = Contamination;
			}
		}
	}

	// create per cluster label mapping and vice versa
	const auto & optClust = *(res.clusterings.mostLikelyClustering);
	
	std::map<std::string, unsigned> clusterPerContig; 
    for (unsigned i = 0; i < optClust.labels.size(); i++)
    {
        unsigned lbl = optClust.labels.at(i);
    	std::string contig = res.seqVectorization.contigs.at(i);
        clusterPerContig[contig] = lbl;
    }

    // find clusters containing clean contigs and, if exist, set others to contamination

    std::set<unsigned> cleanClusters;

    for (auto & t : res.stats.taxonomies)
    {
    	unsigned cluster = clusterPerContig[t.first];
    	auto & tax = t.second;
    	if(tax.state == Clean)
    	{
    		cleanClusters.insert(cluster);
    	}
    }

    if(cleanClusters.size() > 0)
    {
	    for (auto & t : res.stats.taxonomies)
	    {
	    	unsigned cluster = clusterPerContig[t.first];
	    	auto & tax = t.second;
	    	if(cleanClusters.find(cluster) == cleanClusters.end() && tax.state == NA)
	    	{
	    		tax.state = Contamination;
	    	}
	    }
    }

}

void TaxonomyAnnotation::annotateFromFile(ResultContainer & res, const std::string & filename)
{
	initialize(res);
	// line format: contig id<tab>taxonomy

	std::ifstream infile(filename);
	std::string line;

	const auto & contigs = res.stats.includedContigs;

	while(std::getline(infile, line))
	{
		if(line == "")
		{
			continue;
		}

		std::stringstream linestream(line);
		std::string contig;
		std::string taxonomy;
		std::getline(linestream, contig, '\t');
		std::getline(linestream, taxonomy, '\t');
		std::transform(taxonomy.begin(), taxonomy.end(), taxonomy.begin(), ::tolower);

		if(std::find(contigs.begin(), contigs.end(), contig) != contigs.end())
		{
			res.stats.taxonomies[contig].type = Annotated;
			res.stats.taxonomies[contig].identifier = taxonomy;
		}
	}

	updateContaminationConfidences(res);
	updateContaminationState(res);
}

void TaxonomyAnnotation::initialize(ResultContainer & res)
{
	for (const auto & c : res.stats.includedContigs)
	{
		res.stats.taxonomies[c].type = Estimated;
		res.stats.taxonomies[c].identifier = "unknown";
		res.stats.taxonomies[c].confidence = -1.0;
		res.stats.taxonomies[c].state = NA;
	}
}

struct TaxonomyComparison
{
	bool operator()(const Taxonomy* lhs, const Taxonomy* rhs) const
	{
		return lhs->confidence < rhs->confidence;
	}
};

void TaxonomyAnnotation::annotateUnknown(ResultContainer & res)
{
	// 1. for each cluster iterate over contained contigs
	// 2. create priorityqueue with taxonomies of each contig, sorted by confidence
	// 3. assign taxonomy with highest confidence to all unknown contigs

	const auto & optClust = *(res.clusterings.mostLikelyClustering);

	std::map<unsigned, std::priority_queue<Taxonomy*, std::vector<Taxonomy*>, TaxonomyComparison>> sortedTaxByCluster; 
    for (unsigned i = 0; i < optClust.labels.size(); i++)
    {
        unsigned lbl = optClust.labels.at(i);
    	std::string contig = res.seqVectorization.contigs.at(i);
		sortedTaxByCluster[lbl].push(& res.stats.taxonomies.at(contig));    	
    }

    for (auto & it : sortedTaxByCluster)
    {
    	auto & taxonomies = it.second;
    	auto & mostProbableTax = *(taxonomies.top());
    	taxonomies.pop();

    	while(!taxonomies.empty())
    	{
    		auto & tax = *(taxonomies.top());
    		if (tax.identifier == "unknown")
    		{
    			tax.identifier = mostProbableTax.identifier;
    			tax.confidence = mostProbableTax.confidence;
    		}
    		taxonomies.pop();
    	}
    }

    updateContaminationState(res);
}