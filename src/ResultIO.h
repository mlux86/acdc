#pragma once

#include <vector>
#include <string>
#include <eigen3/Eigen/Dense>
#include <json/json.h>

#include "Clustering.h"
#include "ClusterAnalysis.h"
#include "KrakenAdapter.h"

struct ResultContainer
{
	unsigned id;

	std::string fasta;
	std::vector<std::string> fastaLabels;

	ClusterAnalysisResult oneshot;
	std::vector<ClusterAnalysisResult> bootstraps;

	KrakenResult kraken;	
};

class ResultIO
{

private:    
	std::vector<unsigned> numericLabels(const std::vector<std::string> & labels);	
	void filterFasta(const std::string & fasta, const std::set<std::string> contigs, const std::string & exportFilename);	

    Json::Value clusteringResultToJSON(const ClusteringResult & cr);
    Json::Value clusterAnalysisResultToJSON(const ResultContainer & result, const ClusterAnalysisResult & car, bool alignToOneshot);
	void writeResultContainerToJSON(ResultContainer result, const std::string & filename);
    void exportClusteringFastas(const ResultContainer & result);

	std::string outputDir;
    bool krakenEnabled;
    std::vector<std::string> jsonFiles;

public:
	ResultIO(const std::string & outputDir, bool krakenEnabled);
	~ResultIO();
	
	void processResult(const ResultContainer & result);
	void finish();
};
