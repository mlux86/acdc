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
	std::string fasta;
	std::vector<std::string> fastaLabels;

	ClusterAnalysisResult oneshot;
	std::vector<ClusterAnalysisResult> bootstraps;

	KrakenResult kraken;	
};

class ResultIO
{

private:
    static Json::Value matrixToJSON(const Eigen::MatrixXd & mat);
    static Json::Value clusteringResultToJSON(const ClusteringResult & cr);
    static Json::Value clusterAnalysisResultToJSON(const ClusterAnalysisResult & car);
	static void writeResultContainerToJSON(ResultContainer result, const std::string & filename);

	std::string outputDir;
    std::vector<std::string> jsonFiles;

public:
	ResultIO(const std::string & outputDir);
	~ResultIO();
	
	void processResult(const ResultContainer & result);
	void finish();
};
