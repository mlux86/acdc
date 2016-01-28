#pragma once

#include "Clustering.h"
#include "Controller.h"
#include "WebServer.h"
#include "Util.h"
#include "ClusterAnalysis.h"
#include "KrakenAdapter.h"

#include <eigen3/Eigen/Dense>
#include <string>
#include <unordered_map>
#include <map>
#include <mutex>

class DatasetController : public SimpleGetController
{
private:
	const std::string ParamInfo = "info";
	const std::string ParamData = "data";
	const std::string ParamLabels = "labels";
	const std::string ParamReduction = "reduction";
	const std::string LabelsOrig = "orig";
	const std::string LabelsDip = "dip";
	const std::string LabelsConnComp = "cc";
	const std::string LabelsKraken = "kraken";
	const std::string ReductionPca = "pca";
	const std::string ReductionSne = "sne";
	const std::string ParamOneshot = "oneshot";
	const std::string ParamBootstrapId = "bootstrapId";

public:
	DatasetController(const std::string path_);
	virtual ~DatasetController();
    virtual void respond(std::stringstream & response, const std::map<std::string, std::string> params);

};

class StatsController : public SimpleGetController
{
public:
	StatsController(const std::string path_);
	virtual ~StatsController();
	virtual void respond(std::stringstream & response, const std::map<std::string, std::string> params);	

	std::map<unsigned, double> getClusterConfidences(const std::string & dataset, const std::string & clust);
};

struct VisualizationDataEntry
{
	Eigen::MatrixXd dataPca;
	Eigen::MatrixXd dataSne;
	std::vector<std::string> labels;	
	ClusterAnalysisResult clustRes;
};

struct VisualizationData
{
	VisualizationDataEntry oneshot;
	std::vector<VisualizationDataEntry> bootstraps;
	std::vector<std::string> krakenClassification;
};

class VisualizationServer
{

private:
	VisualizationServer();
	~VisualizationServer();
	VisualizationServer(VisualizationServer const &) = delete;
    void operator=(VisualizationServer const &) = delete;

    std::unique_ptr<WebServer> server;
    std::mutex dataMtx;

    std::unordered_map< std::string, std::unique_ptr<VisualizationData>> datasets;

public:
	static VisualizationServer & getInstance();

	void run(const unsigned port);
	void stop();
	void addClustering(
		const std::string & key,
		bool oneshot,
		const Eigen::MatrixXd & dataPca,
		const Eigen::MatrixXd & dataSne, 
		const std::vector<std::string> & labels, 
		const ClusterAnalysisResult & ClusterAnalysisResult);
	const VisualizationData * getClustering(const std::string & name);
	const std::vector<std::string> getClusteringNames();
	void addKrakenResult(const std::string & name, const KrakenResult & krakenResult);
};
