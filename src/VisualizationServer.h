#pragma once

#include "Clustering.h"
#include "Controller.h"
#include "WebServer.h"
#include "Util.h"
#include "ClusterAnalysis.h"

#include <eigen3/Eigen/Dense>
#include <string>
#include <unordered_map>
#include <map>
#include <mutex>

class DatasetController : public SimpleGetController
{

public:
	DatasetController(const std::string path_);
	virtual ~DatasetController();
    virtual void respond(std::stringstream & response, const std::map<std::string, std::string> params);

};

struct VisualizationData
{
	Eigen::MatrixXd dataPca;
	Eigen::MatrixXd dataSne;
	std::vector<std::string> labels;	
	ClusterAnalysisResult clustRes;
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

	std::unordered_map< std::string, std::unique_ptr<VisualizationData> > datasets;

public:
	static VisualizationServer & getInstance();

	void run(const unsigned port);
	void stop();
	void addClustering(
		const std::string & key,
		const Eigen::MatrixXd & dataPca,
		const Eigen::MatrixXd & dataSne, 
		const std::vector<std::string> & labels, 
		const ClusterAnalysisResult & ClusterAnalysisResult);
	const VisualizationData * getClustering(const std::string & name);
	
};
