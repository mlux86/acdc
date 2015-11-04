#pragma once

#include "Clustering.h"
#include "Controller.h"
#include "WebServer.h"

#include <eigen3/Eigen/Dense>
#include <string>
#include <unordered_map>

class HelloController : public SimpleGetController
{

public:
	HelloController(const std::string path_) : SimpleGetController(path_) 
	{

	}

    virtual void respond(std::stringstream & response)
    {
    	response << "Hello.\n"; 
    }

};

class VisualizationServer
{

private:
	VisualizationServer();
	~VisualizationServer();
	VisualizationServer(VisualizationServer const &) = delete;
    void operator=(VisualizationServer const &) = delete;

    WebServer * server;

	std::unordered_map< std::string, std::pair<const Eigen::MatrixXd*, const ClusteringResult*> > datasets;

public:
	static VisualizationServer & getInstance();

	void run(const unsigned port);
	void stop();
	void serveDataset();
	void addClustering(const Eigen::MatrixXd & data, const ClusteringResult & clust, const std::string & name);
	
};