#pragma once

#include "Clustering.h"
#include "Controller.h"
#include "WebServer.h"
#include "Util.h"

#include <eigen3/Eigen/Dense>
#include <string>
#include <unordered_map>
#include <map>
#include <fstream>
#include <streambuf>
#include <mutex>

class VisualizationServer
{

private:
	VisualizationServer();
	~VisualizationServer();
	VisualizationServer(VisualizationServer const &) = delete;
    void operator=(VisualizationServer const &) = delete;

    std::unique_ptr<WebServer> server;
    std::mutex dataMtx;

	// std::unordered_map< std::string, std::pair<const Eigen::MatrixXd *, const ClusteringResult*> > datasets;
	std::unordered_map< std::string, std::pair< std::unique_ptr<Eigen::MatrixXd>, std::unique_ptr<ClusteringResult> > > datasets;

public:
	static VisualizationServer & getInstance();

	void run(const unsigned port);
	void stop();
	void addClustering(const Eigen::MatrixXd & data, const ClusteringResult & clust, const std::string & name);
	std::pair<const Eigen::MatrixXd *, const ClusteringResult*> getClustering(const std::string & name);
	
};

class DatasetController : public SimpleGetController
{

private:
	VisualizationServer * visServer;

public:
	DatasetController(const std::string path_, VisualizationServer * vs) : SimpleGetController(path_), visServer(vs)
	{

	}

    virtual void respond(std::stringstream & response, const std::map<std::string, std::string> params)
    {
    	Json::Value root;

		if (params.find("data") == params.end())
		{
			response << root;
			return;
		}

		const std::string key = params.at("data");

    	try 
    	{
	    	auto dat = visServer->getClustering(key);
	    	root["mat"] = Util::clusteringToJson(*(dat.first), *(dat.second));    	
	    	response << root; 
    	}
    	catch(const std::out_of_range & oor) 
    	{
    		Json::Value root;
    		root["error"] = "No data with key '" + key + "'";
    		response << root;
    	}    	
    	catch(const std::exception & e) 
    	{
    		ELOG << e.what() << "\n";
    	}
    }

};

class StaticController : public SimpleGetController
{

private:
	std::string filename;

public:
	StaticController(const std::string path_, const std::string f) : SimpleGetController(path_), filename(f)
	{

	}

    virtual void respond(std::stringstream & response, const std::map<std::string, std::string> params)
    {
		std::ifstream t(filename);
		std::string str((std::istreambuf_iterator<char>(t)),
		                 std::istreambuf_iterator<char>());    	
		response << str;
    }

};

