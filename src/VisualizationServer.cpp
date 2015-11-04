#include "Logger.h"
#include "VisualizationServer.h"

VisualizationServer::VisualizationServer()
{
}

VisualizationServer::~VisualizationServer()
{
}

VisualizationServer & VisualizationServer::getInstance()
{
	static VisualizationServer instance;
	return instance;
}

void VisualizationServer::run(unsigned port)
{
	std::string path = "/";
    HelloController c(path);

    server = new WebServer(port);
    server->addController(& c);
    server->start();
}

void VisualizationServer::stop()
{
	server->stop();
}

void VisualizationServer::serveDataset()
{
}

void VisualizationServer::addClustering(const Eigen::MatrixXd & data, const ClusteringResult & clust, const std::string & name)
{
	datasets[name] = std::make_pair(&data, &clust);
}
