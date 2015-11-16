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
    DatasetController c("/json", this);
    StaticController s("/", "../scatter.html");
    StaticController cjs("/canvasjs", "../canvasjs.min.js");

    server.reset(new WebServer(port));
    server->addController(& c);
    server->addController(& s);
    server->addController(& cjs);
    server->start();
}

void VisualizationServer::stop()
{
	server->stop();
}

void VisualizationServer::addClustering(const Eigen::MatrixXd & data_, const ClusteringResult & clust_, const std::string & name)
{
    dataMtx.lock();

    Eigen::MatrixXd * data = new Eigen::MatrixXd;
    *data = data_;
    ClusteringResult * clust = new ClusteringResult;
    *clust = clust_;
    datasets[name] = std::make_pair(std::unique_ptr<Eigen::MatrixXd>(data), std::unique_ptr<ClusteringResult>(clust));

    DLOG << "Added data set to visualization server: " << name << "\n";
    
    dataMtx.unlock();
}

std::pair<const Eigen::MatrixXd *, const ClusteringResult*> VisualizationServer::getClustering(const std::string & name) 
{
    auto & v = datasets.at(name);
	return std::make_pair(v.first.get(), v.second.get());
}