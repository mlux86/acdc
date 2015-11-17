#include "Logger.h"
#include "VisualizationServer.h"

DatasetController::DatasetController(const std::string path_) : SimpleGetController(path_)
{
}

DatasetController::~DatasetController()
{
}

void DatasetController::respond(std::stringstream & response, const std::map<std::string, std::string> params)
{
    Json::Value root;

    if (params.find("data") == params.end() || params.find("labels") == params.end() || params.find("reduction") == params.end())
    {
        response << root;
        return;
    }

    const std::string key = params.at("data");
    const std::string labels = params.at("labels");
    const std::string reduction = params.at("reduction");

    try 
    {
        const VisualizationData * vdat = VisualizationServer::getInstance().getClustering(key);

        const Eigen::MatrixXd * shownData;
        if (reduction == "pca")
        {
            shownData = &(vdat->dataPca);
        } else
        {
            shownData = &(vdat->dataSne);
        }

        if (labels == "cc")
        {
            root["mat"] = Util::clusteringToJson(*shownData, vdat->clustRes.resConnComponents.labels, vdat->labels);
        } else if (labels == "dip")
        {
            root["mat"] = Util::clusteringToJson(*shownData, vdat->clustRes.resDipMeans.labels, vdat->labels);
        } else if (labels == "orig")
        {
            root["mat"] = Util::clusteringToJson(*shownData, Util::numericLabels(vdat->labels), vdat->labels);
        }
        
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
    DatasetController c("/json");
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

void VisualizationServer::addClustering(
        const std::string & key,
        const Eigen::MatrixXd & dataPca,
        const Eigen::MatrixXd & dataSne, 
        const std::vector<std::string> & labels,
        const ClusterAnalysisResult & clustRes)
{
    dataMtx.lock();

    std::unique_ptr<VisualizationData> dat(new VisualizationData);
    dat->dataPca = dataPca;
    dat->dataSne = dataSne;
    dat->labels = labels;
    dat->clustRes = clustRes;
    datasets[key] = std::move(dat);    

    DLOG << "Added data set to visualization server: " << key << "\n";
    
    dataMtx.unlock();
}

const VisualizationData * VisualizationServer::getClustering(const std::string & name) 
{
	return datasets.at(name).get();
}