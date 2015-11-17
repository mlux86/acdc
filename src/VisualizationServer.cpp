#include "Logger.h"
#include "VisualizationServer.h"

#include <thread>
#include <chrono>

DatasetController::DatasetController(const std::string path_) : SimpleGetController(path_)
{
}

DatasetController::~DatasetController()
{
}

void DatasetController::respond(std::stringstream & response, const std::map<std::string, std::string> params)
{
    Json::Value root;

    if (params.find(ParamInfo) != params.end())
    {
        Json::Value availData = Json::Value(Json::arrayValue);
        for (const auto & name : VisualizationServer::getInstance().getClusteringNames())
        {
            availData.append(name);
        }
        root["datasets"] = availData;

        Json::Value availLabels = Json::Value(Json::arrayValue);
        availLabels.append(LabelsConnComp);
        availLabels.append(LabelsDip);
        availLabels.append(LabelsOrig);
        root["labels"] = availLabels;

        Json::Value availReductions = Json::Value(Json::arrayValue);
        availReductions.append(ReductionSne);
        availReductions.append(ReductionPca);
        root["reductions"] = availReductions;

        response << root;
        return;
    }

    if (params.find(ParamData) == params.end() || params.find(ParamLabels) == params.end() || params.find(ParamReduction) == params.end())
    {
        response << root;
        return;
    }

    const std::string key = params.at(ParamData);
    const std::string labels = params.at(ParamLabels);
    const std::string reduction = params.at(ParamReduction);

    try 
    {
        const VisualizationData * vdat = VisualizationServer::getInstance().getClustering(key);

        const Eigen::MatrixXd * shownData;
        if (reduction == ReductionPca)
        {
            shownData = &(vdat->dataPca);
        } else
        {
            shownData = &(vdat->dataSne);
        }

        if (labels == LabelsConnComp)
        {
            root["mat"] = Util::clusteringToJson(*shownData, vdat->clustRes.resConnComponents.labels, vdat->labels);
        } else if (labels == LabelsDip)
        {
            root["mat"] = Util::clusteringToJson(*shownData, vdat->clustRes.resDipMeans.labels, vdat->labels);
        } else if (labels == LabelsOrig)
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
    VisualizationData * res = nullptr;

    while (true)
    {
        dataMtx.lock();
        if (datasets.find(name) != datasets.end())
        {
            res = datasets.at(name).get();            
        }
        dataMtx.unlock();
        if(res == nullptr)
        {
            std::this_thread::sleep_for(std::chrono::seconds(1));        
        } else
        {
            break;
        }
    }

	return res;

    // return datasets.at(name).get();
}

const std::vector<std::string> VisualizationServer::getClusteringNames()
{
    dataMtx.lock();
    std::vector<std::string> res;
    for (const auto & it : datasets)
    {
        res.push_back(it.first);
    }
    dataMtx.unlock();
    std::sort(res.begin(), res.end());
    return res;
}