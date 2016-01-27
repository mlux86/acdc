#include "Logger.h"
#include "VisualizationServer.h"

#include <thread>
#include <chrono>

#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/path.hpp"
#include "boost/lexical_cast.hpp"

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
    bool oneshot = (params.find(ParamOneshot) != params.end());

    try 
    {

        const VisualizationData * vdat = VisualizationServer::getInstance().getClustering(key);
        const VisualizationDataEntry * vde;
        if (oneshot)
        {
            vde = &(vdat->oneshot);
        } else
        {
            if (params.find(ParamBootstrapId) == params.end())
            {
                response << root;
                return;
            }            
            unsigned bootstrapId = boost::lexical_cast<unsigned>(params.at(ParamBootstrapId));
            bootstrapId = std::min(bootstrapId, (unsigned)(vdat->bootstraps.size()-1));
            vde = &(vdat->bootstraps[bootstrapId]);
        }

        Eigen::MatrixXd shownData;
        if (reduction == ReductionPca)
        {
            shownData = vde->dataPca;
        } else
        {
            shownData = Util::alignDataset(vdat->oneshot.dataSne, vde->dataSne, vdat->oneshot.labels, vde->labels);
        }
        if (labels == LabelsConnComp)
        {
            root["mat"] = Util::clusteringToJson(shownData, vde->clustRes.resConnComponents.labels, vde->labels);
        } else if (labels == LabelsDip)
        {
            root["mat"] = Util::clusteringToJson(shownData, vde->clustRes.resDipMeans.labels, vde->labels);
        } else if (labels == LabelsOrig)
        {
            root["mat"] = Util::clusteringToJson(shownData, Util::numericLabels(vde->labels), vde->labels);
        }

        root["numBootstraps"] = (unsigned) vdat->bootstraps.size();

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
    std::vector<std::unique_ptr<Controller>> controllers;

    controllers.emplace_back(new DatasetController("/json"));
    controllers.emplace_back(new StatsController("/stats"));

    server.reset(new WebServer(port));

    if (boost::filesystem::is_directory("../assets"))
    {
        boost::filesystem::directory_iterator end_iter;
        for (boost::filesystem::directory_iterator dir_itr("../assets"); dir_itr != end_iter; ++dir_itr)
        {
            if (boost::filesystem::is_regular_file(dir_itr->status()))
            {
                std::string fname = dir_itr->path().filename().string();
                controllers.emplace_back(new StaticController("/" + fname, "../assets/" + fname));
                if (fname == "index.html")
                {
                    controllers.emplace_back(new StaticController("/", "../assets/" + fname));
                }
            }
        }
    }

    for (const auto & c : controllers)
    {
        server->addController(c.get());
    }

    server->start();
}

void VisualizationServer::stop()
{
	server->stop();
}

void VisualizationServer::addClustering(
        const std::string & key,
        bool oneshot,
        const Eigen::MatrixXd & dataPca,
        const Eigen::MatrixXd & dataSne, 
        const std::vector<std::string> & labels,
        const ClusterAnalysisResult & clustRes)
{
    dataMtx.lock();

    if (datasets.find(key) == datasets.end())
    {
        std::unique_ptr<VisualizationData> dat(new VisualizationData);
        datasets[key] = std::move(dat);
    }

    VisualizationData * dat = datasets[key].get();
    
    VisualizationDataEntry vde;
    vde.dataPca = dataPca;
    vde.dataSne = dataSne;
    vde.labels = labels;
    vde.clustRes = clustRes;

    if (oneshot)
    {
        dat->oneshot = vde;
    } else
    {
        dat->bootstraps.push_back(vde);
    }

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

StatsController::StatsController(const std::string path_) : SimpleGetController(path_)
{
}

StatsController::~StatsController()
{
}

void StatsController::respond(std::stringstream & response, const std::map<std::string, std::string> params)
{
    Json::Value root;

    const std::vector<std::string> clusts = {"cc", "dip"};

    for (const auto & dataset : VisualizationServer::getInstance().getClusteringNames())
    {
        
        Json::Value entry;

        for (const auto & clust : clusts)
        {

            const std::map<unsigned, double> confs = getClusterConfidences(dataset, clust);
            Json::Value confsJson = Json::Value(Json::arrayValue);

            for (const auto & it : confs)
            {
                unsigned numClusters = it.first;
                double confidence = it.second;

                Json::Value v;
                v["numClusters"] = numClusters;
                v["confidence"] = confidence;

                confsJson.append(v);
            }

            entry[clust] = confsJson;
        }


        root[dataset] = entry;
    }

    response << root;

}

std::map<unsigned, double> StatsController::getClusterConfidences(const std::string & dataset, const std::string & clust)
{
    if (clust != "cc" && clust != "dip")
    {
        throw std::runtime_error("Unknown clustering!");
    }

    std::map<unsigned, double> confs;

    const VisualizationData & vd = *(VisualizationServer::getInstance().getClustering(dataset));

    unsigned n = 0;

    for (const auto & bs : vd.bootstraps)
    {
        unsigned numClusters;

        if (clust == "dip")
        {
            numClusters = bs.clustRes.resDipMeans.numClusters;
        } else 
        {
            numClusters = bs.clustRes.resConnComponents.numClusters;
        }

        if (confs.find(numClusters) == confs.end())
        {
            confs[numClusters] = 0;
        }
        confs[numClusters]++;
        n++;
    }

    for (auto & it : confs)
    {
        it.second /= n;
    }

    return confs;
}