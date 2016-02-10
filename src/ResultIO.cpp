#include <fstream>

#include "ResultIO.h"

ResultIO::ResultIO(const std::string & dir) : outputDir(dir)
{

}

ResultIO::~ResultIO() 
{

}

Json::Value ResultIO::matrixToJSON(const Eigen::MatrixXd & mat)
{
    unsigned n = mat.rows();
    unsigned dim = mat.cols();

	auto arr = Json::Value(Json::arrayValue);

    for(unsigned i = 0; i < n; i++)
    {
        auto entry = Json::Value(Json::arrayValue);
    	for(unsigned j = 0; j < dim; j++)
    	{
    		entry.append(mat(i, j));
    	}
        arr.append(entry);
    }

	return arr;
}

Json::Value ResultIO::clusteringResultToJSON(const ClusteringResult & cr)
{
	auto root = Json::Value();

	root["numClusters"] = cr.numClusters;

	root["labels"] = Json::Value(Json::arrayValue);
	for (const auto & lbl : cr.labels)
	{
		root["labels"].append(lbl);
	}

	return root;
}

Json::Value ResultIO::clusterAnalysisResultToJSON(const ClusterAnalysisResult & car)
{
	auto root = Json::Value();

	root["dataSne"] = ResultIO::matrixToJSON(car.dataSne);
	root["dataPca"] = ResultIO::matrixToJSON(car.dataPca);

	root["bootstrapIndices"] = Json::Value(Json::arrayValue);
	for (const auto & idx : car.bootstrapIndexes)
	{
		root["bootstrapIndices"].append(idx);
	}

	root["resConnComponents"] = ResultIO::clusteringResultToJSON(car.resConnComponents);
	root["resDipMeans"] = ResultIO::clusteringResultToJSON(car.resDipMeans);

	return root;
}

void ResultIO::writeResultContainerToJSON(ResultContainer result, const std::string & filename)
{
	auto root = Json::Value();

	root["fasta"] = result.fasta;

	root["fastaLabels"] = Json::Value(Json::arrayValue);
	for (auto & lbl : result.fastaLabels)
	{
		root["fastaLabels"].append(lbl);
	}

	root["oneshot"] = ResultIO::clusterAnalysisResultToJSON(result.oneshot);
	
	root["bootstraps"] = Json::Value(Json::arrayValue);
	for (auto & bs : result.bootstraps)
	{
		// align bootstrap data to oneshot data (data points and labels)
		bs.dataSne = Util::alignBootstrap(result.oneshot.dataSne, bs.dataSne, bs.bootstrapIndexes);        
		bs.resConnComponents.labels = Util::alignBootstrapLabels(result.oneshot.resConnComponents.labels, bs.resConnComponents.labels, bs.bootstrapIndexes);
		bs.resDipMeans.labels = Util::alignBootstrapLabels(result.oneshot.resDipMeans.labels, bs.resDipMeans.labels, bs.bootstrapIndexes);

		root["bootstraps"].append(ResultIO::clusterAnalysisResultToJSON(bs));
	}	

	root["krakenLabels"] = Json::Value(Json::arrayValue);
    for (auto & lbl : result.fastaLabels)
    {
        std::string krakenLbl = "unknown";
        if (result.kraken.classification.find(lbl) != result.kraken.classification.end())
        {
            krakenLbl = result.kraken.classification.at(lbl);
        }
        root["krakenLabels"].append(krakenLbl);
    }

	Json::StreamWriterBuilder builder;
	std::unique_ptr<Json::StreamWriter> writer(builder.newStreamWriter());
	std::ofstream ofs(filename, std::ofstream::out);
	ofs << "results['" << result.fasta << "'] = ";
	writer->write(root, &ofs);
	ofs.close();
}

void ResultIO::processResult(const ResultContainer & result)
{
	std::stringstream ss;
	ss << outputDir << "/" << jsonFiles.size() << ".js";
	std::string fname = ss.str();
	ResultIO::writeResultContainerToJSON(result, fname);
	jsonFiles.push_back(fname);
}

void ResultIO::finish()
{
}