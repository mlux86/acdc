#include <fstream>
#include <map>
#include <unordered_map>
#include <set>
#include <memory>

#include "MatrixUtil.h"
#include "MLUtil.h"
#include "ResultIO.h"
#include "SequenceUtil.h"
#include "Clustering.h"

ResultIO::ResultIO(const std::string & dir, bool krakenEnabled_) : outputDir(dir), krakenEnabled(krakenEnabled_)
{

}

ResultIO::~ResultIO()
{

}

std::vector<unsigned> ResultIO::numericLabels(const std::vector<std::string> & labels)
{
    unsigned n = labels.size();
    std::vector<unsigned> v(n, 0);
    std::unordered_map<std::string, unsigned> lblMap;
    unsigned cnt = 0;
    unsigned i = 0;
    for (const auto & lbl : labels)
    {
        if(lblMap.find(lbl) == lblMap.end())
        {
            lblMap[lbl] = cnt++;
        }
        v[i++] = lblMap[lbl];
    }
    return v;
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

    root["outlierClusters"] = Json::Value(Json::arrayValue);
    for (const auto & lbl : cr.outlierClusters)
    {
        root["outlierClusters"].append(lbl);
    }

	return root;
}

Json::Value ResultIO::clusterAnalysisResultToJSON(const ResultContainer & result, const ClusterAnalysisResult & car, bool alignToOneshot)
{
	auto root = Json::Value();

    if (alignToOneshot)
    {
       root["dataSne"] = MatrixUtil::matrixToJSON(MLUtil::alignBootstrap(result.oneshot.dataSne, car.dataSne, car.bootstrapIndexes));
    } else
    {
	   root["dataSne"] = MatrixUtil::matrixToJSON(car.dataSne);
    }
	root["dataPca"] = MatrixUtil::matrixToJSON(car.dataPca);

	root["bootstrapIndexes"] = Json::Value(Json::arrayValue);
	for (const auto & idx : car.bootstrapIndexes)
	{
		root["bootstrapIndexes"].append(idx);
	}

    root["isMultiModal"] = car.isMultiModal;

    root["numClustPca"] = car.numClustPca;
	root["clustsPca"] = Json::Value(Json::arrayValue);
    for (auto & c : car.clustsPca)
    {
        root["clustsPca"].append(clusteringResultToJSON(c));
    }

    root["numClustSne"] = car.numClustSne;
    root["clustsSne"] = Json::Value(Json::arrayValue);
    for (auto & c : car.clustsSne)
    {
        root["clustsSne"].append(clusteringResultToJSON(c));
    }

    root["hasSeparatedComponents"] = car.hasSeparatedComponents;
    root["numClustCC"] = car.numClustCC;
    root["clustCC"] = clusteringResultToJSON(car.clustCC);

	return root;
}

void ResultIO::writeResultContainerToJSON(ResultContainer result, const std::string & filename)
{
	auto root = Json::Value();

	root["fasta"] = result.fasta;
	root["id"] = result.id;

	root["fastaLabels"] = Json::Value(Json::arrayValue);
	for (auto & lbl : result.fastaLabels)
	{
		root["fastaLabels"].append(lbl);
	}

	root["oneshot"] = clusterAnalysisResultToJSON(result, result.oneshot, false);

	root["bootstraps"] = Json::Value(Json::arrayValue);
	for (auto & bs : result.bootstraps)
	{
		root["bootstraps"].append(clusterAnalysisResultToJSON(result, bs, true));
	}

    if (krakenEnabled)
    {
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
        root["krakenBacterialBackground"] = result.kraken.bacterialBackground;
    }

    root["contains16S"] = Json::Value(Json::arrayValue);
    for (const std::string & s : result._16S)
    {
        root["contains16S"].append(!s.empty());
    }

    root["stats"] = Json::Value();
    root["stats"]["numBasepairs"] = (unsigned)result.stats.numBasepairs;
    root["stats"]["gcContent"] = result.stats.gcContent;
    root["stats"]["contigGcContent"] = Json::Value();
    for (const auto & entry : result.stats.contigGcContent)
    {
        root["stats"]["contigGcContent"][entry.first] = entry.second;
    }
    root["stats"]["contigLength"] = Json::Value();
    for (const auto & entry : result.stats.contigLength)
    {
        root["stats"]["contigLength"][entry.first] = (unsigned)entry.second;
    }

	Json::StreamWriterBuilder builder;
	std::unique_ptr<Json::StreamWriter> writer(builder.newStreamWriter());
	std::ofstream ofs(filename, std::ofstream::out | std::ofstream::app);
	ofs << "results['" << result.fasta << "'] = ";
	writer->write(root, &ofs);
	ofs << ";" << std::endl;
	ofs.close();
}

void ResultIO::export16S(const ResultContainer & result)
{
    for (unsigned i = 0; i < result._16S.size(); i++)
    {
        const std::string & s = result._16S.at(i);
        if (!s.empty())
        {
            std::stringstream ss;
            ss << outputDir << "/export/" << result.id << "-" << i << ".16s";
            std::ofstream ofs(ss.str(), std::ofstream::out);
            ofs << s << std::endl;
            ofs.close();
        }
    }
}

void ResultIO::exportClusteringFastas(const ResultContainer & result)
{

	// first, we need the kraken result as actual labels
	std::vector<std::string> krakenLabels;
    for (auto & lbl : result.fastaLabels)
    {
        std::string krakenLbl = "unknown";
        if (result.kraken.classification.find(lbl) != result.kraken.classification.end())
        {
            krakenLbl = result.kraken.classification.at(lbl);
        }
        krakenLabels.push_back(krakenLbl);
    }

    // start export

	std::map<unsigned, std::set<std::string>> contigIdsCC; // export a set of contigs per label
	std::map<unsigned, std::set<std::string>> contigIdsKraken; // export a set of contigs per label

    for (unsigned i = 0; i < result.fastaLabels.size(); ++i)
    {
        contigIdsCC[result.oneshot.clustCC.labels.at(i)].insert(result.fastaLabels.at(i));
        contigIdsKraken[numericLabels(krakenLabels).at(i)].insert(result.fastaLabels.at(i));
    }

    // CC
    for (auto & kv : contigIdsCC)
    {
        unsigned lbl = kv.first;
        auto & contigs = kv.second;
        std::stringstream ss;
        ss << outputDir << "/export/" << result.id << "-cc-" << lbl << ".fasta";
        SequenceUtil::exportFilteredFasta(result.fasta, contigs, ss.str());
    }

    // kraken
    if (krakenEnabled)
    {
        for (auto & kv : contigIdsKraken)
        {
            unsigned lbl = kv.first;
            auto & contigs = kv.second;
            std::stringstream ss;
            ss << outputDir << "/export/" << result.id << "-kraken-" << lbl << ".fasta";
            SequenceUtil::exportFilteredFasta(result.fasta, contigs, ss.str());
        }

        std::stringstream ss;
        ss << outputDir << "/export/" << result.id << ".oneshot.kraken";
        std::ofstream ofs(ss.str(), std::ofstream::out);
        for (auto & lbl : krakenLabels)
        {
            ofs << lbl << std::endl;
        }
        ofs.close();
    }

    // dip sne & dip pca

    unsigned maxK = 5;
    for (unsigned k = 0; k < maxK; k++)
    {
        std::map<unsigned, std::set<std::string>> contigIdsSne; // export a set of contigs per label
        std::map<unsigned, std::set<std::string>> contigIdsPca; // export a set of contigs per label

        for (unsigned i = 0; i < result.fastaLabels.size(); ++i)
        {
            contigIdsSne[result.oneshot.clustsSne.at(k).labels.at(i)].insert(result.fastaLabels.at(i));
            contigIdsPca[result.oneshot.clustsPca.at(k).labels.at(i)].insert(result.fastaLabels.at(i));
        }

        for (auto & kv : contigIdsSne)
        {
        	unsigned lbl = kv.first;
        	auto & contigs = kv.second;
        	std::stringstream ss;
        	ss << outputDir << "/export/" << result.id << "-dip-dataSne-" << (k+1) << "-" << lbl << ".fasta";
        	SequenceUtil::exportFilteredFasta(result.fasta, contigs, ss.str());
        }

        for (auto & kv : contigIdsPca)
        {
            unsigned lbl = kv.first;
            auto & contigs = kv.second;
            std::stringstream ss;
            ss << outputDir << "/export/" << result.id << "-dip-dataPca-" << (k+1) << "-" << lbl << ".fasta";
            SequenceUtil::exportFilteredFasta(result.fasta, contigs, ss.str());
        }
    }


}

void ResultIO::processResult(const ResultContainer & result)
{
	std::stringstream ss;
	ss << outputDir << "/" << "data.js";
	std::string fname = ss.str();
	exportClusteringFastas(result);
    export16S(result);
	writeResultContainerToJSON(result, fname);

    // ss.str("");
    // ss << outputDir << "/export/" << result.id<< ".oneshot.orig";
    // MatrixUtil::saveMatrix(result.oneshot.dataOrig, ss.str(), '\t');
    // ss.str("");
    // ss << outputDir << "/export/" << result.id<< ".oneshot.sne";
    // MatrixUtil::saveMatrix(result.oneshot.dataSne, ss.str(), '\t');
    // ss.str("");
    // ss << outputDir << "/export/" << result.id <<  ".oneshot.pca";
    // MatrixUtil::saveMatrix(result.oneshot.dataPca, ss.str(), '\t');
    // ss.str("");
    // ss << outputDir << "/export/" << result.id << ".oneshot.contigs";
    // std::ofstream ofs(ss.str(), std::ofstream::out);
    // for (auto & lbl : result.fastaLabels)
    // {
    //     ofs << lbl << std::endl;
    // }
    // ofs.close();
}
