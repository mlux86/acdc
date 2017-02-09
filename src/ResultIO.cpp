#include <fstream>
#include <map>
#include <unordered_map>
#include <set>
#include <memory>

#include <seqan/alignment_free.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

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

void ResultIO::exportClusteringInfo(const ResultContainer & result, const std::string & filename)
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

    auto root = Json::Value();    

    // CC
    root["cc"] = Json::Value();
    for (auto & kv : contigIdsCC)
    {
        unsigned lbl = kv.first;
        std::string lblStr = std::to_string(lbl);
        auto & contigs = kv.second;
        root["cc"][lblStr] = Json::Value(Json::arrayValue);
        for (const auto & c : contigs)
        {
            root["cc"][lblStr].append(c);
        }
    }

    // kraken
    root["kraken"] = Json::Value();
    if (krakenEnabled)
    {
        for (auto & kv : contigIdsKraken)
        {
            unsigned lbl = kv.first;
            std::string lblStr = std::to_string(lbl);
            auto & contigs = kv.second;
            root["kraken"][lblStr] = Json::Value(Json::arrayValue);
            for (const auto & c : contigs)
            {
                root["kraken"][lblStr].append(c);
            }
        }
    }

    // dip sne & dip pca

    root["dip"] = Json::Value();
    root["dip"]["dataSne"] = Json::Value();
    root["dip"]["dataPca"] = Json::Value();

    unsigned maxK = 5;
    for (unsigned k = 0; k < maxK; k++)
    {
        std::string kStr = std::to_string(k+1);

        std::map<unsigned, std::set<std::string>> contigIdsSne; // export a set of contigs per label
        std::map<unsigned, std::set<std::string>> contigIdsPca; // export a set of contigs per label

        for (unsigned i = 0; i < result.fastaLabels.size(); ++i)
        {
            contigIdsSne[result.oneshot.clustsSne.at(k).labels.at(i)].insert(result.fastaLabels.at(i));
            contigIdsPca[result.oneshot.clustsPca.at(k).labels.at(i)].insert(result.fastaLabels.at(i));
        }

        root["dip"]["dataSne"][kStr] = Json::Value();
        for (auto & kv : contigIdsSne)
        {
            unsigned lbl = kv.first;
            std::string lblStr = std::to_string(lbl);
            auto & contigs = kv.second;
            root["dip"]["dataSne"][kStr][lblStr] = Json::Value(Json::arrayValue);
            for (const auto & c : contigs)
            {
                root["dip"]["dataSne"][kStr][lblStr].append(c);
            }
        }

        root["dip"]["dataPca"][kStr] = Json::Value();
        for (auto & kv : contigIdsPca)
        {
            unsigned lbl = kv.first;
            std::string lblStr = std::to_string(lbl);
            auto & contigs = kv.second;
            root["dip"]["dataPca"][kStr][lblStr] = Json::Value(Json::arrayValue);
            for (const auto & c : contigs)
            {
                root["dip"]["dataPca"][kStr][lblStr].append(c);
            }
        }
    }

    Json::StreamWriterBuilder builder;
    std::unique_ptr<Json::StreamWriter> writer(builder.newStreamWriter());
    std::ofstream ofs(filename, std::ofstream::out | std::ofstream::app);
    ofs << "clusterinfo['" << result.fasta << "'] = ";
    writer->write(root, &ofs);
    ofs << ";" << std::endl;
    ofs.close();

}

void ResultIO::exportContigJS(const std::string & fastaFilename)
{
    seqan::SeqFileIn seqFileIn(fastaFilename.c_str());
    seqan::StringSet<seqan::CharString> ids;
    seqan::StringSet<seqan::String<seqan::Iupac> > seqs;
    seqan::readRecords(ids, seqs, seqFileIn);
    unsigned n = seqan::length(ids);

    std::stringstream ss;
    ss << outputDir << "/contigs.js";
    std::ofstream ofs(ss.str(), std::ofstream::out | std::ofstream::app);

    for (unsigned i = 0; i < n; i++)
    {
        std::string id;
        move(id, ids[i]);
        std::string seq;
        move(seq, seqs[i]);
        ofs << "inputcontigs['" << id << "'] = '" << seq << "'" << std::endl;
    }

    ofs.close();    
}

void ResultIO::processResult(const ResultContainer & result)
{
    export16S(result);

    exportContigJS(result.fasta);

    std::stringstream ss;
    ss << outputDir << "/" << "clusterinfo.js";
    std::string fname = ss.str();    
    exportClusteringInfo(result, fname);

	std::stringstream ss2;
	ss2 << outputDir << "/" << "data.js";
	std::string fname2 = ss2.str();
	writeResultContainerToJSON(result, fname2);

}
