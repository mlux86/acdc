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
#include "Opts.h"
#include "IOUtil.h"

#include <yaml-cpp/yaml.h>

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

Json::Value ResultIO::contaminationDetectionResultToJSON(const ResultContainer & result)
{
	auto root = Json::Value();

    root["dataSne"] = MatrixUtil::matrixToJSON(result.dataSne);
	root["dataPca"] = MatrixUtil::matrixToJSON(result.dataPca);

    root["confidenceCC"] = result.contaminationAnalysis.confidenceCC;
    root["confidenceDip"] = result.contaminationAnalysis.confidenceDip;
    root["contaminationStatus"] = result.contaminationAnalysis.state;

    root["numClustPca"] = result.clusterings.numClustPca;
	root["clustsPca"] = Json::Value(Json::arrayValue);
    for (auto & c : result.clusterings.clustsPca)
    {
        root["clustsPca"].append(clusteringResultToJSON(c));
    }

    root["numClustSne"] = result.clusterings.numClustSne;
    root["clustsSne"] = Json::Value(Json::arrayValue);
    for (auto & c : result.clusterings.clustsSne)
    {
        root["clustsSne"].append(clusteringResultToJSON(c));
    }

    root["numClustCC"] = result.clusterings.numClustCC;
    root["clustCC"] = clusteringResultToJSON(result.clusterings.clustCC);

	return root;
}

void ResultIO::writeResultContainerToJSON(ResultContainer result, const std::string & filename)
{
	auto root = Json::Value();

	root["fasta"] = result.fasta;
	root["id"] = result.id;

	root["fastaLabels"] = Json::Value(Json::arrayValue);
	for (auto & lbl : result.seqVectorization.contigs)
	{
		root["fastaLabels"].append(lbl);
	}

	root["oneshot"] = contaminationDetectionResultToJSON(result);

    if (krakenEnabled)
    {
    	root["krakenLabels"] = Json::Value(Json::arrayValue);
        for (auto & lbl : result.seqVectorization.contigs)
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
    for (auto & lbl : result.seqVectorization.contigs)
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


    for (unsigned i = 0; i < result.seqVectorization.contigs.size(); ++i)
    {
        contigIdsCC[result.clusterings.clustCC.labels.at(i)].insert(result.seqVectorization.contigs.at(i));
        contigIdsKraken[numericLabels(krakenLabels).at(i)].insert(result.seqVectorization.contigs.at(i));
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

        for (unsigned i = 0; i < result.seqVectorization.contigs.size(); ++i)
        {
            contigIdsSne[result.clusterings.clustsSne.at(k).labels.at(i)].insert(result.seqVectorization.contigs.at(i));
            contigIdsPca[result.clusterings.clustsPca.at(k).labels.at(i)].insert(result.seqVectorization.contigs.at(i));
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

std::vector<unsigned> ResultIO::contigsIndicesWith16S(const ResultContainer & result)
{
    std::vector<std::string> contigs16s;
    for (unsigned i = 0; i < result._16S.size(); i++)
    {
        if(result._16S.at(i).size() > 0)
        {
            contigs16s.push_back(result.seqVectorization.contigs.at(i));
        }
    }

    std::vector<unsigned> contigIndices;

    for (unsigned i = 0; i < result.stats.contigs.size(); i++)
    {
        if(std::find(contigs16s.begin(), contigs16s.end(), result.stats.contigs.at(i)) != contigs16s.end())
        {
            contigIndices.push_back(i);
        }
    }

    return contigIndices;
}

void ResultIO::writeYAML(const ResultContainer & result, const std::string & filename)
{
    std::vector<Fixed> sneCol0 = IOUtil::columnToFixed(result.dataSne, 0);
    std::vector<Fixed> sneCol1 = IOUtil::columnToFixed(result.dataSne, 1);
    std::vector<Fixed> pcaCol0 = IOUtil::columnToFixed(result.dataPca, 0);
    std::vector<Fixed> pcaCol1 = IOUtil::columnToFixed(result.dataPca, 1);

    YAML::Emitter out;
    out << YAML::BeginMap;
    out << YAML::Key << "cli_call" << YAML::Value << Opts::cliCall();
    out << YAML::Key << "acdc_parameters" << YAML::Value << Opts::parameters();
    out << YAML::Key << "contigs" << YAML::Value << YAML::Flow << result.stats.contigs;
    out << YAML::Key << "visualizations" << YAML::Value 
        << YAML::BeginMap 
            << YAML::Key << "sne" << YAML::Value 
                << YAML::BeginMap 
                    << YAML::Key << "x1" << YAML::Value << YAML::Flow << sneCol0 
                    << YAML::Key << "x2" << YAML::Value << YAML::Flow << sneCol1
                << YAML::EndMap
            << YAML::Key << "pca" << YAML::Value 
                << YAML::BeginMap 
                    << YAML::Key << "x1" << YAML::Value << YAML::Flow << pcaCol0 
                    << YAML::Key << "x2" << YAML::Value << YAML::Flow << pcaCol1
                << YAML::EndMap
        << YAML::EndMap;
    out << YAML::Key << "contigs_16s" << YAML::Value << YAML::Flow << contigsIndicesWith16S(result);
    out << YAML::Key << "confidence_cc" << YAML::Value << result.contaminationAnalysis.confidenceCC;
    out << YAML::Key << "confidence_dip" << YAML::Value << result.contaminationAnalysis.confidenceDip;
    out << YAML::Key << "contamination_state" << YAML::Value << result.contaminationAnalysis.state;
    out << YAML::EndMap;

    std::ofstream ofs(filename, std::ofstream::out);
    ofs << out.c_str();
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

    std::stringstream ss3;
    ss3 << outputDir << "/" << "result.yaml";
    std::string fname3 = ss3.str();
    writeYAML(result, fname3);

}
