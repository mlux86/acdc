#include <fstream>
#include <map>
#include <unordered_map>
#include <set>
#include <boost/filesystem.hpp>
#include <seqan/alignment_free.h> 
#include <seqan/sequence.h> 
#include <seqan/seq_io.h>

#include "MatrixUtil.h"
#include "MLUtil.h"
#include "ResultIO.h"

ResultIO::ResultIO(const std::string & dir) : outputDir(dir)
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

void ResultIO::filterFasta(const std::string & fasta, const std::set<std::string> contigs, const std::string & exportFilename)
{

    std::stringstream ss;
    
    seqan::SeqFileIn seqFileIn(fasta.c_str());
    seqan::StringSet<seqan::CharString> ids;
    seqan::StringSet<seqan::String<seqan::Iupac> > seqs;

    seqan::readRecords(ids, seqs, seqFileIn);

    unsigned n = seqan::length(ids);

    for (unsigned i = 0; i < n; i++)
    {
        std::string id;
        move(id, ids[i]);
        if (std::find(contigs.begin(), contigs.end(), id) != contigs.end())
        {
            std::string seq;
            move(seq, seqs[i]);
            ss << '>' << id << '\n' << seq << '\n';
        }
    }

    std::ofstream ofs(exportFilename, std::ofstream::out);
    ofs << ss.str();
    ofs.close();
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

	root["dataSne"] = MatrixUtil::matrixToJSON(car.dataSne);
	root["dataPca"] = MatrixUtil::matrixToJSON(car.dataPca);

	root["bootstrapIndices"] = Json::Value(Json::arrayValue);
	for (const auto & idx : car.bootstrapIndexes)
	{
		root["bootstrapIndices"].append(idx);
	}

	root["resConnComponents"] = clusteringResultToJSON(car.resConnComponents);
	root["resDipMeans"] = clusteringResultToJSON(car.resDipMeans);

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

	root["oneshot"] = clusterAnalysisResultToJSON(result.oneshot);
	
	root["bootstraps"] = Json::Value(Json::arrayValue);
	for (auto & bs : result.bootstraps)
	{
		// align bootstrap data to oneshot data (data points and labels)
		bs.dataSne = MLUtil::alignBootstrap(result.oneshot.dataSne, bs.dataSne, bs.bootstrapIndexes);        
		bs.resConnComponents.labels = MLUtil::alignBootstrapLabels(result.oneshot.resConnComponents.labels, bs.resConnComponents.labels, bs.bootstrapIndexes);
		bs.resDipMeans.labels = MLUtil::alignBootstrapLabels(result.oneshot.resDipMeans.labels, bs.resDipMeans.labels, bs.bootstrapIndexes);

		root["bootstraps"].append(clusterAnalysisResultToJSON(bs));
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
	std::ofstream ofs(filename, std::ofstream::out | std::ofstream::app);
	ofs << "results['" << result.fasta << "'] = ";
	writer->write(root, &ofs);
	ofs << ";" << std::endl;
	ofs.close();
}

void ResultIO::exportClusteringFastas(const ResultContainer & result)
{
    // create export directory

	boost::filesystem::path exportPath (outputDir + "/export");
	boost::system::error_code returnedError;
	boost::filesystem::create_directories(exportPath, returnedError);
	if (returnedError)
	{
		throw std::runtime_error("Could not create output directory, aborting.");
	}

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
	std::map<unsigned, std::set<std::string>> contigIdsDip; // export a set of contigs per label
	std::map<unsigned, std::set<std::string>> contigIdsKraken; // export a set of contigs per label

    for (unsigned i = 0; i < result.fastaLabels.size(); ++i)
    {        
        contigIdsCC[result.oneshot.resConnComponents.labels.at(i)].insert(result.fastaLabels.at(i));
        contigIdsDip[result.oneshot.resDipMeans.labels.at(i)].insert(result.fastaLabels.at(i));
        contigIdsKraken[numericLabels(krakenLabels).at(i)].insert(result.fastaLabels.at(i));
    }	

    // connected components
    for (auto & kv : contigIdsCC) 
    {
    	unsigned lbl = kv.first;
    	auto & contigs = kv.second;
    	std::stringstream ss;
    	ss << outputDir << "/export/" << result.id << "-cc-" << lbl << ".fasta";
    	filterFasta(result.fasta, contigs, ss.str());
    }

    // dip means
    for (auto & kv : contigIdsDip) 
    {
    	unsigned lbl = kv.first;
    	auto & contigs = kv.second;
    	std::stringstream ss;
    	ss << outputDir << "/export/" << result.id << "-dip-" << lbl << ".fasta";
    	filterFasta(result.fasta, contigs, ss.str());
    }

    // kraken
    for (auto & kv : contigIdsKraken) 
    {
    	unsigned lbl = kv.first;
    	auto & contigs = kv.second;
    	std::stringstream ss;
    	ss << outputDir << "/export/" << result.id << "-kraken-" << lbl << ".fasta";
    	filterFasta(result.fasta, contigs, ss.str());
    }

}

void ResultIO::processResult(const ResultContainer & result)
{
	std::stringstream ss;
	ss << outputDir << "/" << "data.js";
	std::string fname = ss.str();
	exportClusteringFastas(result);
	writeResultContainerToJSON(result, fname);
	jsonFiles.push_back(fname);
}

void ResultIO::finish()
{
}