#include <fstream>
#include <map>
#include <unordered_map>
#include <set>
#include <memory>

#include <seqan/alignment_free.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

#include <boost/algorithm/string.hpp>

#include "MatrixUtil.h"
#include "MLUtil.h"
#include "ResultIO.h"
#include "SequenceUtil.h"
#include "Clustering.h"
#include "Opts.h"
#include "IOUtil.h"

ResultIO::ResultIO(const std::string & dir, bool krakenEnabled_) : outputDir(dir), krakenEnabled(krakenEnabled_)
{

}

ResultIO::~ResultIO()
{

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

    ofs << "inputcontigs['" << fastaFilename << "'] = new Array();" << std::endl;
    for (unsigned i = 0; i < n; i++)
    {
        std::string id;
        move(id, ids[i]);
        std::string seq;
        move(seq, seqs[i]);
        ofs << "inputcontigs['" << fastaFilename << "']['" << id << "'] = '" << seq << "'" << std::endl;
    }
    ofs.close();
}

std::vector<unsigned> ResultIO::contigIndicesVisualization(const ResultContainer & result)
{

    std::map<std::string, unsigned> idxMap;
    for (unsigned i = 0; i < result.stats.includedContigs.size(); i++)
    {
        idxMap[result.stats.includedContigs[i]] = i;
    }

    std::vector<unsigned> contigIndices;

    for (auto contig : result.seqVectorization.contigs)
    {
        contigIndices.push_back(idxMap[contig]);
    }

    return contigIndices;
}

std::vector<unsigned> ResultIO::labelsPerContig(const ResultContainer & result, const ClusteringResult & clust)
{

    std::map<std::string, unsigned> contigLabels; // contig -> cluster label

    for (unsigned i = 0; i < clust.labels.size(); i++)
    {
        unsigned lbl = clust.labels.at(i);
        std::string contig = result.seqVectorization.contigs.at(i);
        if (contigLabels.find(contig) == contigLabels.end())
        {
            contigLabels[contig] = lbl;
        }
    }

    std::vector<unsigned> labelsPerContig;

    for (auto contig : result.stats.includedContigs)
    {
        labelsPerContig.push_back(contigLabels.at(contig));
    }

    return labelsPerContig;
}

void ResultIO::clusteringToYAML(YAML::Emitter & out, const ResultContainer & result, const std::vector<ClusteringResult> & clusts)
{
    out << YAML::BeginSeq;
    for (auto & clust : clusts)
    {
        out << YAML::BeginMap;
        out << YAML::Key << "assignment" << YAML::Value;

        out << YAML::BeginMap
                << YAML::Key << "k" << YAML::Value << clust.numClusters
                << YAML::Key << "labels" << YAML::Value << YAML::Flow << labelsPerContig(result, clust)
                << YAML::Key << "outlierClusters" << YAML::Value << YAML::Flow << clust.outlierClusters
            << YAML::EndMap;

        out << YAML::EndMap;
    }
    out << YAML::EndSeq;
}

void ResultIO::writeYAML(const ResultContainer & result, std::ostream & os)
{
    std::vector<Fixed> sneCol0 = IOUtil::columnToFixed(result.dataSne, 0);
    std::vector<Fixed> sneCol1 = IOUtil::columnToFixed(result.dataSne, 1);
    std::vector<Fixed> pcaCol0 = IOUtil::columnToFixed(result.dataPca, 0);
    std::vector<Fixed> pcaCol1 = IOUtil::columnToFixed(result.dataPca, 1);

    std::vector<unsigned long> contigLengths;
    for (const std::string & c : result.stats.includedContigs) contigLengths.push_back(result.stats.contigLength.at(c));
    std::vector<double> contigGCs;
    for (const std::string & c : result.stats.includedContigs) contigGCs.push_back(result.stats.contigGcContent.at(c));

    std::vector<std::string> krakenContigs;
    for (const std::string & c : result.stats.includedContigs)
    {
        if (result.kraken.classification.find(c) != result.kraken.classification.end())
        {
            krakenContigs.push_back(result.kraken.classification.at(c));
        } else
        {
            krakenContigs.push_back("unknown");
        }
    }

    YAML::Emitter out;
    out << YAML::BeginMap;
    out << YAML::Key << "acdc_parameters" << YAML::Value << Opts::parameters();
    out << YAML::Key << "input_fasta" << YAML::Value << IOUtil::absoluteFilepath(result.fasta);
    out << YAML::Key << "fasta_stats" << YAML::Value
        << YAML::BeginMap
            << YAML::Key << "included_contigs" << YAML::Value << YAML::Flow << result.stats.includedContigs
            << YAML::Key << "excluded_contigs" << YAML::Value << YAML::Flow << result.stats.discardedContigs
            << YAML::Key << "num_basepairs" << YAML::Value << result.stats.numBasepairs
            << YAML::Key << "gc_content" << YAML::Value << result.stats.gcContent
            << YAML::Key << "contig_lengths" << YAML::Value << YAML::Flow << contigLengths
            << YAML::Key << "contig_gc" << YAML::Value << YAML::Flow  << contigGCs
        << YAML::EndMap;
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
            << YAML::Key << "contig_labels" << YAML::Value  << YAML::Flow << contigIndicesVisualization(result)
        << YAML::EndMap;
    out << YAML::Key << "confidence_cc" << YAML::Value << result.contaminationAnalysis.confidenceCC;
    out << YAML::Key << "confidence_dip" << YAML::Value << result.contaminationAnalysis.confidenceDip;
    out << YAML::Key << "contamination_state" << YAML::Value << result.contaminationAnalysis.state;
    out << YAML::Key << "cluster_estimates" << YAML::Value
        << YAML::BeginMap;
            if (result.clusterings.mostLikelyClustering != nullptr)
            {
            out << YAML::Key << "most_likely_clustering" << YAML::Value
                << YAML::BeginMap
                    << YAML::Key << "method" << YAML::Value << result.clusterings.mostLikelyClusteringName
                    << YAML::Key << "estimated_k" << YAML::Value << result.clusterings.mostLikelyClustering->numClusters;
                out << YAML::EndMap;
            }
            out << YAML::Key << "cc" << YAML::Value
                << YAML::BeginMap
                    << YAML::Key << "method" << YAML::Value << result.clusterings.estimatorCC->name()
                    << YAML::Key << "parameters" << YAML::Value << result.clusterings.estimatorCC->parameters()
                    << YAML::Key << "estimated_k" << YAML::Value << (result.clusterings.numClustCC)
                    << YAML::Key << "assignments" << YAML::Value; clusteringToYAML(out, result, result.clusterings.clustsCC);
                out << YAML::EndMap
            << YAML::Key << "validity_sne" << YAML::Value
                << YAML::BeginMap
                    << YAML::Key << "method" << YAML::Value << result.clusterings.estimatorClusterValidity->name()
                    << YAML::Key << "parameters" << YAML::Value << result.clusterings.estimatorClusterValidity->parameters()
                    << YAML::Key << "dimensionality_reduction" << YAML::Value << "tsne"
                    << YAML::Key << "estimated_k" << YAML::Value << result.clusterings.numClustSne
                    << YAML::Key << "assignments" << YAML::Value; clusteringToYAML(out, result, result.clusterings.clustsSne);
                out << YAML::EndMap
            << YAML::Key << "validity_pca" << YAML::Value
                << YAML::BeginMap
                    << YAML::Key << "method" << YAML::Value << result.clusterings.estimatorClusterValidity->name()
                    << YAML::Key << "parameters" << YAML::Value << result.clusterings.estimatorClusterValidity->parameters()
                    << YAML::Key << "dimensionality_reduction" << YAML::Value << "pca"
                    << YAML::Key << "estimated_k" << YAML::Value << result.clusterings.numClustPca
                    << YAML::Key << "assignments" << YAML::Value; clusteringToYAML(out, result, result.clusterings.clustsPca);
                out << YAML::EndMap
        << YAML::EndMap;
    out << YAML::Key << "kraken" << YAML::Value
        << YAML::BeginMap
            << YAML::Key << "enabled" << YAML::Value << krakenEnabled
            << YAML::Key << "bacterial_background" << YAML::Value << result.kraken.bacterialBackground
            << YAML::Key << "classification" << YAML::Value << YAML::Flow << krakenContigs
        << YAML::EndMap;
    out << YAML::Key << "rnammer" << YAML::Value
        << YAML::BeginMap
            << YAML::Key << "sixteen_s_per_point" << YAML::Value << result._16S
        << YAML::EndMap;
    if (result.contaminationAnalysis.state != "clean")
    {
        out << YAML::Key << "taxonomies" << YAML::Value << YAML::BeginMap;
        for (const auto & contig : result.stats.includedContigs)
        {
            if (result.stats.taxonomies.find(contig) != result.stats.taxonomies.end())
            {
                out << YAML::Key << contig << YAML::Value << result.stats.taxonomies.at(contig);
            }
        }
        out << YAML::EndMap;
    }
    out << YAML::EndMap;

    os << out.c_str();
}

void ResultIO::processResult(const ResultContainer & result)
{
    exportContigJS(result.fasta);

    // generate yaml string
    std::stringstream ss;
    writeYAML(result, ss);
    std::string yamlstr = ss.str();

    // write yaml file
    std::string outFileName = result.fasta;
    boost::replace_all(outFileName, "/", "_");
    std::ofstream ofs(outputDir + "/" + outFileName + ".yaml", std::ofstream::out);
    ofs << yamlstr;
    ofs.close();

    // write yaml for javascript
    boost::replace_all(yamlstr, "\n", "\\n");
    std::ofstream ofs2(outputDir + "/" + "clusterinfo.js", std::ofstream::out | std::ofstream::app);
    ofs2 << "clusterinfo['" << result.fasta << "'] = '";
    ofs2 << yamlstr;
    ofs2 << "';" << std::endl;
    ofs2.close();
}
