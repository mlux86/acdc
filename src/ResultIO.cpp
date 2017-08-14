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

    for (unsigned i = 0; i < result.stats.includedContigs.size(); i++)
    {
        if(std::find(contigs16s.begin(), contigs16s.end(), result.stats.includedContigs.at(i)) != contigs16s.end())
        {
            contigIndices.push_back(i);
        }
    }

    return contigIndices;
}

std::vector<unsigned> ResultIO::contigIndicesOfDR(const ResultContainer & result)
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
                << YAML::Key << "k" << YAML::Value << (clust.numClusters + clust.outlierClusters.size())
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
            << YAML::Key << "contig_labels" << YAML::Value  << YAML::Flow << contigIndicesOfDR(result)
        << YAML::EndMap;
    out << YAML::Key << "contigs_16s" << YAML::Value << YAML::Flow << contigsIndicesWith16S(result);
    out << YAML::Key << "confidence_cc" << YAML::Value << result.contaminationAnalysis.confidenceCC;
    out << YAML::Key << "confidence_dip" << YAML::Value << result.contaminationAnalysis.confidenceDip;
    out << YAML::Key << "contamination_state" << YAML::Value << result.contaminationAnalysis.state;
    out << YAML::Key << "cluster_estimates" << YAML::Value 
        << YAML::BeginMap 
            << YAML::Key << "cc" << YAML::Value 
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
    out << YAML::EndMap;

    os << out.c_str();
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

void ResultIO::processResult(const ResultContainer & result)
{
    export16S(result);

    exportContigJS(result.fasta);

    // generate yaml string
    std::stringstream ss;
    writeYAML(result, ss);
    std::string yamlstr = ss.str();

    // write yaml file
    std::ofstream ofs(outputDir + "/" + result.fasta + ".yaml", std::ofstream::out);
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
