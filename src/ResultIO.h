#pragma once

#include <vector>
#include <string>
#include <yaml-cpp/yaml.h>

#include "ContaminationDetection.h"
#include "Decontamination.h"
#include "KrakenAdapter.h"
#include "SequenceUtil.h"
#include "Clustering.h"

// forward declarations
struct ClusteringResult;
struct ContaminationDetectionSummary;

// Contains the full result of one sample analysis
struct ResultContainer
{
	// Unique ID of the sample
	unsigned id;

	// Fasta file analyzed
	std::string fasta;

	// Vectorized representation of result
	SequenceVectorizationResult seqVectorization;

	// 2D t-SNE representation of the data
	Eigen::MatrixXd dataSne;
	
	// 2D PCA representation of the data
	Eigen::MatrixXd dataPca;

	// Result of contamination analysis
	ContaminationDetectionSummary contaminationAnalysis;

	// Result of clusterings for decontamination
	DecontaminationResult clusterings;

	// Kraken results
	KrakenResult kraken;

	// 16S genes
	// _16S[i].empty() => no 16S gene found in oneshot data point i
	// otherwise       => 16S sequence _16s[i] found in oneshot data point i
	std::vector<std::string> _16S;

    // Statistics
    SequenceStats stats;
};

// Processing of result files (i.e. I/O, writing assets, JS/YAML export, etc.)
class ResultIO
{

private:
	// Converts a vector of strings to a vector unsigned where each unique string corresponds to a unique unsigned value
	std::vector<unsigned> numericLabels(const std::vector<std::string> & labels);

	// Write 16S sequences
	void export16S(const ResultContainer & result);

	// Export fasta files for cluster exporting
    void exportClusteringInfo(const ResultContainer & result, const std::string & filename);

	void exportContigJS(const std::string & fastaFilename);

	std::vector<unsigned> contigIndicesOfDR(const ResultContainer & result);

	std::vector<unsigned> labelsPerContig(const ResultContainer & result, const ClusteringResult & clust);

	std::vector<unsigned> contigsIndicesWith16S(const ResultContainer & result);

	void clusteringToYAML(YAML::Emitter & out, const ResultContainer & result, const std::vector<ClusteringResult> & clusts);

	void writeYAML(const ResultContainer & result, std::ostream & os);

    // Output directory to write to
	std::string outputDir;

	// is Kraken enabled?
    bool krakenEnabled;

public:
	ResultIO(const std::string & outputDir, bool krakenEnabled);
	~ResultIO();

	// Processes a ResultContainer object to serialize it and write other information
	void processResult(const ResultContainer & result);

};
