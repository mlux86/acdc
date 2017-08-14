#pragma once

#include <vector>
#include <map>
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
	std::map<unsigned, std::string> _16S;

    // Statistics
    SequenceStats stats;
};

// Processing of result files (i.e. I/O, writing assets, JS/YAML export, etc.)
class ResultIO
{

private:
	// exports fasta contigs to a JS file to enable filtering per cluster
	void exportContigJS(const std::string & fastaFilename);

	// helper: creates a list of contig labels per data point for visualization labels
	std::vector<unsigned> contigIndicesVisualization(const ResultContainer & result);

	// helper: converts a labelling from a ClusteringResult (label per data point) to one label per contig
	std::vector<unsigned> labelsPerContig(const ResultContainer & result, const ClusteringResult & clust);

	// helper: writes a ClusteringResult object to a YAML strean
	void clusteringToYAML(YAML::Emitter & out, const ResultContainer & result, const std::vector<ClusteringResult> & clusts);

	// write YAML output 
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
