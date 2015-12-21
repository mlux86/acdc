#include "Logger.h"
#include "ThreadPool.h"
#include "ClusterAnalysis.h"
#include "BarnesHutSNEAdapter.h"
#include "TarjansAlgorithm.h"
#include "Util.h"
#include "VisualizationServer.h"

ClusterAnalysis::ClusterAnalysis()
{

}

ClusterAnalysis::~ClusterAnalysis()
{

}

ClusterAnalysisResult ClusterAnalysis::bootstrapTask(const unsigned taskId, const Eigen::MatrixXd & dataOrig, const std::vector<std::string> & labelsOrig, const Opts & opts, const std::vector<unsigned> indices)
{
	Eigen::MatrixXd data = Eigen::MatrixXd::Zero(indices.size(), dataOrig.cols());
	std::vector<std::string> labels(indices.size());	
	for (unsigned i = 0; i < indices.size(); i++)
	{
		data.row(i) = dataOrig.row(indices[i]);
		labels[i] = labelsOrig[indices[i]];
	}
	auto res = ClusterAnalysis::analyze(data, opts);

	auto dataPca = Util::pca(data, opts.tsneDim());
	
	std::stringstream ss;
	ss << "boot" << taskId;

	VisualizationServer::getInstance().addClustering(ss.str(), dataPca, res.first, labels, res.second);
	
	return res.second;
}

std::pair<Eigen::MatrixXd, ClusterAnalysisResult> ClusterAnalysis::analyze(const Eigen::MatrixXd & data, const Opts & opts)
{
	ClusterAnalysisResult res;

	VLOG << "Running t-SNE...\n";
	auto datSne = BarnesHutSNEAdapter::runBarnesHutSNE(data, opts);

	VLOG << "Counting connected components...\n";
	Eigen::MatrixXd affinities = Util::knnAffinityMatrix(datSne, 7, false);
	TarjansAlgorithm ta;
	res.resConnComponents = ta.run(affinities);

	VLOG << "Running dipMeans...\n";
	// res.resDipMeans = Clustering::dipMeans(datSne, 0, 0.01, 5);

	return std::make_pair(datSne, res);
}

std::vector<ClusterAnalysisResult> ClusterAnalysis::analyzeBootstraps(const Eigen::MatrixXd & data, const std::vector<std::string> & labels, const Opts & opts)
{
	const std::vector< std::vector<unsigned> > bootstrapIndices = Util::stratifiedSubsamplingIndices(data.rows(), opts.numBootstraps(), opts.bootstrapRatio());

	ThreadPool pool(opts.numThreads());

	std::vector< std::future<ClusterAnalysisResult> > futures;

	unsigned i = 0;
	for (const auto & indices : bootstrapIndices)
	{
		futures.push_back(
			pool.enqueue(&ClusterAnalysis::bootstrapTask, i++, data, labels, opts, indices)
		);
	}

	std::vector<ClusterAnalysisResult> results;

	for (auto & fut : futures)
	{
		results.push_back(fut.get());
	}

	return results;
}