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

std::pair<Eigen::MatrixXd, ClusterAnalysisResult> ClusterAnalysis::analyze(const Eigen::MatrixXd & data, const Opts & opts)
{
	ClusterAnalysisResult res;

	VLOG << "Running t-SNE...\n";
	auto reduced = BarnesHutSNEAdapter::runBarnesHutSNE(data, opts);

	VLOG << "Counting connected components...\n";
	Eigen::MatrixXd affinities = Util::knnAffinityMatrix(reduced, 9, false);
	TarjansAlgorithm ta;
	res.resConnComponents = ta.run(affinities);

	VLOG << "Running dipMeans...\n";
	res.resDipMeans = Clustering::dipMeans(reduced, 0, 0.01, 5);

	return std::make_pair(reduced, res);
}

ClusterAnalysisResult ClusterAnalysis::bootstrapTask(const unsigned taskId, const Eigen::MatrixXd & dataOrig, const Opts & opts, const std::vector<unsigned> indices)
{
	Eigen::MatrixXd data = Eigen::MatrixXd::Zero(indices.size(), dataOrig.cols());
	for (unsigned i = 0; i < indices.size(); i++)
	{
		data.row(i) = dataOrig.row(indices[i]);
	}
	auto res = ClusterAnalysis::analyze(data, opts);

	std::stringstream ss;
	ss << "boot" << taskId;

	VisualizationServer::getInstance().addClustering(ss.str(), res.first, res.second);
	
	return res.second;
}

std::vector<ClusterAnalysisResult> ClusterAnalysis::analyzeBootstraps(const Eigen::MatrixXd & data, const Opts & opts)
{
	const std::vector< std::vector<unsigned> > bootstrapIndices = Util::stratifiedSubsamplingIndices(data.rows(), opts.numBootstraps(), opts.bootstrapRatio());

	ThreadPool pool(opts.numThreads());

	std::vector< std::future<ClusterAnalysisResult> > futures;

	unsigned i = 0;
	for (const auto & indices : bootstrapIndices)
	{
		futures.push_back(
			pool.enqueue(&ClusterAnalysis::bootstrapTask, i++, data, opts, indices)
		);
	}

	std::vector<ClusterAnalysisResult> results;

	for (auto & fut : futures)
	{
		results.push_back(fut.get());
	}

	return results;
}