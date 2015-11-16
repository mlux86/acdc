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

ClusterAnalysisResult ClusterAnalysis::analyze(const Eigen::MatrixXd & data, const Opts & opts, const std::string & descPrefix)
{
	VLOG << "Running t-SNE...\n";
	Eigen::MatrixXd reduced = BarnesHutSNEAdapter::runBarnesHutSNE(data, opts);
	
	ClusterAnalysisResult res;

	VLOG << "Counting connected components...\n";
	Eigen::MatrixXd affinities = Util::knnAffinityMatrix(reduced, 9, false);
	TarjansAlgorithm ta;
	res.resConnComponents = ta.run(affinities);

	VLOG << "Running dipMeans...\n";
	res.resDipMeans = Clustering::dipMeans(reduced, 0, 0.01, 5);

	VisualizationServer::getInstance().addClustering(reduced, res.resConnComponents, descPrefix + "-cc");
	VisualizationServer::getInstance().addClustering(reduced, res.resDipMeans, descPrefix + "-dm");

	return res;
}

ClusterAnalysisResult ClusterAnalysis::bootstrapTask(const unsigned taskId, const Eigen::MatrixXd & dataOrig, const Opts & opts, const std::vector<unsigned> indices)
{
	Eigen::MatrixXd data = Eigen::MatrixXd::Zero(indices.size(), dataOrig.cols());
	for (unsigned i = 0; i < indices.size(); i++)
	{
		data.row(i) = dataOrig.row(indices[i]);
	}
	std::stringstream ss;
	ss << "bootstrap" << taskId;
	return ClusterAnalysis::analyze(data, opts, ss.str());
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