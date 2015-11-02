#include "easylogging++.h"

#include "ThreadPool.h"
#include "ClusterAnalysis.h"
#include "BarnesHutSNEAdapter.h"
#include "TarjansAlgorithm.h"
#include "Util.h"

ClusterAnalysis::ClusterAnalysis()
{

}

ClusterAnalysis::~ClusterAnalysis()
{

}

ClusterAnalysisResult ClusterAnalysis::analyze(const Eigen::MatrixXd & data, const Opts & opts)
{
	VLOG(1) << "Running t-SNE...";
	Eigen::MatrixXd reduced = BarnesHutSNEAdapter::runBarnesHutSNE(data, opts);
	
	ClusterAnalysisResult res;

	VLOG(1) << "Counting connected components...";
	Eigen::MatrixXd affinities = Util::knnAffinityMatrix(reduced, 9, false);
	TarjansAlgorithm ta;
	res.resConnComponents = ta.run(affinities);

	VLOG(1) << "Running dipMeans...";
	res.resDipMeans = Clustering::dipMeans(reduced);

	return res;
}

ClusterAnalysisResult ClusterAnalysis::bootstrapTask(const Eigen::MatrixXd & dataOrig, const Opts & opts, const std::vector<unsigned> indices)
{
	Eigen::MatrixXd data = Eigen::MatrixXd::Zero(indices.size(), dataOrig.cols());
	for (unsigned i = 0; i < indices.size(); i++)
	{
		data.row(i) = dataOrig.row(indices[i]);
	}
	return ClusterAnalysis::analyze(data, opts);
}

std::vector<ClusterAnalysisResult> ClusterAnalysis::analyzeBootstraps(const Eigen::MatrixXd & data, const Opts & opts, const unsigned numBootstraps)
{
	const std::vector< std::vector<unsigned> > bootstrapIndices = Util::stratifiedSubsamplingIndices(data.rows(), numBootstraps, 0.8);

	unsigned concurentThreadsSupported = std::thread::hardware_concurrency();
	if (concurentThreadsSupported == 0)
	{
		VLOG(2) << "Could not detect number of cores. Defaulting to one thread.";
		concurentThreadsSupported = 1;
	} else
	{
		VLOG(2) << "Detected " << concurentThreadsSupported << " cores.";
	}

	ThreadPool pool(concurentThreadsSupported);

	std::vector< std::future<ClusterAnalysisResult> > futures;

	for (const auto & indices : bootstrapIndices)
	{
		futures.push_back(
			pool.enqueue(&ClusterAnalysis::bootstrapTask, data, opts, indices)
		);
	}

	std::vector<ClusterAnalysisResult> results;

	for (auto & fut : futures)
	{
		results.push_back(fut.get());
	}

	return results;
}