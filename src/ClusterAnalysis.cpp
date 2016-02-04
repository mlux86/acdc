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

ClusterAnalysisResult ClusterAnalysis::bootstrapTask(const Eigen::MatrixXd & dataOrig, const Opts & opts, const std::vector<unsigned> indices)
{
	Eigen::MatrixXd data = Eigen::MatrixXd::Zero(indices.size(), dataOrig.cols());

	for (unsigned i = 0; i < indices.size(); i++)
	{
		data.row(i) = dataOrig.row(indices[i]);
	}

	auto res = ClusterAnalysis::analyze(data, opts);

	res.bootstrapIndexes = indices;

	return res;
}

ClusterAnalysisResult ClusterAnalysis::analyze(const Eigen::MatrixXd & data, const Opts & opts)
{
	ClusterAnalysisResult res;

	VLOG << "PCA...\n";
	res.dataPca = Util::pca(data, opts.tsneDim());

	VLOG << "Running t-SNE...\n";
	res.dataSne = BarnesHutSNEAdapter::runBarnesHutSNE(data, opts);

	VLOG << "Counting connected components...\n";
	Eigen::MatrixXd affinities = Util::knnAffinityMatrix(res.dataSne, 7, false);
	TarjansAlgorithm ta;
	res.resConnComponents = ta.run(affinities);

	VLOG << "Running dipMeans...\n";
	// res.resDipMeans = Clustering::dipMeans(datSne, 0, 0.01, 5);

	return res;
}

std::vector<ClusterAnalysisResult> ClusterAnalysis::analyzeBootstraps(const Eigen::MatrixXd & data, const Opts & opts)
{
	const std::vector< std::vector<unsigned> > bootstrapIndexes = Util::stratifiedSubsamplingIndices(data.rows(), opts.numBootstraps(), opts.bootstrapRatio());

	ThreadPool pool(opts.numThreads());

	std::vector< std::future<ClusterAnalysisResult> > futures;

	for (const auto & indices : bootstrapIndexes)
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