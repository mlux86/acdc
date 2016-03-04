#include "Logger.h"
#include "ThreadPool.h"
#include "ClusterAnalysis.h"
#include "BarnesHutSNEAdapter.h"
#include "TarjansAlgorithm.h"
#include "MatrixUtil.h"
#include "MLUtil.h"

ClusterAnalysis::ClusterAnalysis()
{

}

ClusterAnalysis::~ClusterAnalysis()
{

}

std::vector< std::vector<unsigned> > ClusterAnalysis::stratifiedSubsamplingIndices(const unsigned n, const unsigned k, const double ratio)
{
    std::vector<unsigned> indices(n);
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), std::mt19937{std::random_device{}()});

    unsigned m = (unsigned) ((double)n * ratio);
    unsigned step = (n - m) / k;

    std::vector< std::vector<unsigned> > result;

    for (unsigned i = 0; i < k; i++)
    {
        unsigned from = i*step;
        unsigned to = i*step + m;
        if (i == k-1)
        {
            to = std::max(to, n);
        }

        std::vector<unsigned> e(indices.begin()+from, indices.begin()+to);
        result.push_back(e); 
    }

    return result;
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

	res.dataOrig = data;

	VLOG << "PCA...\n";
	res.dataPca = MLUtil::pca(data, opts.tsneDim());

	VLOG << "Running t-SNE...\n";
	res.dataSne = BarnesHutSNEAdapter::runBarnesHutSNE(data, opts);

	VLOG << "Counting connected components...\n";
	Eigen::MatrixXd affinities = MLUtil::knnAffinityMatrix(res.dataSne, 7, false);
	TarjansAlgorithm ta;
	res.resConnComponents = ta.run(affinities);

	VLOG << "Running dipMeans...\n";
	res.resDipMeans = Clustering::dipMeans(res.dataSne, 0, 0.01, 5);

	return res;
}

std::vector<ClusterAnalysisResult> ClusterAnalysis::analyzeBootstraps(const Eigen::MatrixXd & data, const Opts & opts)
{
	ThreadPool pool(opts.numThreads());

	std::vector< std::future<ClusterAnalysisResult> > futures;

	// add oneshot task
	futures.push_back(
		pool.enqueue(&ClusterAnalysis::analyze, data, opts)
	);

	// add bootstrap tasks
	const std::vector< std::vector<unsigned> > bootstrapIndexes = ClusterAnalysis::stratifiedSubsamplingIndices(data.rows(), opts.numBootstraps(), opts.bootstrapRatio());
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