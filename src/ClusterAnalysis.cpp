#include "Logger.h"
#include "ThreadPool.h"
#include "ClusterAnalysis.h"
#include "BarnesHutSNEAdapter.h"
#include "MLUtil.h"
#include "Opts.h"
#include "SequenceVectorizer.h"

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

ClusterAnalysisResult ClusterAnalysis::bootstrapTask(const SequenceVectorizationResult & svr, const std::vector<unsigned> indices)
{
	// generate new data matrices from bootstrap indices and cluster them

	Eigen::MatrixXd data = Eigen::MatrixXd::Zero(indices.size(), svr.data.cols());
	std::vector<std::string> contigs(indices.size());

	for (unsigned i = 0; i < indices.size(); i++)
	{
		data.row(i) = svr.data.row(indices[i]);
		contigs[i] = svr.contigs[indices[i]];
	}

	auto res = ClusterAnalysis::analyze(data, contigs, svr.contigSizes);

	res.bootstrapIndexes = indices;

	return res;
}

ClusterAnalysisResult ClusterAnalysis::analyze(const Eigen::MatrixXd & data, const std::vector<std::string> & contigs, const std::map<std::string, unsigned> & contigSizes)
{
	ClusterAnalysisResult res;

	res.dataOrig = data;

	VLOG << "PCA..." << std::endl;
	res.dataPca = MLUtil::pca(data, Opts::tsneDim());

	VLOG << "Running t-SNE..." << std::endl;
	res.dataSne = BarnesHutSNEAdapter::runBarnesHutSNE(data);

	Clustering cSne(res.dataSne, contigs, contigSizes);
	Clustering cPca(res.dataPca, contigs, contigSizes);

	VLOG << "Testing for multi-modality..." << std::endl;
	res.isMultiModal = cSne.isMultiModal(0, 0.001) ||
		               cPca.isMultiModal(0, 0.001);

	VLOG << "Clustering PCA...\n";
	auto resPca = cPca.estimateK(5);
	res.numClustPca = !res.isMultiModal ? 1 : resPca.first;
	res.clustsPca = resPca.second;	

	VLOG << "Clustering t-SNE...\n";
	auto resSne = cSne.estimateK(5);
	res.numClustSne = !res.isMultiModal ? 1 : resSne.first;
	res.clustsSne = resSne.second;

	VLOG << "Clustering connected components" << std::endl;
	res.clustCC = cSne.connComponents(9);
	res.numClustCC = res.clustCC.numClusters;
	res.hasSeparatedComponents = res.numClustCC > 1;

	return res;
}

std::vector<ClusterAnalysisResult> ClusterAnalysis::analyzeBootstraps(const SequenceVectorizationResult & svr)
{
	ThreadPool pool(Opts::numThreads());

	std::vector< std::future<ClusterAnalysisResult> > futures;

	// add oneshot task (always first element in returned vector)
	futures.push_back(
		pool.enqueue(&ClusterAnalysis::analyze, svr.data, svr.contigs, svr.contigSizes)
	);

	// add bootstrap tasks
	const std::vector< std::vector<unsigned> > bootstrapIndexes = ClusterAnalysis::stratifiedSubsamplingIndices(svr.data.rows(), Opts::numBootstraps(), Opts::bootstrapRatio());
	for (const auto & indices : bootstrapIndexes)
	{
		futures.push_back(
			pool.enqueue(&ClusterAnalysis::bootstrapTask, svr, indices)
		);
	}

	std::vector<ClusterAnalysisResult> results;

	for (auto & fut : futures)
	{
		results.push_back(fut.get());
	}

	return results;
}