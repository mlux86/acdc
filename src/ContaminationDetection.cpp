#include <random>
#include "Logger.h"
#include "ThreadPool.h"
#include "ClusteringUtil.h"
#include "ContaminationDetection.h"
#include "BarnesHutSNEAdapter.h"
#include "MLUtil.h"
#include "Opts.h"
#include "SequenceVectorizer.h"

ContaminationDetection::ContaminationDetection()
{

}

ContaminationDetection::~ContaminationDetection()
{

}

std::vector< std::vector<unsigned> > ContaminationDetection::stratifiedSubsamplingIndices(const unsigned n, const unsigned k, const double ratio)
{
    std::vector<unsigned> indices(n);
    std::iota(indices.begin(), indices.end(), 0);

	auto rng = std::default_random_engine {};
    std::shuffle(indices.begin(), indices.end(), rng);

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

ContaminationDetectionResult ContaminationDetection::bootstrapTask(const SequenceVectorizationResult & svr, const std::vector<unsigned> indices, const unsigned bootstrapId)
{
	// generate new data matrices from bootstrap indices and cluster them

	Eigen::MatrixXd data = Eigen::MatrixXd::Zero(indices.size(), svr.data.cols());
	std::vector<std::string> contigs(indices.size());

	for (unsigned i = 0; i < indices.size(); i++)
	{
		data.row(i) = svr.data.row(indices[i]);
		contigs[i] = svr.contigs[indices[i]];
	}

	auto res = ContaminationDetection::analyze(data, contigs, svr.contigSizes, Opts::randomSeed() + bootstrapId);

	res.bootstrapIndexes = indices;

	return res;
}

ContaminationDetectionResult ContaminationDetection::analyze(const Eigen::MatrixXd & data, const std::vector<std::string> & contigs, const std::map<std::string, unsigned> & contigSizes, const unsigned randomSeed)
{
	ContaminationDetectionResult res;

	res.dataOrig = data;

	VLOG << "PCA..." << std::endl;
	res.dataPca = MLUtil::pca(data, Opts::tsneDim());

	VLOG << "Running t-SNE..." << std::endl;
	res.dataSne = BarnesHutSNEAdapter::runBarnesHutSNE(data, randomSeed);

	VLOG << "Testing for multi-modality..." << std::endl;
	res.isMultiModal = ClusteringUtil::isMultiModal(res.dataSne, 0, 0.001) ||
		               ClusteringUtil::isMultiModal(res.dataPca, 0, 0.001);

	VLOG << "Clustering connected components" << std::endl;
	std::unique_ptr<ClusterPostProcessing> cpp(new ClusterPostProcessing(contigs, contigSizes));
	ConnectedComponentsEstimator cce(cpp, 9);
	auto clustCC = cce.estimateK(res.dataSne);
	res.hasSeparatedComponents = clustCC.first > 1;

	return res;
}

std::vector<ContaminationDetectionResult> ContaminationDetection::analyzeBootstraps(const SequenceVectorizationResult & svr)
{
	ThreadPool pool(Opts::numThreads());

	std::vector< std::future<ContaminationDetectionResult> > futures;

	// add oneshot task (always first element in returned vector)
	futures.push_back(
		pool.enqueue(&ContaminationDetection::analyze, svr.data, svr.contigs, svr.contigSizes, Opts::randomSeed())
	);

	// add bootstrap tasks
	const std::vector< std::vector<unsigned> > bootstrapIndexes = ContaminationDetection::stratifiedSubsamplingIndices(svr.data.rows(), Opts::numBootstraps(), Opts::bootstrapRatio());
	unsigned bootstrapId = 1;
	for (const auto & indices : bootstrapIndexes)
	{
		futures.push_back(
			pool.enqueue(&ContaminationDetection::bootstrapTask, svr, indices, bootstrapId)
		);
		bootstrapId++;
	}

	std::vector<ContaminationDetectionResult> results;

	for (auto & fut : futures)
	{
		results.push_back(fut.get());
	}

	return results;
}

ContaminationDetectionSummary ContaminationDetection::summarizeBootstraps(const std::vector<ContaminationDetectionResult> & bootstraps)
{
	ContaminationDetectionSummary cds; 

	cds.confidenceDip = 0;
	cds.confidenceCC = 0;

	for (auto & bs : bootstraps)
	{
		if (bs.isMultiModal)
		{
			cds.confidenceDip++;
		}

		if (bs.hasSeparatedComponents)
		{
			cds.confidenceCC++;
		}		
	}

	cds.confidenceDip /= bootstraps.size();
	cds.confidenceCC /= bootstraps.size();

	if (cds.confidenceDip < 0.25 && cds.confidenceCC < 0.25)
	{
		cds.state = "clean";
	} else if (cds.confidenceDip > 0.75 || cds.confidenceCC > 0.75)
	{
		cds.state = "contaminated";
	} else
	{
		cds.state = "warning";
	}

	return(cds);
}