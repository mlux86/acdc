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

ContaminationDetectionResult ContaminationDetection::bootstrapTask(const SequenceVectorizationResult & svr, const std::vector<unsigned> indices)
{
	// generate new data matrices from bootstrap indices and cluster them

	Eigen::MatrixXd data = Eigen::MatrixXd::Zero(indices.size(), svr.data.cols());
	std::vector<std::string> contigs(indices.size());

	for (unsigned i = 0; i < indices.size(); i++)
	{
		data.row(i) = svr.data.row(indices[i]);
		contigs[i] = svr.contigs[indices[i]];
	}

	auto res = ContaminationDetection::analyze(data, contigs, svr.contigSizes);

	res.bootstrapIndexes = indices;

	return res;
}

ContaminationDetectionResult ContaminationDetection::analyze(const Eigen::MatrixXd & data, const std::vector<std::string> & contigs, const std::map<std::string, unsigned> & contigSizes)
{
	ContaminationDetectionResult res;

	res.dataOrig = data;

	VLOG << "PCA..." << std::endl;
	res.dataPca = MLUtil::pca(data, Opts::tsneDim());

	VLOG << "Running t-SNE..." << std::endl;
	res.dataSne = BarnesHutSNEAdapter::runBarnesHutSNE(data);

	VLOG << "Testing for multi-modality..." << std::endl;
	res.isMultiModal = ClusteringUtil::isMultiModal(res.dataSne, 0, 0.001) ||
		               ClusteringUtil::isMultiModal(res.dataPca, 0, 0.001);

	VLOG << "Clustering connected components" << std::endl;
	Clustering cSne(res.dataSne, contigs, contigSizes);
	auto clustCC = cSne.connComponents(9);
	res.hasSeparatedComponents = clustCC.numClusters > 1;

	return res;
}

std::vector<ContaminationDetectionResult> ContaminationDetection::analyzeBootstraps(const SequenceVectorizationResult & svr)
{
	ThreadPool pool(Opts::numThreads());

	std::vector< std::future<ContaminationDetectionResult> > futures;

	// add oneshot task (always first element in returned vector)
	futures.push_back(
		pool.enqueue(&ContaminationDetection::analyze, svr.data, svr.contigs, svr.contigSizes)
	);

	// add bootstrap tasks
	const std::vector< std::vector<unsigned> > bootstrapIndexes = ContaminationDetection::stratifiedSubsamplingIndices(svr.data.rows(), Opts::numBootstraps(), Opts::bootstrapRatio());
	for (const auto & indices : bootstrapIndexes)
	{
		futures.push_back(
			pool.enqueue(&ContaminationDetection::bootstrapTask, svr, indices)
		);
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
		cds.status = "clean";
	} else if (cds.confidenceDip > 0.75 || cds.confidenceCC > 0.75)
	{
		cds.status = "contaminated";
	} else
	{
		cds.status = "warning";
	}

	return(cds);
}