#include "easylogging++.h"
INITIALIZE_EASYLOGGINGPP

#include "SequenceVectorizer.h"
#include "Opts.h"
#include "ClusterAnalysis.h"

// logging: VERBOSE 1: progress
//          VERBOSE 2: parameters and more than progress
//          VERBOSE 3: t-SNE and more

int main(int argc, char const *argv[])
{
	START_EASYLOGGINGPP(argc, argv);	

	std::unique_ptr<Opts> opts;

	try 
	{
		opts.reset(new Opts(argc, argv));
	}
	catch(const std::exception & e) 
	{
		std::cerr << e.what() << '\n';
		return EXIT_FAILURE;
	}

	if (opts->needsHelp())
	{
		std::cout << opts->helpDesc() << std::endl;
		return EXIT_SUCCESS;
	}

	if (opts->inputFASTA().empty())
	{
		std::cerr << "No input FASTA file given (--input-fasta,-i), aborting." << std::endl;
		return EXIT_FAILURE;
	}


	VLOG(1) << "Vectorizing contigs...";
	SequenceVectorizer sv(*opts);
	auto dat = sv.vectorize();

	// ClusterAnalysisResult car = ClusterAnalysis::analyze(dat.first, *opts);
	auto results = ClusterAnalysis::analyzeBootstraps(dat.first, *opts, 8);

	for (const auto & res : results)
	{
		VLOG(1) << "CC k=" << res.resConnComponents.numClusters;
		VLOG(1) << "DM k=" << res.resDipMeans.numClusters; 
	}

	return EXIT_SUCCESS;
}