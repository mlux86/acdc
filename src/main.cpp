#include "easylogging++.h"
INITIALIZE_EASYLOGGINGPP

#include "SequenceVectorizer.h"
#include "BarnesHutSNEAdapter.h"
#include "Util.h"
#include "Opts.h"
#include "Clustering.h"
#include "TarjansAlgorithm.h"

// logging: VERBOSE 1: progress
//          VERBOSE 2: parameters and more than progress
//          VERBOSE 3: t-SNE and more

int main(int argc, char const *argv[])
{
	START_EASYLOGGINGPP(argc, argv);	

	// std::unique_ptr<Opts> opts;

	// try 
	// {
	// 	opts.reset(new Opts(argc, argv));
	// }
	// catch(const std::exception & e) 
	// {
	// 	std::cerr << e.what() << '\n';
	// 	return EXIT_FAILURE;
	// }

	// if (opts->needsHelp())
	// {
	// 	std::cout << opts->helpDesc() << std::endl;
	// 	return EXIT_SUCCESS;
	// }

	// if (opts->inputFASTA().empty())
	// {
	// 	std::cerr << "No input FASTA file given (--input-fasta,-i), aborting." << std::endl;
	// 	return EXIT_FAILURE;
	// }


	// VLOG(1) << "Vectorizing contigs...";
	// SequenceVectorizer sv(*opts);
	// auto dat = sv.vectorize();

	// VLOG(1) << "Running t-SNE...";
	// Eigen::MatrixXd reduced = BarnesHutSNEAdapter::runBarnesHutSNE(dat.first, *opts);
	
	// Util::saveMatrix(reduced, "/tmp/reduced", ' ');

	// VLOG(1) << "Counting connected components...";
	// Eigen::MatrixXd affinities = Util::knnAffinityMatrix(reduced, 9, false);
	// TarjansAlgorithm ta;
	// ClusteringResult res = ta.run(affinities);
	// VLOG(1) << "Found " << res.numClusters << " clusters.";


	Eigen::MatrixXd mat = Util::loadMatrix("/tmp/meh", ' ');

	ClusteringResult kmeansRes = Clustering::kMeans(mat, 2);

	Util::saveMatrix(kmeansRes.labels, "/tmp/labels", ' ');

	return EXIT_SUCCESS;
}