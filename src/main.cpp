#include "SequenceVectorizer.h"
#include "Opts.h"
#include "ClusterAnalysis.h"
#include "Logger.h"
#include "TarjansAlgorithm.h"
#include "BarnesHutSNEAdapter.h"
#include "KrakenAdapter.h"
#include "RnammerAdapter.h"
#include "ResultIO.h"
#include "IOUtil.h"
#include "MLUtil.h"
#include "HierarchicalClustering.h"

#include <boost/filesystem.hpp>
#include <string>
#include <fstream>
#include <streambuf>
#include <future>
#include <map>

#include <eigen3/Eigen/Dense>
#include "MatrixUtil.h"

int main(int argc, char *argv[])
{
	std::unique_ptr<Opts> opts;

	std::string banner = "";

	// build program arguments

	try
	{
		opts.reset(new Opts(argc, argv));

		std::ifstream ifs(opts->sharePath() + "/banner.txt");
		banner = std::string((std::istreambuf_iterator<char>(ifs)),
		                 std::istreambuf_iterator<char>());
	}
	catch(const std::exception & e)
	{
		std::cerr << e.what() << std::endl;
		return EXIT_FAILURE;
	}

	// show help

	if (opts->needsHelp())
	{
		std::cout << banner << std::endl;
		std::cout << opts->helpDesc() << std::endl;
		return EXIT_SUCCESS;
	}

	// configure logger

	switch (opts->logLevel())
	{
		case -1: Logger::getInstance().setLevel(Off); break;
		case 0: Logger::getInstance().setLevel(Error); break;
		case 1: Logger::getInstance().setLevel(Info); break;
		case 2: Logger::getInstance().setLevel(Verbose); break;
		case 3: Logger::getInstance().setLevel(Debug); break;
		default: throw std::runtime_error("Loglevel undefined!");
	}

	// check for input files

	if (opts->inputFASTAs().empty())
	{
		ELOG << "No input FASTA file(s) given (--input-fasta,-i)." << std::endl;
		return EXIT_FAILURE;
	}

	// setup Kraken

	KrakenAdapter krk(*opts);
	bool krakenExists = krk.krakenExists();
	if (!krakenExists)
	{
		ELOG << "Kraken not found! It will be disabled.\n- Please make sure that the folders containing the 'kraken' and 'kraken-translate' executables is in your $PATH.\n- Please make sure to supply a database using the --kraken-db switch." << std::endl;
	}

	// setup Rnammer 
	bool rmrExists = RnammerAdapter::rnammerExists();
	if (!rmrExists)
	{
		ELOG << "Rnammer not found! It will be disabled.\n- Please make sure that the folders containing the 'rnammer' executable is in your $PATH." << std::endl;
	}

	// create output directory
	boost::filesystem::path outPath (opts->outputDir());
	boost::system::error_code returnedError;
	boost::filesystem::create_directories(outPath, returnedError);
	if (returnedError)
	{
		ELOG << "Could not create output directory, aborting." << std::endl;
		return EXIT_FAILURE;
	}

	// copy result assets

	try
	{
		IOUtil::copyDir(boost::filesystem::path(opts->sharePath() + "/assets"), outPath, true);
	} catch(const boost::filesystem::filesystem_error & e)
	{
		ELOG << e.what() << std::endl;
		return EXIT_FAILURE;
	}

	// remove old data.js file if necessary

	std::string dataFile = opts->outputDir() + "/data.js";
	if (boost::filesystem::exists(boost::filesystem::path(dataFile)))
	{
		std::remove(dataFile.c_str());
	}

    // process files

	ResultIO rio(opts->outputDir(), krakenExists);
	unsigned idCnt = 1;

	for (const auto & fasta : opts->inputFASTAs())
	{
		if (!boost::filesystem::is_regular_file(boost::filesystem::path(fasta)))
		{
			ELOG << "File '" << fasta << "' does not exist or is not a regular file! Skipping..." << std::endl;
			continue;
		}

		ILOG << "Processing file: " << fasta << std::endl;
		ResultContainer result;

		result.id = idCnt++;
		result.fasta = fasta;

		try
		{
			if (krakenExists)
			{
				ILOG << "Running Kraken..." << std::endl;
				result.kraken = krk.runKraken(fasta);
			}

			ILOG << "Vectorizing contigs..." << std::endl;
			SequenceVectorizer sv(fasta, *opts);
			auto svr = sv.vectorize();
			result.fastaLabels = svr.contigs;

			if (rmrExists)
			{
				ILOG << "Running Rnammer..." << std::endl;
				result.contains16S = RnammerAdapter::mark16S(fasta, svr);
			}

			ILOG << "Clustering..." << std::endl;
			result.bootstraps = ClusterAnalysis::analyzeBootstraps(svr, *opts);
			result.oneshot = result.bootstraps.at(0);
			result.bootstraps.erase(result.bootstraps.begin());

			rio.processResult(result);
		} catch(const std::exception & e)
		{
			ELOG << "An error occurred processing this file: " << e.what() << std::endl;
		}
	}

	return EXIT_SUCCESS;
}
