#include "SequenceVectorizer.h"
#include "Opts.h"
#include "ClusterAnalysis.h"
#include "Logger.h"
#include "TarjansAlgorithm.h"
#include "BarnesHutSNEAdapter.h"
#include "KrakenAdapter.h"
#include "ResultIO.h"

#include <boost/filesystem.hpp>
#include <string>
#include <fstream>
#include <streambuf>
#include <future>

int main(int argc, char const *argv[])
{

	std::unique_ptr<Opts> opts;

	std::string banner = "";

	try 
	{
		opts.reset(new Opts(argc, argv));

		std::ifstream ifs("../banner.txt"); // TODO dynamic path
		banner = std::string((std::istreambuf_iterator<char>(ifs)),
		                 std::istreambuf_iterator<char>());
	}
	catch(const std::exception & e) 
	{
		ELOG << e.what() << '\n';
		return EXIT_FAILURE;
	}

	if (opts->needsHelp())
	{
		std::cout << banner << std::endl;
		std::cout << opts->helpDesc() << std::endl;
		return EXIT_SUCCESS;
	}

	switch (opts->logLevel())
	{
		case -1: Logger::getInstance().setLevel(Off); break;
		case 0: Logger::getInstance().setLevel(Error); break;
		case 1: Logger::getInstance().setLevel(Info); break;
		case 2: Logger::getInstance().setLevel(Verbose); break;
		case 3: Logger::getInstance().setLevel(Debug); break;
		default: throw std::runtime_error("Loglevel undefined!");
	}

	if (opts->inputFASTAs().empty())
	{
		ELOG << "No input FASTA file(s) given (--input-fasta,-i), aborting.\n";
		return EXIT_FAILURE;
	}


	boost::filesystem::path outPath (opts->outputDir());
	boost::system::error_code returnedError;
	boost::filesystem::create_directories(outPath, returnedError);
	if (returnedError)
	{
		ELOG << "Could not create output directory, aborting.\n";
		return EXIT_FAILURE;
	}
	try 
	{
		Util::copyDir(boost::filesystem::path("../assets"), outPath);
	} catch(const boost::filesystem::filesystem_error & e)
	{
		ELOG << e.what() << "\n";
		return EXIT_FAILURE;
	}

	ResultIO rio(opts->outputDir());

	for (const auto & fasta : opts->inputFASTAs())
	{
		ResultContainer result;

		result.fasta = fasta;

		ILOG << "Running Kraken...\n";
		result.kraken = KrakenAdapter::runKraken(fasta, *opts);

		ILOG << "Vectorizing contigs...\n";
		SequenceVectorizer sv(fasta, *opts);
		auto dat = sv.vectorize();
		result.fastaLabels = dat.second;

		ILOG << "One-shot analysis...\n";
		result.oneshot = ClusterAnalysis::analyze(dat.first, *opts);

		ILOG << "Bootstrap analysis...\n";
		result.bootstraps = ClusterAnalysis::analyzeBootstraps(dat.first, *opts); 

		rio.processResult(result);
	}

	rio.finish();

	return EXIT_SUCCESS;
}