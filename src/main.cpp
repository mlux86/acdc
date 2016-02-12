#include "SequenceVectorizer.h"
#include "Opts.h"
#include "ClusterAnalysis.h"
#include "Logger.h"
#include "TarjansAlgorithm.h"
#include "BarnesHutSNEAdapter.h"
#include "KrakenAdapter.h"
#include "ResultIO.h"
#include "IOUtil.h"

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
		std::cerr << e.what() << std::endl;
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
		IOUtil::copyDir(boost::filesystem::path("../assets"), outPath, true);
	} catch(const boost::filesystem::filesystem_error & e)
	{
		ELOG << e.what() << "\n";
		return EXIT_FAILURE;
	}

	std::string dataFile = opts->outputDir() + "/data.js";
	if (boost::filesystem::exists(boost::filesystem::path(dataFile)))
	{
		std::remove(dataFile.c_str());
	}

	ResultIO rio(opts->outputDir());
	unsigned idCnt = 1;

	for (const auto & fasta : opts->inputFASTAs())
	{
		if (!boost::filesystem::is_regular_file(boost::filesystem::path(fasta)))
		{
			ELOG << "File '" << fasta << "' does not exist or is not a regular file! Skipping...\n";
			continue;  
		}

		ILOG << "Processing file: " << fasta << "\n";
		ResultContainer result;

		result.id = idCnt++;
		result.fasta = fasta;

		try 
		{
			ILOG << "  Running Kraken ... "; ILOG.flush();
			result.kraken = KrakenAdapter::runKraken(fasta, *opts);

			ILOG << "vectorizing contigs ... "; ILOG.flush();
			SequenceVectorizer sv(fasta, *opts);
			auto dat = sv.vectorize();
			result.fastaLabels = dat.second;

			ILOG << "one-shot analysis ... "; ILOG.flush();
			result.oneshot = ClusterAnalysis::analyze(dat.first, *opts);

			ILOG << "bootstrap analysis ... "; ILOG.flush();
			result.bootstraps = ClusterAnalysis::analyzeBootstraps(dat.first, *opts); 
			
			ILOG << " writing result ... "; ILOG.flush();
			rio.processResult(result);
			ILOG << "done!\n";
		} catch(const std::runtime_error & e)
		{
			ELOG << e.what() << "\n";
			return EXIT_FAILURE;
		}		
	}

	rio.finish();

	return EXIT_SUCCESS;
}