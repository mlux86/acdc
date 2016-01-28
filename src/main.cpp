#include "SequenceVectorizer.h"
#include "Opts.h"
#include "ClusterAnalysis.h"
#include "Logger.h"
#include "VisualizationServer.h"
#include "TarjansAlgorithm.h"
#include "BarnesHutSNEAdapter.h"
#include "KrakenAdapter.h"

#include "Controller.h"
#include "WebServer.h"

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

	const unsigned port = opts->port();
	std::thread vst([port]()
	{
		VisualizationServer::getInstance().run(port);
	});


	for (const auto & fasta : opts->inputFASTAs())
	{
		ILOG << "Running Kraken...\n";
		auto krakenResult = KrakenAdapter::runKraken(fasta, *opts);

		ILOG << "Vectorizing contigs...\n";
		SequenceVectorizer sv(fasta, *opts);
		auto dat = sv.vectorize();

		ILOG << "One-shot analysis...\n";
		auto res = ClusterAnalysis::analyze(dat.first, *opts);
		auto dataPca = Util::pca(dat.first, opts->tsneDim());
		VisualizationServer::getInstance().addClustering(fasta, true, dataPca, res.first, dat.second, res.second);
		VisualizationServer::getInstance().addKrakenResult(fasta, krakenResult);

		ILOG << "Bootstrap analysis...\n";
		auto results = ClusterAnalysis::analyzeBootstraps(fasta, dat.first, dat.second, *opts);
	}

	// ILOG << "Waiting for visualization server to stop...\n";
	// VisualizationServer::getInstance().stop();
	vst.join();

	return EXIT_SUCCESS;
}