#include "crow.h"

#include "SequenceVectorizer.h"
#include "Opts.h"
#include "ClusterAnalysis.h"
#include "Logger.h"

#include <string>
#include <fstream>
#include <streambuf>

int main(int argc, char const *argv[])
{
	crow::SimpleApp app;
	CROW_ROUTE(app, "/")([](){
        return "Hello world";
    });	
    app.port(18080).multithreaded();

	std::unique_ptr<Opts> opts;

	std::string banner = "";

	try 
	{
		opts.reset(new Opts(argc, argv));

		std::ifstream ifs("../banner.txt");
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

	if (opts->inputFASTA().empty())
	{
		ELOG << "No input FASTA file given (--input-fasta,-i), aborting.\n";
		return EXIT_FAILURE;
	}

	ILOG << "Vectorizing contigs...\n";
	SequenceVectorizer sv(*opts);
	auto dat = sv.vectorize();

	ILOG << "Starting bootstrap analysis...\n";
	auto results = ClusterAnalysis::analyzeBootstraps(dat.first, *opts);

	for (const auto & res : results)
	{
		ILOG << "CC k=" << res.resConnComponents.numClusters << "\n";
		ILOG << "DM k=" << res.resDipMeans.numClusters << "\n"; 
	}

	return EXIT_SUCCESS;
}