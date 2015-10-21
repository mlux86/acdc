#include <iostream>

#include "SequenceVectorizer.h"
#include "BarnesHutSNEBridge.h"
#include "Util.h"
#include "Opts.h"

int main(int argc, char const *argv[])
{

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


	SequenceVectorizer sv(opts->windowKmerLength(), opts->windowWidth(), opts->windowStep());

	auto dat = sv.vectorize(opts->inputFASTA());

	auto reduced = BarnesHutSNEBridge::runBarnesHutSNE(dat.first, *opts);

	Util::saveMatrix(reduced, "/tmp/reduced", ' ');

	return 0;
}