#include "Opts.h"

#include <boost/program_options.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <sstream>

Opts::Opts(int argc, char const *argv[])
{
	initialize(argc, argv);
}

Opts::~Opts()
{
}

void Opts::initialize(int argc, char const *argv[])
{
	boost::program_options::options_description description("Usage");

	description.add_options()
	    ("help,h", "Display this help message")
	    ("version,v", "Display the version number")
	    ("input-fasta,i", boost::program_options::value<std::string>()->default_value(""), "Input FASTA file")
	    ("tsne-dimension,d", boost::program_options::value<unsigned>()->default_value(2), "T-SNE dimension")
	    ("tsne-perplexity,p", boost::program_options::value<unsigned>(), "T-SNE perplexity (overrides automatic estimation)")
	    ("tsne-theta,t", boost::program_options::value<float>()->default_value(0.5), "T-SNE parameter 'theta' of the underlying Barnes-Hut approximation")
	    ("window-kmer-length,k", boost::program_options::value<unsigned>()->default_value(4), "Length of the k-mers in the sequence vectorizer window")
	    ("window-width,w", boost::program_options::value<unsigned>()->default_value(1000), "Width of the sliding sequence vectorizer window")
	    ("window-step,s", boost::program_options::value<unsigned>()->default_value(500), "Step of the sliding sequence vectorizer window")
	    ;
	
	boost::program_options::variables_map vm;
	boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(description).run(), vm);
	boost::program_options::notify(vm);

	std::stringstream ss;
	ss << description;
	_helpDesc = ss.str();

	_needsHelp = vm.count("help");
	_inputFASTA = vm["input-fasta"].as<std::string>();
	boost::algorithm::trim(_inputFASTA);
	_tsneDim = vm["tsne-dimension"].as<unsigned>();
	if (vm.count("tsne-perplexity"))
	{
		_tsnePerplexity = vm["tsne-perplexity"].as<unsigned>();
	}
	_tsneTheta = vm["tsne-theta"].as<float>();
	_windowKmerLength = vm["window-kmer-length"].as<unsigned>();
	_windowWidth = vm["window-width"].as<unsigned>();
	_windowStep = vm["window-step"].as<unsigned>();	
}

bool Opts::needsHelp() const
{
	return _needsHelp;
}

std::string Opts::helpDesc() const
{
	return _helpDesc;
}

std::string Opts::inputFASTA() const
{
	return _inputFASTA;
}

unsigned Opts::tsneDim() const
{
	return _tsneDim;
}

unsigned Opts::tsnePerplexity() const
{
	return _tsnePerplexity;
}

float Opts::tsneTheta() const
{
	return _tsneTheta;
}

unsigned Opts::windowKmerLength() const
{
	return _windowKmerLength;
}

unsigned Opts::windowWidth() const
{
	return _windowWidth;
}

unsigned Opts::windowStep() const
{
	return _windowStep;
}
