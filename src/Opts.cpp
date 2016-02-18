#include "Opts.h"
#include "OptsAccumulator.h"
#include "Logger.h"
#include "IOUtil.h"

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <sstream>
#include <thread>

Opts::Opts(int argc, char *argv[])
{
	initialize(argc, argv);
}

Opts::~Opts()
{
}

void Opts::initialize(int argc, char *argv[])
{
	boost::program_options::options_description description("Usage", 250);

	unsigned threads = std::thread::hardware_concurrency();
	if (threads == 0)
	{
		std::cerr << "Could not detect number of cores. Defaulting to one thread." << std::endl;
		threads = 1;
	} 

	// find path of binary
    boost::filesystem::path fullPath(boost::filesystem::initial_path<boost::filesystem::path>());
    fullPath = boost::filesystem::system_complete(boost::filesystem::path(argv[0]));
    _execPath = boost::filesystem::canonical(boost::filesystem::complete(fullPath)).parent_path().string();
    _sharePath = _execPath + "/../share/acdc";

	description.add_options()
	    ("help,h", "Display this help message")
	    ("verbose,v", accumulator<int>(&_logLevel)->implicit_value(1), "Verbose output (use -vv for more or -vvv for maximum verbosity)")
	    ("quiet,q", "No output")
	    ("input-fasta,i", boost::program_options::value<std::string>(), "Input FASTA file")
	    ("input-list,I", boost::program_options::value<std::string>(), "File with a list of input FASTA files, one per line")
	    ("tsne-dimension,d", boost::program_options::value<unsigned>()->default_value(2), "T-SNE dimension")
	    ("tsne-pca-dimension,u", boost::program_options::value<unsigned>()->default_value(50), "T-SNE initial PCA dimension")
	    ("tsne-perplexity,p", boost::program_options::value<unsigned>(), "T-SNE perplexity (overrides automatic estimation)")
	    ("tsne-theta,t", boost::program_options::value<float>()->default_value(0.5), "T-SNE parameter 'theta' of the underlying Barnes-Hut approximation")
	    ("min-contig-length,m", boost::program_options::value<unsigned>()->default_value(2000), "Minimal length of contigs to consider for analysis")
	    ("window-kmer-length,k", boost::program_options::value<unsigned>()->default_value(4), "Length of the k-mers in the sequence vectorizer window")
	    ("window-width,w", boost::program_options::value<unsigned>(), "Width of the sliding sequence vectorizer window (overrides automatic estimation using number of target points)")
	    ("window-step,s", boost::program_options::value<unsigned>(), "Step of the sliding sequence vectorizer window (overrides automatic estimation using number of target points)")
	    ("target-num-points,n", boost::program_options::value<unsigned>()->default_value(500), "Approximate number of target points for estimating window parameters")
	    ("num-bootstraps,b", boost::program_options::value<unsigned>()->default_value(10), "Number of bootstraps")
	    ("bootstrap-ratio,r", boost::program_options::value<double>()->default_value(0.75), "Bootstrap subsampling ratio")
	    ("num-threads,T", boost::program_options::value<unsigned>()->default_value(threads), "Number of threads for bootstrap analysis  (default: detect number of cores)")
	    ("output-dir,o", boost::program_options::value<std::string>()->default_value("./results"), "Result output directory")
	    ("kraken-db,K", boost::program_options::value<std::string>()->default_value(""), "Database to use for Kraken classification")
	    ;
	
	boost::program_options::variables_map vm;
	boost::program_options::parsed_options parsed = boost::program_options::command_line_parser(argc, argv).options(description).allow_unregistered().run(); 
	boost::program_options::store(parsed, vm);
	boost::program_options::notify(vm);

	std::stringstream ss;
	ss << description;
	_helpDesc = ss.str();

	if(vm.count("quiet")) 
	{
		_logLevel = -1;
	} else
	{
		_logLevel = std::min(_logLevel, 3);
	}

	_needsHelp = vm.count("help");
	if (vm.count("input-list"))
	{
		std::string fastaListFile = vm["input-list"].as<std::string>();
		boost::algorithm::trim(fastaListFile);	
		_inputFASTAs = IOUtil::fileLinesToVec(fastaListFile);
	} else if (vm.count("input-fasta"))
	{
		std::string singleFasta = vm["input-fasta"].as<std::string>();
		boost::algorithm::trim(singleFasta);
		_inputFASTAs.push_back(singleFasta);
	}
	_tsneDim = vm["tsne-dimension"].as<unsigned>();
	_tsnePcaDim = vm["tsne-pca-dimension"].as<unsigned>();
	if (vm.count("tsne-perplexity"))
	{
		_tsnePerplexity = vm["tsne-perplexity"].as<unsigned>();
	}
	_tsneTheta = vm["tsne-theta"].as<float>();
	_minContigLength = vm["min-contig-length"].as<unsigned>();
	_minContigLength = std::max(100u, _minContigLength);
	_windowKmerLength = vm["window-kmer-length"].as<unsigned>();
	if (vm.count("window-width"))
	{
		_windowWidth = vm["window-width"].as<unsigned>();
	}
	if (vm.count("window-step"))
	{
		_windowStep = vm["window-step"].as<unsigned>();
	}	
	_targetNumPoints = vm["target-num-points"].as<unsigned>();
	_numThreads = vm["num-threads"].as<unsigned>();
	_numThreads = std::max(1u, _numThreads);
	_numBootstraps = vm["num-bootstraps"].as<unsigned>();
	_numBootstraps = std::max(1u, _numBootstraps);
	_bootstrapRatio = vm["bootstrap-ratio"].as<double>();
	_outputDir = vm["output-dir"].as<std::string>();
	_krakenDb = vm["kraken-db"].as<std::string>();
}

std::string Opts::execPath() const
{
	return _execPath;
}

std::string Opts::sharePath() const
{
	return _sharePath;
}

bool Opts::needsHelp() const
{
	return _needsHelp;
}

unsigned Opts::logLevel() const
{
	return _logLevel;
}

std::string Opts::helpDesc() const
{
	return _helpDesc;
}

std::vector<std::string> Opts::inputFASTAs() const
{
	return _inputFASTAs;
}

unsigned Opts::tsneDim() const
{
	return _tsneDim;
}

unsigned Opts::tsnePcaDim() const
{
	return _tsnePcaDim;
}

unsigned Opts::tsnePerplexity() const
{
	return _tsnePerplexity;
}

float Opts::tsneTheta() const
{
	return _tsneTheta;
}

unsigned Opts::minContigLength() const
{
	return _minContigLength;
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

unsigned Opts::targetNumPoints() const
{
	return _targetNumPoints;
}

unsigned Opts::numThreads() const
{
	return _numThreads;
}

unsigned Opts::numBootstraps() const
{
	return _numBootstraps;
}

double Opts::bootstrapRatio() const
{
	return _bootstrapRatio;
}

std::string Opts::outputDir() const
{
	return _outputDir;
}

std::string Opts::krakenDb() const
{
	return _krakenDb;
}