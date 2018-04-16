#include "Opts.h"
#include "OptsAccumulator.h"
#include "Logger.h"
#include "IOUtil.h"
#include "version.h"

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/program_options.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/join.hpp>
#include <sstream>
#include <thread>

Opts::Opts()
{
}

Opts::~Opts()
{
}

Opts & Opts::getInstance()
{
	static Opts instance;
	return instance;
}

void Opts::initializeOnce(int argc, char *argv[])
{
	Opts & opts = getInstance();

	boost::program_options::options_description description("Usage", 250);

	unsigned threads = std::thread::hardware_concurrency();
	if (threads == 0)
	{
		std::cerr << "Could not detect number of cores. Defaulting to one thread." << std::endl;
		threads = 1;
	}

	// find path of binary

	char result[4096];
	unsigned count = readlink( "/proc/self/exe", result, 4096);
	std::string selfExe(result, (count > 0) ? count : 0);

	// find execution and asset paths

    boost::filesystem::path path(selfExe);
    boost::filesystem::path fullPath = boost::filesystem::canonical(boost::filesystem::system_complete(path));
    opts._execPath = fullPath.parent_path().string();
    opts._sharePath = opts._execPath + "/../share/acdc";

    std::vector<std::string> params(argv+1, argv+argc);
    opts._cliCall = fullPath.string() + " " + boost::algorithm::join(params, " ");

    // build boost program arguments

	description.add_options()
	    ("help,h", "Display this help message")
	    ("version,V", "Display the software version")
	    ("verbose,v", accumulator<int>(&(opts._logLevel))->implicit_value(1), "Verbose output (use -vv for more or -vvv for maximum verbosity)")
	    ("quiet,q", "No output")
	    ("random-seed,R", boost::program_options::value<unsigned>(), "Random seed number for reproducible results") 
	    ("input-fasta,i", boost::program_options::value<std::string>(), "Input FASTA file")
	    ("input-list,I", boost::program_options::value<std::string>(), "File with a list of input FASTA files, one per line")
	    ("tsne-dimension,d", boost::program_options::value<unsigned>()->default_value(2), "T-SNE dimension")
	    ("tsne-pca-dimension,u", boost::program_options::value<unsigned>()->default_value(50), "T-SNE initial PCA dimension")
	    ("tsne-perplexity,p", boost::program_options::value<unsigned>(), "T-SNE perplexity (overrides automatic estimation)")
	    ("tsne-theta,t", boost::program_options::value<float>()->default_value(0.5), "T-SNE parameter 'theta' of the underlying Barnes-Hut approximation")
	    ("min-contig-length,m", boost::program_options::value<unsigned>()->default_value(0), "Minimal length of contigs to consider for analysis")
	    ("window-kmer-length,k", boost::program_options::value<unsigned>()->default_value(4), "Length of the k-mers in the sequence vectorizer window")
	    ("window-width,w", boost::program_options::value<unsigned>(), "Width of the sliding sequence vectorizer window (overrides automatic estimation using number of target points)")
	    ("window-step,s", boost::program_options::value<unsigned>(), "Step of the sliding sequence vectorizer window (overrides automatic estimation using number of target points)")
	    ("target-num-points,n", boost::program_options::value<unsigned>()->default_value(1000), "Approximate number of target points for estimating window parameters") 
	    ("num-bootstraps,b", boost::program_options::value<unsigned>()->default_value(10), "Number of bootstraps")
	    ("bootstrap-ratio,r", boost::program_options::value<double>()->default_value(0.75), "Bootstrap subsampling ratio")
	    ("num-threads,T", boost::program_options::value<unsigned>()->default_value(threads), "Number of threads for bootstrap analysis  (default: detect number of cores)")
	    ("output-dir,o", boost::program_options::value<std::string>()->default_value("./results"), "Result output directory")
	    ("kraken-db,K", boost::program_options::value<std::string>()->default_value(""), "Database to use for Kraken classification")
	    ("aggressive-threshold,a", boost::program_options::value<unsigned>()->default_value(5000), "Aggressive threshold: Treat clusters having a bp size below this threshold as outliers. (Default = 0 = aggressive mode disabled)")
	    ("taxonomy-file,x", boost::program_options::value<std::string>()->default_value(""), "File with external taxonomy information")
	    ("target-taxonomy,X", boost::program_options::value<std::string>()->default_value(""), "Target taxonomy identifier")
	    ;

	boost::program_options::variables_map vm;
	boost::program_options::parsed_options parsed = boost::program_options::command_line_parser(argc, argv).options(description).allow_unregistered().run();
	boost::program_options::store(parsed, vm);
	boost::program_options::notify(vm);

	// initialize members from variables map

	std::stringstream ss;
	ss << description;
	opts._helpDesc = ss.str();

	if(vm.count("quiet"))
	{
		opts._logLevel = -1;
	} else
	{
		opts._logLevel = std::min(opts._logLevel, 3);
	}

	opts._needsVersion = vm.count("version");
	opts._version = ACDC_VERSION;
	opts._needsHelp = vm.count("help");
	
	if(vm.count("random-seed"))
	{
		opts._randomSeed = vm["random-seed"].as<unsigned>();
	} else
	{
		opts._randomSeed = time(NULL);
	}
	std::srand(Opts::randomSeed());
	
	if (vm.count("input-list"))
	{
		std::string fastaListFile = vm["input-list"].as<std::string>();
		boost::algorithm::trim(fastaListFile);
		opts._inputFASTAs = IOUtil::fileLinesToVec(fastaListFile);
	} else if (vm.count("input-fasta"))
	{
		std::string singleFasta = vm["input-fasta"].as<std::string>();
		boost::algorithm::trim(singleFasta);
		opts._inputFASTAs.push_back(singleFasta);
	}
	opts._tsneDim = vm["tsne-dimension"].as<unsigned>();
	opts._tsnePcaDim = vm["tsne-pca-dimension"].as<unsigned>();
	if (vm.count("tsne-perplexity"))
	{
		opts._tsnePerplexity = vm["tsne-perplexity"].as<unsigned>();
	}
	opts._tsneTheta = vm["tsne-theta"].as<float>();
	opts._minContigLength = vm["min-contig-length"].as<unsigned>();
	opts._minContigLength = std::max(100u, opts._minContigLength);
	opts._windowKmerLength = vm["window-kmer-length"].as<unsigned>();
	if (vm.count("window-width"))
	{
		opts._windowWidth = vm["window-width"].as<unsigned>();
	}
	if (vm.count("window-step"))
	{
		opts._windowStep = vm["window-step"].as<unsigned>();
	}
	opts._targetNumPoints = vm["target-num-points"].as<unsigned>();
	opts._numThreads = vm["num-threads"].as<unsigned>();
	opts._numThreads = std::max(1u, opts._numThreads);
	opts._numBootstraps = vm["num-bootstraps"].as<unsigned>();
	opts._numBootstraps = std::max(1u, opts._numBootstraps);
	opts._bootstrapRatio = vm["bootstrap-ratio"].as<double>();
	opts._outputDir = vm["output-dir"].as<std::string>();
	opts._krakenDb = vm["kraken-db"].as<std::string>();
	opts._aggressiveThreshold = vm["aggressive-threshold"].as<unsigned>();
	opts._taxonomyFile = vm["taxonomy-file"].as<std::string>();
	opts._targetTaxonomy = vm["target-taxonomy"].as<std::string>();
}

std::map <std::string, std::string> Opts::parameters()
{
	std::map <std::string, std::string> params;
	params["execPath"] = Opts::getInstance()._execPath;
	params["sharePath"] = Opts::getInstance()._sharePath;
	params["cliCall"] = Opts::getInstance()._cliCall;
	params["randomSeed"] = std::to_string(Opts::getInstance()._randomSeed);
	params["needsHelp"] = std::to_string(Opts::getInstance()._needsHelp);
	params["needsVersion"] = std::to_string(Opts::getInstance()._needsVersion);
	params["logLevel"] = std::to_string(Opts::getInstance()._logLevel);
	params["tsneDim"] = std::to_string(Opts::getInstance()._tsneDim);
	params["tsnePcaDim"] = std::to_string(Opts::getInstance()._tsnePcaDim);
	params["tsnePerplexity"] = std::to_string(Opts::getInstance()._tsnePerplexity);
	params["tsneTheta"] = std::to_string(Opts::getInstance()._tsneTheta);
	params["minContigLength"] = std::to_string(Opts::getInstance()._minContigLength);
	params["windowKmerLength"] = std::to_string(Opts::getInstance()._windowKmerLength);
	params["windowWidth"] = std::to_string(Opts::getInstance()._windowWidth);
	params["windowStep"] = std::to_string(Opts::getInstance()._windowStep);
	params["targetNumPoints"] = std::to_string(Opts::getInstance()._targetNumPoints);
	params["numThreads"] = std::to_string(Opts::getInstance()._numThreads);
	params["numBootstraps"] = std::to_string(Opts::getInstance()._numBootstraps);
	params["bootstrapRatio"] = std::to_string(Opts::getInstance()._bootstrapRatio);
	params["outputDir"] = Opts::getInstance()._outputDir;
	params["krakenDb"] = Opts::getInstance()._krakenDb;
	params["aggressiveThreshold"] = std::to_string(Opts::getInstance()._aggressiveThreshold);
	params["taxonomyFile"] = Opts::getInstance()._taxonomyFile;
	params["targetTaxonomy"] = Opts::getInstance()._targetTaxonomy;
	return params;
}

std::string Opts::execPath()
{
	return Opts::getInstance()._execPath;
}

std::string Opts::sharePath()
{
	return Opts::getInstance()._sharePath;
}

std::string Opts::cliCall()
{
	return Opts::getInstance()._cliCall;
}

bool Opts::needsHelp()
{
	return Opts::getInstance()._needsHelp;
}

bool Opts::needsVersion()
{
	return Opts::getInstance()._needsVersion;
}

std::string Opts::version()
{
	return Opts::getInstance()._version;
}

unsigned Opts::logLevel()
{
	return Opts::getInstance()._logLevel;
}

std::string Opts::helpDesc()
{
	return Opts::getInstance()._helpDesc;
}

unsigned Opts::randomSeed()
{
	return Opts::getInstance()._randomSeed;
}

std::vector<std::string> Opts::inputFASTAs()
{
	return Opts::getInstance()._inputFASTAs;
}

unsigned Opts::tsneDim()
{
	return Opts::getInstance()._tsneDim;
}

unsigned Opts::tsnePcaDim()
{
	return Opts::getInstance()._tsnePcaDim;
}

unsigned Opts::tsnePerplexity()
{
	return Opts::getInstance()._tsnePerplexity;
}

float Opts::tsneTheta()
{
	return Opts::getInstance()._tsneTheta;
}

unsigned Opts::minContigLength()
{
	return Opts::getInstance()._minContigLength;
}

unsigned Opts::windowKmerLength()
{
	return Opts::getInstance()._windowKmerLength;
}

unsigned Opts::windowWidth()
{
	return Opts::getInstance()._windowWidth;
}

unsigned Opts::windowStep()
{
	return Opts::getInstance()._windowStep;
}

unsigned Opts::targetNumPoints()
{
	return Opts::getInstance()._targetNumPoints;
}

unsigned Opts::numThreads()
{
	return Opts::getInstance()._numThreads;
}

unsigned Opts::numBootstraps()
{
	return Opts::getInstance()._numBootstraps;
}

double Opts::bootstrapRatio()
{
	return Opts::getInstance()._bootstrapRatio;
}

std::string Opts::outputDir()
{
	return Opts::getInstance()._outputDir;
}

std::string Opts::krakenDb()
{
	return Opts::getInstance()._krakenDb;
}

unsigned Opts::aggressiveThreshold()
{
	return Opts::getInstance()._aggressiveThreshold;
}

std::string Opts::taxonomyFile()
{
	return Opts::getInstance()._taxonomyFile;
}

std::string Opts::targetTaxonomy()
{
	return Opts::getInstance()._targetTaxonomy;
}