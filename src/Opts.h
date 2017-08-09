#pragma once

#include <vector>
#include <string>
#include <map>

// Program arguments data class
// See implementation method initializeOnce() for description of arguments
class Opts
{
private:
	std::string _execPath = "";
	std::string _sharePath = "";

	std::string _cliCall = "";

	bool _needsHelp = false;
	std::string _helpDesc = "";
	int _logLevel = 0;

	std::vector<std::string> _inputFASTAs;

	unsigned _tsneDim = 0;
	unsigned _tsnePcaDim = 0;
	unsigned _tsnePerplexity = 0;
	float _tsneTheta = 0.0;

	unsigned _minContigLength = 0;
	unsigned _windowKmerLength = 0;
	unsigned _windowWidth = 0;
	unsigned _windowStep = 0;

	unsigned _targetNumPoints = 0;

	unsigned _numThreads = 0;
	unsigned _numBootstraps = 0;
	double _bootstrapRatio = 0;

	std::string _outputDir = "";

	std::string _krakenDb = "";

	unsigned _aggressiveThreshold = 0;

	// Singleton object, keep constructors and copy operations private
	Opts();
	~Opts();
	Opts(Opts const &) = delete;
    void operator=(Opts const &) = delete;	
	static Opts & getInstance();

public:
	// Call this method with program arguments once to initialize the singleton object
	static void initializeOnce(int argc, char *argv[]);
	static std::map <std::string, std::string> parameters();

	static std::string execPath();
	static std::string sharePath();
	static std::string cliCall();
	static bool needsHelp();
	static unsigned logLevel();
	static std::string helpDesc();
	static std::vector<std::string> inputFASTAs();
	static unsigned tsneDim();
	static unsigned tsnePcaDim();
	static unsigned tsnePerplexity();
	static float tsneTheta();
	static unsigned minContigLength();
	static unsigned windowKmerLength();
	static unsigned windowWidth();
	static unsigned windowStep();
	static unsigned targetNumPoints();
	static unsigned numThreads();
	static unsigned numBootstraps();
	static double bootstrapRatio();
	static std::string outputDir();
	static std::string krakenDb();
	static unsigned aggressiveThreshold();	
};
