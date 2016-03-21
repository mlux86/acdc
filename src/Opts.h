#pragma once

#include <vector>
#include <string>

class Opts
{
private:
	std::string _execPath = "";
	std::string _sharePath = "";

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

	Opts();
	~Opts();
	Opts(Opts const &) = delete;
    void operator=(Opts const &) = delete;	
	static Opts & getInstance();

public:
	static void initializeOnce(int argc, char *argv[]);

	static std::string execPath();
	static std::string sharePath();
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
