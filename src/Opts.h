#pragma once

#include <vector>
#include <string>

#include "Util.h"

class Opts
{
private:
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

	unsigned _port = 0;

	std::string _krakenScript = "";

	void initialize(int argc, char const *argv[]);


public:
	Opts(int argc, char const *argv[]);
	~Opts();

	bool needsHelp() const;
	unsigned logLevel() const;
	std::string helpDesc() const;
	std::vector<std::string> inputFASTAs() const;
	unsigned tsneDim() const;
	unsigned tsnePcaDim() const;
	unsigned tsnePerplexity() const;
	float tsneTheta() const;
	unsigned minContigLength() const;
	unsigned windowKmerLength() const;
	unsigned windowWidth() const;
	unsigned windowStep() const;
	unsigned targetNumPoints() const;
	unsigned numThreads() const;
	unsigned numBootstraps() const;
	double bootstrapRatio() const;
	unsigned port() const;
	std::string krakenScript() const;
	
};
