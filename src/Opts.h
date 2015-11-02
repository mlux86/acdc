#ifndef __Opts__
#define __Opts__

#include <string>

class Opts
{
private:
	bool _needsHelp = false;
	std::string _helpDesc = "";
	int _logLevel = 0;

	std::string _inputFASTA = "";

	unsigned _tsneDim = 0;
	unsigned _tsnePcaDim = 0;
	unsigned _tsnePerplexity = 0;
	float _tsneTheta = 0.0;

	unsigned _windowKmerLength = 0;
	unsigned _windowWidth = 0;
	unsigned _windowStep = 0;

	unsigned _targetNumPoints = 0;

	void initialize(int argc, char const *argv[]);


public:
	Opts(int argc, char const *argv[]);
	~Opts();

	bool needsHelp() const;
	unsigned logLevel() const;
	std::string helpDesc() const;
	std::string inputFASTA() const;
	unsigned tsneDim() const;
	unsigned tsnePcaDim() const;
	unsigned tsnePerplexity() const;
	float tsneTheta() const;
	unsigned windowKmerLength() const;
	unsigned windowWidth() const;
	unsigned windowStep() const;
	unsigned targetNumPoints() const;
	
};

#endif 
