#ifndef OPTS_H
#define OPTS_H

#include <string>

class Opts
{
private:
	bool _needsHelp = false;
	std::string _helpDesc = "";

	std::string _inputFASTA = "";

	unsigned _tsneDim;
	unsigned _tsnePerplexity;
	float _tsneTheta;

	unsigned _windowKmerLength;
	unsigned _windowWidth;
	unsigned _windowStep;

	void initialize(int argc, char const *argv[]);


public:
	Opts(int argc, char const *argv[]);
	~Opts();

	bool needsHelp() const;
	std::string helpDesc() const;
	std::string inputFASTA() const;
	unsigned tsneDim() const;
	unsigned tsnePerplexity() const;
	float tsneTheta() const;
	unsigned windowKmerLength() const;
	unsigned windowWidth() const;
	unsigned windowStep() const;
	
};

#endif 
