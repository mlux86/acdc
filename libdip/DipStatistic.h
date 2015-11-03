#pragma once

#include <vector>

struct DipResult
{
	double dip;
	double p;
};

class DipStatistic
{

private:
	unsigned seed;

public:
	DipStatistic();
	DipStatistic(unsigned seed);
	~DipStatistic();

	DipResult calculate(std::vector<double> pdf, const unsigned numBootstraps = 100);

};
