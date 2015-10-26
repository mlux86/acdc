#include "DipStatistic.h"

extern "C"
{
	#include "dip.h"
}

#include <algorithm>

#include <stdlib.h>
#include <time.h>

DipStatistic::DipStatistic()
{
	DipStatistic(time(NULL));
}

DipStatistic::DipStatistic(unsigned seed_)
{
	seed = seed_;
}

DipStatistic::~DipStatistic()
{
}

DipResult DipStatistic::calculate(std::vector<double> pdf, const unsigned numBootstraps)
{
	DipResult result;

	unsigned n = pdf.size();

	std::sort(pdf.begin(), pdf.end());
	{
		int discard_lohi[4];
		int discard_ifault;
		int * discard_gcm = new int[n];
		int * discard_lcm = new int[n];
		int * discard_mn = new int[n];
		int * discard_mj = new int[n];
		diptest(pdf.data(), n, &(result.dip), discard_lohi, &discard_ifault, discard_gcm, discard_lcm, discard_mn, discard_mj);
		delete[] discard_gcm;
		delete[] discard_lcm;
		delete[] discard_mn;
		delete[] discard_mj;
	}


	unsigned bootDipsSmaller = 0;
	for (unsigned i = 0; i < numBootstraps; i++)
	{

		srand(seed + i);

		// get dip from generated uniform [0,1] reference pdf set
		std::vector<double> bootPdf(n);
		for (unsigned j = 0; j < n; j++)
		{
			bootPdf[j] = ((double) rand() / (RAND_MAX)) + 1;
		}

		std::sort(bootPdf.begin(), bootPdf.end());
		double bootDip;
		{
			int discard_lohi[4];
			int discard_ifault;
			int * discard_gcm = new int[n];
			int * discard_lcm = new int[n];
			int * discard_mn = new int[n];
			int * discard_mj = new int[n];
			diptest(bootPdf.data(), n, &bootDip, discard_lohi, &discard_ifault, discard_gcm, discard_lcm, discard_mn, discard_mj);
			delete[] discard_gcm;
			delete[] discard_lcm;
			delete[] discard_mn;
			delete[] discard_mj;			
		}

		if (result.dip < bootDip)
		{
			bootDipsSmaller++;
		}
	}

	// calculate p-value
	result.p = (double)bootDipsSmaller / (double)numBootstraps;
	
	return result;
}