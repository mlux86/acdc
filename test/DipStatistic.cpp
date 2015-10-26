#include "catch.h"

#include "Util.h"
#include "DipStatistic.h"

#include <algorithm>
#include <iostream>

TEST_CASE("Dip statistic", "[dip]")
{

	Eigen::MatrixXd unimodal = Util::loadMatrix("../data/unimodal.txt", ' ');
	Eigen::MatrixXd bimodal = Util::loadMatrix("../data/bimodal.txt", ' ');

	std::vector<double> pdfUni;
	for (unsigned i = 0; i < unimodal.rows(); i++)
	{
		pdfUni.push_back(unimodal(i, 0));
	}
	std::random_shuffle(pdfUni.begin(), pdfUni.end());

	std::vector<double> pdfBi;
	for (unsigned i = 0; i < bimodal.rows(); i++)
	{
		pdfBi.push_back(bimodal(i, 0));
	}
	std::random_shuffle(pdfBi.begin(), pdfBi.end());

	DipStatistic ds(42);

	DipResult resUni = ds.calculate(pdfUni, 100);
	REQUIRE(abs(resUni.dip - 0.0203084) < 1e-7);
	REQUIRE(abs(resUni.p - 0.91) < 1e-2);

	DipResult resBi = ds.calculate(pdfBi, 100);
	REQUIRE(abs(resBi.dip - 0.0395693) < 1e-7);
	REQUIRE(resBi.p < 1e-10);

}
