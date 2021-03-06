#include "catch.h"

#include "MatrixUtil.h"
#include "DipStatistic.h"

#include <algorithm>
#include <iostream>

TEST_CASE("Dip statistic", "[dip]")
{

	Eigen::MatrixXd unimodal = MatrixUtil::loadMatrix("../testdata/unimodal.txt", ' ');
	Eigen::MatrixXd bimodal = MatrixUtil::loadMatrix("../testdata/bimodal.txt", ' ');

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

	DipStatistic ds;

	DipResult resUni = ds.calculate(pdfUni, 100);
	REQUIRE(fabs(resUni.dip - 0.0203084) < 1e-7);
	REQUIRE(resUni.p > 0.01);

	DipResult resBi = ds.calculate(pdfBi, 100);
	REQUIRE(fabs(resBi.dip - 0.0395693) < 1e-7);
	REQUIRE(resUni.p > 0.01);

}
