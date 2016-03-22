#include "catch.h"

#include "MatrixUtil.h"
#include "MLUtil.h"
#include "TarjansAlgorithm.h"
#include "Clustering.h"

#include <algorithm>
#include <iostream>

TEST_CASE("Tarjans algorithm, two clusters", "[tarjan]")
{
	Eigen::MatrixXd clustData = MatrixUtil::loadMatrix("../testdata/twoclusters.txt", '\t');
	Eigen::MatrixXd clustDataLabels = MatrixUtil::loadMatrix("../testdata/twoclusters.labels.txt", '\t');

	Eigen::MatrixXd affinities = MLUtil::knnAffinityMatrix(clustData, 7, false);
	TarjansAlgorithm ta;
	auto ccRes = ta.run(affinities);

	Eigen::MatrixXd labels(ccRes.labels.size(), 1);
	for (unsigned i = 0; i < ccRes.labels.size(); ++i)
	{
		labels(i, 0) = ccRes.labels.at(i);
	}

	REQUIRE(ccRes.numClusters == 2);
	REQUIRE(labels.isApprox(clustDataLabels));

}

TEST_CASE("Tarjans algorithm, chess board", "[tarjan]")
{
	Eigen::MatrixXd clustData = MatrixUtil::loadMatrix("../testdata/chessboard.txt", '\t');
	Eigen::MatrixXd clustDataLabels = MatrixUtil::loadMatrix("../testdata/chessboard.labels.txt", '\t');

	Eigen::MatrixXd affinities = MLUtil::knnAffinityMatrix(clustData, 49, true);
	TarjansAlgorithm ta;
	auto ccRes = ta.run(affinities);

	Eigen::MatrixXd labels(ccRes.labels.size(), 1);
	for (unsigned i = 0; i < ccRes.labels.size(); ++i)
	{
		labels(i, 0) = ccRes.labels.at(i);
	}

	REQUIRE(ccRes.numClusters == 64);
	REQUIRE(labels.isApprox(clustDataLabels));


	affinities = MLUtil::knnAffinityMatrix(clustData, 50, true); // all clusters contain only 50 points, sho if neighborhood is larger, only one big cluster should "survive"
	ccRes = ta.run(affinities);

	labels = Eigen::MatrixXd(ccRes.labels.size(), 1);
	for (unsigned i = 0; i < ccRes.labels.size(); ++i)
	{
		labels(i, 0) = ccRes.labels.at(i);
	}

	REQUIRE(ccRes.numClusters == 1);
	REQUIRE(!labels.isApprox(clustDataLabels));

}