#include "catch.h"

#include "MatrixUtil.h"
#include "MLUtil.h"
#include "TarjansAlgorithm.h"

#include <algorithm>
#include <iostream>

TEST_CASE("Tarjans algorithm, two clusters", "[tarjan]")
{
	Eigen::MatrixXd clustData = MatrixUtil::loadMatrix("../data/twoclusters.txt", ' ');
	Eigen::MatrixXd clustDataLabels = MatrixUtil::loadMatrix("../data/twoclusters.labels.txt", ' ');

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

TEST_CASE("Tarjans algorithm, two clusters & one mini cluster", "[tarjan]")
{
	Eigen::MatrixXd clustData = MatrixUtil::loadMatrix("../data/threeclusters.txt", ' ');
	Eigen::MatrixXd clustDataLabels = MatrixUtil::loadMatrix("../data/threeclusters.labels.txt", ' ');

	Eigen::MatrixXd affinities = MLUtil::knnAffinityMatrix(clustData, 7, false);
	TarjansAlgorithm ta;
	auto ccRes = ta.run(affinities);

	Eigen::MatrixXd labels(ccRes.labels.size(), 1);
	for (unsigned i = 0; i < ccRes.labels.size(); ++i)
	{
		labels(i, 0) = ccRes.labels.at(i);
	}

	REQUIRE(ccRes.numClusters == 3);
	REQUIRE(labels.isApprox(clustDataLabels));


	affinities = MLUtil::knnAffinityMatrix(clustData, 8, false); // one more, but cluster is only 8 points large, not enough for 8 neighbors
	ccRes = ta.run(affinities);

	labels = Eigen::MatrixXd(ccRes.labels.size(), 1);
	for (unsigned i = 0; i < ccRes.labels.size(); ++i)
	{
		labels(i, 0) = ccRes.labels.at(i);
	}

	REQUIRE(ccRes.numClusters == 2);
	REQUIRE(!labels.isApprox(clustDataLabels));

}