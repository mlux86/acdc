#include "easylogging++.h"

#include "Clustering.h"
#include "Util.h"

#include <Eigen/Eigenvalues>
#include <math.h>  
#include <numeric>  
#include <algorithm>
#include <limits>

Clustering::Clustering()
{
}

Clustering::~Clustering()
{
}

ClusteringResult Clustering::spectralClustering(Eigen::MatrixXd & adjacencies)
{
	unsigned n = adjacencies.rows();
	if (n != adjacencies.cols())
	{
		throw std::runtime_error("Adjacency matrix must be quadratic!");
	}

	Eigen::MatrixXd degreeMat = adjacencies.rowwise().sum().asDiagonal();

	// check for zero degree points and construct a smaller matrix if necessary
	std::vector<unsigned> zeroDegreeIndexes;
	for (unsigned i = 0; i < n; i++)
	{
		if (degreeMat(i,i) == 0)
		{
			zeroDegreeIndexes.push_back(i);
		}
	}
	// reverse iterate because largest row to remove is at end
	// otherwise, removing smaller rows will shift indexes
	for (auto iter = zeroDegreeIndexes.rbegin(); iter != zeroDegreeIndexes.rend(); ++iter) 
	{
		auto idx = *iter;
		Util::matrixRemoveRow(adjacencies, idx);
		Util::matrixRemoveRow(degreeMat, idx);
		Util::matrixRemoveColumn(adjacencies, idx);
		Util::matrixRemoveColumn(degreeMat, idx);
		n--;
	}	 

	// square root of a diagonal matrix
	for (unsigned i = 0; i < n; i++) 
	{
		degreeMat(i, i) = 1.0 / sqrt(degreeMat(i,i));
	}

	Eigen::MatrixXd laplacian = degreeMat * adjacencies * degreeMat;
	Eigen::MatrixXd eye = Eigen::MatrixXd::Identity(n, n);
	laplacian += 0.000001 * eye;

	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes(laplacian, false);

	Eigen::VectorXd eigvals = saes.eigenvalues();

	ClusteringResult res;

	res.numClusters = 0;
	for (unsigned i = 0; i< eigvals.size(); i++)
	{
		if (eigvals(i) > 1-0.0001 && eigvals(i) < 1+0.0001)
		{
			res.numClusters++;
		}
	}

	return res;
}

ClusteringResult Clustering::kMeans(const Eigen::MatrixXd & data, unsigned k, unsigned numBootstraps)
{

	ClusteringResult bestRes;
	double minMse = std::numeric_limits<double>::max();

	for (unsigned i = 0; i < numBootstraps; i++)
	{
		auto tmp = Clustering::kMeansIter(data, k);
		VLOG(2) << tmp.second;
		if (tmp.second < minMse)
		{
			minMse = tmp.second;
			bestRes = tmp.first;
		}
	}

	return bestRes;
}

std::pair<ClusteringResult, double> Clustering::kMeansIter(const Eigen::MatrixXd & data, unsigned k)
{
	unsigned n = data.rows();
	unsigned dim = data.cols();

	ClusteringResult res;
	res.labels = Eigen::MatrixXd(n, 1);
	res.numClusters = k;

	// initialize random means
	Eigen::MatrixXd means = Eigen::MatrixXd::Zero(k, dim);
	std::vector<unsigned> tmpIdx(n);
	std::iota(tmpIdx.begin(), tmpIdx.end(), 0);
	std::random_shuffle(tmpIdx.begin(), tmpIdx.end());
	for (unsigned j = 0; j < k; j++)
	{
		means.row(j) = data.row(tmpIdx[j]);
	}

	Eigen::MatrixXd oldMeans;
	unsigned numEpsilonChanges = 0;
	double mse = 0;

	// start iterations
	while (numEpsilonChanges < 10)
	{

		// assignment
		mse = 0;
		for (unsigned i = 0; i < n; i++)
		{
			double minNorm = std::numeric_limits<double>::max();
			for (unsigned j = 0; j < k; j++)
			{
				double norm = (data.row(i) - means.row(j)).squaredNorm();
				if (norm < minNorm)
				{
					minNorm = norm;
					res.labels(i) = j;
				}
				mse += norm;
			}
		}
		mse /= n;

		// recalculate means
		oldMeans = means;
		means = Eigen::MatrixXd::Zero(k, dim);
		for (unsigned j = 0; j < k; j++)
		{
			unsigned nk = 0;
			for (unsigned i = 0; i < n; i++)
			{
				if (res.labels(i) == j)
				{
					means.row(j) += data.row(i);
					nk++;
				}
			}
			means.row(j) /= nk;
		}

		// check for convergence
		double meansDiff = (oldMeans - means).cwiseAbs().sum();		
		if (meansDiff < 1e-10)
		{
			numEpsilonChanges++;
		} else
		{
			numEpsilonChanges = 0;
		}

	}
	
	return std::make_pair(res, mse);
}
