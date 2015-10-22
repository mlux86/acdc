#include "easylogging++.h"

#include "Clustering.h"
#include "Util.h"

#include <Eigen/Eigenvalues>
#include <math.h>  

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
