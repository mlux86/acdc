#include "Kmeans.h"

#include <math.h>  
#include <numeric> 
#include <ctime>
#include <cstdlib>
#include <algorithm>
#include <chrono>
#include <random>
#include <limits>

Kmeans::Kmeans(const unsigned k_) : k(k_)
{
}

Kmeans::~Kmeans()
{
}

void Kmeans::initSample()
{
	init = "sample";
}

void Kmeans::initPlusPlus()
{
	init = "plusplus";
}

void Kmeans::initMeans(const Eigen::MatrixXd & initMeans)
{
	init = "means";
	means = initMeans;
}

Eigen::MatrixXd Kmeans::getMeans()
{
	return means;
}

void Kmeans::setNumBootstraps(const unsigned n)
{
	numBootstraps = n;
}

ClusteringResult Kmeans::run(const Eigen::MatrixXd & data)
{

    ClusteringResult bestRes;
    double minMse = std::numeric_limits<double>::max();

    for (unsigned i = 0; i < numBootstraps; i++)
    {
        auto tmp = iteration(data);
        if (tmp.second < minMse)
        {
            minMse = tmp.second;
            bestRes = tmp.first;
        }
    }

    return bestRes;
}

std::pair<ClusteringResult, double> Kmeans::iteration(const Eigen::MatrixXd & data)
{
    unsigned n = data.rows();
    unsigned dim = data.cols();

    ClusteringResult res;
    res.labels = Eigen::VectorXd::Zero(n);
    res.numClusters = k;

    if (init == "sample")
    {
    	sampleInit(data);
    } else if (init == "plusplus")
    {
    	plusPlusInit(data);
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

void Kmeans::sampleInit(const Eigen::MatrixXd & data)
{
    unsigned n = data.rows();
    unsigned dim = data.cols();

    std::srand(unsigned(std::time(0)));

	std::vector<unsigned> indexes(n);
	std::iota(indexes.begin(), indexes.end(), 0);
	std::random_shuffle(indexes.begin(), indexes.end());

	means = Eigen::MatrixXd::Zero(k, dim);
	for (unsigned i = 0; i < k; i++)
	{
		means.row(i) = data.row(indexes[i]);
	}
}

void Kmeans::plusPlusInit(const Eigen::MatrixXd & data)
{
    unsigned n = data.rows();
    unsigned dim = data.cols();
    
    std::vector<unsigned> meanIndexes;

    // choose first mean uniformly at random
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 generator(seed);
    std::uniform_int_distribution<unsigned> distn(0, n - 1);
    meanIndexes.push_back(distn(generator));

    // iteratively add more means
    for (unsigned j = 0; j < k - 1; j++)
    {

        // calculate probability distribution of data based on distances
        std::vector<double> probs(n, 2); // default probabilities of 2
        for (auto idx : meanIndexes) // set already chosen indexes to zero probability
        {
            probs[idx] = 0;
        }
        
        // calculate shortest distances to closes already chosen mean
        std::vector<double> dists;
        double sumDist = 0;
        for (unsigned i = 0; i < n; i++)
        {
            double minDist = std::numeric_limits<double>::max();
            for (auto idx : meanIndexes)
            {
                double d = (data.row(i) - data.row(idx)).squaredNorm();
                if (d < minDist)
                {
                    minDist = d;
                }
            }
            dists.push_back(minDist);
            sumDist += minDist;
        }
        
        // populate probabilities		
        for (unsigned i = 0; i < n; i++)
        {
            if (probs[i] > 1.0) // uninitialized
            {
                probs[i] = dists[i] / sumDist;
            }
        }
        
        // calculate cumulative probabilities
        std::vector<double> cumProbs(n);
        cumProbs[0] = probs[0];
        for (unsigned i = 1; i < n; i++)
        {
            cumProbs[i] = cumProbs[i - 1] + probs[i];
        }
        
        // select next mean index based on cumulative probabilities
        std::uniform_real_distribution<double> unif(0.0, 1.0);
        double r = unif(generator);
        for (unsigned i = 0; i < n; i++)
        {
            if (r <= cumProbs[i])
            {
                meanIndexes.push_back(i);
                break;
            }
        }

    }

    means = Eigen::MatrixXd::Zero(k, dim);
    for (unsigned j = 0; j < k; j++)
    {
        means.row(j) = data.row(meanIndexes[j]);
    }
}