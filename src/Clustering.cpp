#include "Logger.h"
#include "Clustering.h"
#include "DipStatistic.h"
#include "Util.h"
#include "Kmeans.h"

#include <Eigen/Eigenvalues>
#include <math.h>  
#include <numeric>  
#include <algorithm>
#include <chrono>
#include <random>
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
        if (degreeMat(i, i) == 0)
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
        degreeMat(i, i) = 1.0 / sqrt(degreeMat(i, i));
    }

    Eigen::MatrixXd laplacian = degreeMat * adjacencies * degreeMat;
    Eigen::MatrixXd eye = Eigen::MatrixXd::Identity(n, n);
    laplacian += 0.000001 * eye;

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes(laplacian, false);

    Eigen::VectorXd eigvals = saes.eigenvalues();

    ClusteringResult res;

    res.numClusters = 0;
    for (unsigned i = 0; i < eigvals.size(); i++)
    {
        if (eigvals(i) > 1 - 0.0001 && eigvals(i) < 1 + 0.0001)
        {
            res.numClusters++;
        }
    }

    return res;
}

ClusteringResult Clustering::dipMeans(const Eigen::MatrixXd& data, double alpha, double splitThreshold)
{
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 generator(seed);

    unsigned n = data.rows();
    unsigned dim = data.cols();
    
    Eigen::MatrixXd dist = Util::pdist(data);

    unsigned k = 1;
    Eigen::VectorXd labels = Eigen::VectorXd::Zero(n);

    Eigen::MatrixXd means = data.colwise().sum() / n;

    while (true)
    {

        std::vector<double> scores;

        for (unsigned j = 0; j < k; j++) // calculate scores for clusters c_j
        {
            // find indexes of cluster j
            std::vector<unsigned> clusterIdx;
            for (unsigned i = 0; i < n; i++)
            {
                if (labels(i) == j)
                {
                    clusterIdx.push_back(i);
                }
            }

            unsigned numSplitViewers = 0;
            double splitScore = 0.0;
            DipStatistic ds;
            // for each cluster member, see if it is a split viewer
            for (auto idx : clusterIdx)
            {
                std::vector<double> viewerDists;
                for (auto idx2 : clusterIdx)
                {
                    if (idx != idx2)
                    {
                        viewerDists.push_back(dist(idx, idx2));
                    }
                }
                auto dipRes = ds.calculate(viewerDists, 1000);
                if (dipRes.p <= alpha)
                {
                    numSplitViewers++;
                    splitScore += dipRes.dip;
                }
            }
            splitScore /= numSplitViewers;
            double splitPerc = (double)numSplitViewers / (double)clusterIdx.size(); 

            if (splitPerc >= splitThreshold)
            {
                scores.push_back(splitScore);
            } else
            {
                scores.push_back(0.0);
            }
        }

        // to split or not to split...

        auto maxScoreIter = std::max_element(scores.begin(), scores.end());
        if (*maxScoreIter > 0.0)
        {
            unsigned splitLabel = std::distance(scores.begin(), maxScoreIter);

            DLOG << "split " << splitLabel << " with score=" << scores[splitLabel] << "\n";

            // assemble data matrix with cluster members, split means

            Eigen::MatrixXd members(0, dim);
            for (unsigned i = 0; i < n; i++)
            {
                if (labels(i) == splitLabel)
                {
                    members.conservativeResize(members.rows()+1, dim);
                    members.row(members.rows()-1) = data.row(i);
                }
            }

            Eigen::VectorXd oldMean = means.row(splitLabel);
            std::uniform_int_distribution<unsigned> distn(0, members.rows() - 1);

            Eigen::MatrixXd locallyOptimizedMeans;

            double minMse = std::numeric_limits<double>::max();
            for (unsigned i = 0; i < 5; i++)
            {
                unsigned randIdx = distn(generator);
                Eigen::MatrixXd splitMeans = Eigen::MatrixXd(2, dim);
                splitMeans.row(0) = members.row(randIdx);
                splitMeans.row(1) = oldMean - (members.row(randIdx).transpose() - oldMean);

                Kmeans km(2);
                km.initMeans(splitMeans);
                auto tmp = km.iteration(members);
                if (tmp.second < minMse)
                {
                    minMse = tmp.second;
                    locallyOptimizedMeans = km.getMeans();
                }
            }

            // refresh means
            k++;
            Util::matrixRemoveRow(means, splitLabel);
            Eigen::MatrixXd newMeans(k, dim);
            newMeans << means, 
                        locallyOptimizedMeans;


            // finally ... run k-means for refinement on full data set            
            Kmeans km(k);
            km.initMeans(newMeans);
            ClusteringResult res = km.run(data);
            labels = res.labels;
            means = km.getMeans();

        } else
        {
            break;
        }
    }

    ClusteringResult res;
    res.numClusters = k;
    res.labels = labels;
    return res;       
}
