#include "Logger.h"
#include "Clustering.h"
#include "DipStatistic.h"
#include "MLUtil.h"
#include "MatrixUtil.h"
#include "Kmeans.h"
#include "HierarchicalClustering.h"

#include <math.h>
#include <numeric>
#include <algorithm>
#include <chrono>
#include <random>
#include <limits>
#include <vector>

Clustering::Clustering()
{
}

Clustering::~Clustering()
{
}

bool Clustering::isMultiModal(const Eigen::MatrixXd & data, double alpha, double splitThreshold)
{
    unsigned n = data.rows();

    Eigen::MatrixXd dist = MLUtil::pdist(data);

    DipStatistic ds;
    double splitPerc = 0.0;
    // for each cluster member, see if it is a split viewer
    #pragma omp parallel for shared(splitPerc)
    for (unsigned i = 0; i < n; i++)
    {
        if (splitPerc < splitThreshold)
        {
            std::vector<double> viewerDists;
            for (unsigned j = 0; j < n; j++)
            {
                if (i != j)
                {
                    viewerDists.push_back(dist(i, j));
                }
            }
            auto dipRes = ds.calculate(viewerDists, 1000);
            if (dipRes.p <= alpha)
            {
                splitPerc += 1.0 / (double)n;
            }
        }
    }

    return splitPerc >= splitThreshold;
}

double Clustering::daviesBouldin(const Eigen::MatrixXd & data, const std::vector<unsigned> & labels)
{    
    unsigned n = data.rows();
    unsigned dim = data.cols();

    auto uniqueLabels = labels;
    auto it = std::unique(uniqueLabels.begin(), uniqueLabels.end());
    uniqueLabels.resize(std::distance(uniqueLabels.begin(), it));
    std::sort(uniqueLabels.begin(), uniqueLabels.end());

    unsigned numClusters = uniqueLabels.size();
    std::vector<unsigned> clusterSizes(numClusters);

    // compute class wise means
    Eigen::MatrixXd mu(numClusters, dim);
    for (unsigned j = 0; j < numClusters; j++)
    {
        unsigned lbl = uniqueLabels.at(j);        

        // filter elemens with labels == lbl
        Eigen::MatrixXd members(0, dim);
        for (unsigned i = 0; i < n; i++)
        {
            if (lbl == labels.at(i))
            {
                members.conservativeResize(members.rows()+1, dim);
                members.row(members.rows()-1) = data.row(i);
            }
        }
        clusterSizes[j] = members.rows();
        mu.row(j) = members.colwise().sum() / members.rows();
    }

    // compute data to data and intra-cluster distances
    Eigen::MatrixXd tmp(n+numClusters, dim);
    tmp << data, mu;
    Eigen::MatrixXd distances = MLUtil::pdist(tmp);
    Eigen::MatrixXd dc = distances.topRightCorner(n, numClusters); // data to mean distances
    Eigen::MatrixXd cc = distances.bottomRightCorner(numClusters, numClusters); // intra-cluster    

    double db = 0;

    for (unsigned i = 0; i < numClusters; i++)
    {

        unsigned lbl = uniqueLabels.at(i);

        double di = 0; // mean intra cluster distance for i
        for (unsigned j = 0; j < n; j++)
        {
            if (lbl == labels.at(j))
            {
                di += dc(j, i);
            }
        }
        di /= clusterSizes.at(i);

        double s = 0;

        for (unsigned j = 0; j < numClusters; j++)
        {
            unsigned lbl2 = uniqueLabels.at(j);
            if (lbl == lbl2)
            {
                continue;
            }

            double dj = 0; // mean intra cluster distance for i
            for (unsigned k = 0; k < n; k++)
            {
                if (lbl2 == labels.at(k))
                {
                    dj += dc(k, j);
                }
            }
            dj /= clusterSizes.at(j);

            double s_ = (di + dj) / cc(i, j);

            if (s_ > s)
            {
                s = s_;
            }
        }

        db += s;

    }

    db /= numClusters;

    return db;
}

ClusteringResult Clustering::estimateK(const Eigen::MatrixXd & data, unsigned maxK)
{
    ClusteringResult res;

    Eigen::MatrixXd z = HierarchicalClustering::linkage(data);

    double minDb = std::numeric_limits<double>::max();
    for (unsigned k = 2; k < maxK; k++)
    {
        std::vector<unsigned> labels = HierarchicalClustering::cluster(z, k);
        double db = Clustering::daviesBouldin(data, labels);

        if (db < minDb)
        {
            minDb = db;
            res.labels = labels;
            res.numClusters = k;
        }
    }

    return res;
}
