#include "Logger.h"
#include "Clustering.h"
#include "DipStatistic.h"
#include "MLUtil.h"
#include "MatrixUtil.h"
#include "Kmeans.h"
#include "HierarchicalClustering.h"
#include "TarjansAlgorithm.h"

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

ClusteringResult Clustering::connComponents(const Eigen::MatrixXd & data, unsigned knnK, const std::vector<std::string> & contigs)
{
    Eigen::MatrixXd aff = MLUtil::knnAffinityMatrix(data, knnK, true);
    TarjansAlgorithm ta;
    auto res = ta.run(aff);
    Clustering::postprocess(res, contigs); 
    return res;
}

bool Clustering::isMultiModal(const Eigen::MatrixXd & data, double alpha, double splitThreshold)
{
    unsigned n = data.rows();

    Eigen::MatrixXd dist = MLUtil::pdist(data);

    DipStatistic ds;
    double splitPerc = 0.0;
    // for each cluster member, see if it is a split viewer
    // #pragma omp parallel for shared(splitPerc)
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

std::pair<unsigned, std::vector<ClusteringResult>> Clustering::estimateK(const Eigen::MatrixXd & data, unsigned maxK, const std::vector<std::string> & contigs)
{
    std::vector<ClusteringResult> results(maxK);

    results[0].numClusters = 1;
    std::fill(results[0].labels.begin(),results[0].labels.end(), 1);

    Eigen::MatrixXd z = HierarchicalClustering::linkage(data);

    double minDb = std::numeric_limits<double>::max();
    unsigned optK = 0;

    for (unsigned k = 2; k <= maxK; k++)
    {
        results[k-1].numClusters = k;
        results[k-1].labels = HierarchicalClustering::cluster(z, k);        
        Clustering::postprocess(results[k-1], contigs);

        double db = Clustering::daviesBouldin(data, results[k-1].labels);

        if (db < minDb)
        {
            minDb = db;
            optK = results[k-1].numClusters;
        }
    }

    return std::make_pair(optK, results);
}

void Clustering::postprocess(ClusteringResult & cr, const std::vector<std::string> & contigs)
{
    unsigned n = cr.labels.size();

    if (n != contigs.size())
    {
        throw std::runtime_error("Number of labels must match number of contigs.");
    }

    auto uniqueContigs = contigs;
    std::sort(uniqueContigs.begin(), uniqueContigs.end());
    auto it = std::unique(uniqueContigs.begin(), uniqueContigs.end());
    uniqueContigs.resize(std::distance(uniqueContigs.begin(), it));

    for (auto contig : uniqueContigs)
    {
        std::map<unsigned, unsigned> labelSizes;
        for (unsigned i = 0; i < n; i++)            
        {
            if (contigs.at(i) != contig)
            {
                continue;
            }
            unsigned lbl = cr.labels.at(i);
            if(labelSizes.count(lbl) == 0)
            {
                labelSizes[lbl] = 0;
            }
            labelSizes[lbl]++;
        }

        if (labelSizes.size() > 1) // multiple clusters for one contig
        {
            unsigned newLbl = 0;
            unsigned maxClustSize = 0;
            for (const auto it : labelSizes)
            {
                if (it.second > maxClustSize)
                {
                    newLbl = it.first;
                    maxClustSize = it.second;
                }
            }   
            for (unsigned i = 0; i < n; i++)            
            {
                if (contigs.at(i) == contig)
                {
                    cr.labels[i] = newLbl;
                }
            }
        }

    }

    // count clusters
    auto uniqueLabels = cr.labels;
    std::sort(uniqueLabels.begin(), uniqueLabels.end());
    auto it2 = std::unique(uniqueLabels.begin(), uniqueLabels.end());
    uniqueLabels.resize(std::distance(uniqueLabels.begin(), it2));
    cr.numClusters = uniqueLabels.size();

}