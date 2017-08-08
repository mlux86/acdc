#include "Clustering.h"
#include "DipStatistic.h"
#include "MLUtil.h"
#include "ClusteringUtil.h"
#include "HierarchicalClustering.h"
#include "TarjansAlgorithm.h"
#include "Opts.h"

#include <algorithm>
#include <limits>
#include <vector>

Clustering::Clustering(const Eigen::MatrixXd & data_, const std::vector<std::string> & contigs_, const std::map<std::string, unsigned> & contigSizes_) : data(data_), contigs(contigs_), contigSizes(contigSizes_)
{
}

Clustering::~Clustering()
{
}

ClusteringResult Clustering::connComponents(unsigned knnK)
{
    Eigen::MatrixXd aff = MLUtil::knnAffinityMatrix(data, knnK, true);
    TarjansAlgorithm ta;
    auto res = ta.run(aff);
    postprocess(res); 
    return res;
}

std::pair<unsigned, std::vector<ClusteringResult>> Clustering::estimateK(unsigned maxK)
{
    std::vector<ClusteringResult> results(maxK);

    // label trivial clustering for k == 1
    results[0].numClusters = 1;
    results[0].labels.resize(data.rows());
    std::fill(results[0].labels.begin(),results[0].labels.end(), 1);

    // compute linkage only once
    Eigen::MatrixXd z = HierarchicalClustering::linkage(data);

    // find minimum Davies Bouldin index
    double minDb = std::numeric_limits<double>::max();
    unsigned optK = 0;

    for (unsigned k = 2; k <= maxK; k++)
    {
        results[k-1].numClusters = k;
        results[k-1].labels = HierarchicalClustering::cluster(z, k);        
        postprocess(results[k-1]);

        double db = ClusteringUtil::daviesBouldin(data, results[k-1].labels);

        if (db < minDb)
        {
            minDb = db;
            optK = results[k-1].numClusters;
        }
    }

    return std::make_pair(optK, results);
}

void Clustering::postprocess(ClusteringResult & cr)
{
    unsigned n = cr.labels.size();

    if (n != contigs.size())
    {
        throw std::runtime_error("Number of labels must match number of contigs.");
    }

    // first, assign points in contigs which occur in distinct clusters, assign the cluster with the most points

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

        if (labelSizes.size() > 1) // multiple clusters for one contig, re-assign label
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

    // count clusters because they might have changed
    auto uniqueLabels = cr.labels;
    std::sort(uniqueLabels.begin(), uniqueLabels.end());
    auto it2 = std::unique(uniqueLabels.begin(), uniqueLabels.end());
    uniqueLabels.resize(std::distance(uniqueLabels.begin(), it2));
    cr.numClusters = uniqueLabels.size();



    // second, throw out outlier clusters if aggressive mode is on

    if (Opts::aggressiveThreshold() > 0)
    {
        // for each label, add up sizes of contained contigs and see if they are below threshold
        std::map<unsigned, unsigned> clusterSizes; // cluster label => num nucleotides in there

        for (auto contig : uniqueContigs)
        {
            for (unsigned i = 0; i < n; i++)
            {
                if (contigs.at(i) == contig)
                {
                    clusterSizes[cr.labels.at(i)] += contigSizes.at(contig);
                    break;
                }
            }
        }

        // mark outlier clusters if the cluster size is below threshold
        for (auto & it : clusterSizes)
        {
            unsigned lbl = it.first;
            unsigned clusterSize = it.second;

            if (clusterSize < Opts::aggressiveThreshold())
            {
                cr.outlierClusters.push_back(lbl);
                cr.numClusters--;
            }
        }
    }

}