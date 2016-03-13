#include "HierarchicalClustering.h"
#include "fastclusterAcdc.h"
#include "MLUtil.h"

#include <eigen3/Eigen/Dense>

void HierarchicalClustering::clusterNum(const Eigen::MatrixXd & linkage, std::vector<unsigned> & labels, unsigned idx, unsigned clsnum)
{
    unsigned n = labels.size();

    unsigned c1 = linkage(idx, 0);
    unsigned c2 = linkage(idx, 1);

    if (c1 < n)
    {
        labels[c1] = clsnum;
    } else
    {
        HierarchicalClustering::clusterNum(linkage, labels, c1-n, clsnum);
    }

    if (c2 < n)
    {
        labels[c2] = clsnum;
    } else
    {
        HierarchicalClustering::clusterNum(linkage, labels, c2-n, clsnum);
    }
}

std::vector<unsigned> HierarchicalClustering::cluster(const Eigen::MatrixXd & linkage, unsigned maxK)
{
    unsigned n = linkage.rows()+1;
    std::vector<unsigned> labels(n);

    unsigned i; 
    unsigned clsnum = 1; 
    for (unsigned k = n-maxK; k < n-1; k++) 
    { 
        // left tree 
        i = linkage(k, 0);
        if (i < n) 
        { 
            labels[i] = clsnum; 
            clsnum++; 
        } else if (i < 2*n-maxK) 
        { 
            HierarchicalClustering::clusterNum(linkage, labels, i-n, clsnum); 
            clsnum++; 
        } 

        // right tree 
        i = linkage(k, 1);
        if (i < n) 
        { 
            labels[i] = clsnum; 
            clsnum++; 
        } else if (i < 2*n-maxK) 
        { 
            HierarchicalClustering::clusterNum(linkage, labels, i-n, clsnum); 
            clsnum++; 
        } 
    } 

    return labels; 
}

Eigen::MatrixXd HierarchicalClustering::linkage(const Eigen::MatrixXd & data)
{

    unsigned n = data.rows();
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> dist = MLUtil::condensedPdist(data);
    double * const z = new double[(n-1) * 4];
    fc_linkage(dist.data(), z, n, ClusterMethods::FC_METHOD_METR_WARD);

    Eigen::MatrixXd res(n-1, 4);
    for (unsigned i = 0; i < n-1; i++)
    {
        for (unsigned j = 0; j < 4; j++)  
        {
            res(i, j) = *(z + i*4+j);
        }
    }

    delete[] z;

    return res;
}