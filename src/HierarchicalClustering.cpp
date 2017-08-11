#include "HierarchicalClustering.h"
#include "fastclusterAcdc.h"
#include "MLUtil.h"

#include <Eigen/Dense>

HierarchicalClustering::HierarchicalClustering()
{

}

HierarchicalClustering::~HierarchicalClustering()
{

}

void HierarchicalClustering::clusterNum(const Eigen::MatrixXd & linkage, std::vector<unsigned> & labels, unsigned idx, unsigned clsnum)
{
    unsigned n = labels.size();

    unsigned c1 = linkage(idx, 0);
    unsigned c2 = linkage(idx, 1);

    if (c1 < n) // left leaf, mark label
    {
        labels[c1] = clsnum; 
    } else
    {
        HierarchicalClustering::clusterNum(linkage, labels, c1-n, clsnum); // recurse left tree
    }

    if (c2 < n) // right leaf, mark label
    {
        labels[c2] = clsnum;
    } else
    {
        HierarchicalClustering::clusterNum(linkage, labels, c2-n, clsnum); // recurse right tree
    }
}

void HierarchicalClustering::initialize(const Eigen::MatrixXd & data)
{
    linkageMat = HierarchicalClustering::linkage(data);
}

std::vector<unsigned> HierarchicalClustering::cluster(unsigned maxK)
{
    unsigned n = linkageMat.rows()+1;
    std::vector<unsigned> labels(n);

    unsigned i; 
    unsigned clsnum = 1; 
    for (unsigned k = n-maxK; k < n-1; k++) 
    { 
        // left tree 
        i = linkageMat(k, 0);
        if (i < n) // is leaf, mark label
        { 
            labels[i] = clsnum; 
            clsnum++; 
        } else if (i < 2*n-maxK) 
        { 
            HierarchicalClustering::clusterNum(linkageMat, labels, i-n, clsnum); // recurse left tree
            clsnum++; 
        } 

        // right tree 
        i = linkageMat(k, 1);
        if (i < n) // is leaf, mark label
        { 
            labels[i] = clsnum; 
            clsnum++; 
        } else if (i < 2*n-maxK) 
        { 
            HierarchicalClustering::clusterNum(linkageMat, labels, i-n, clsnum); // recurse right tree
            clsnum++; 
        } 
    } 

    return labels; 
}

Eigen::MatrixXd HierarchicalClustering::linkage(const Eigen::MatrixXd & data)
{
    // compute pairwise distances
    unsigned n = data.rows();
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> dist = MLUtil::condensedPdist(data);
    
    // call fastcluster linkage routine
    double * const z = new double[(n-1) * 4];
    fc_linkage(dist.data(), z, n, ClusterMethods::FC_METHOD_METR_WARD);

    // assemble Eigen object from returned linkage
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

std::string HierarchicalClustering::name()
{
    return "HierarchicalClustering";
}

std::map<std::string, std::string> HierarchicalClustering::parameters()
{
    std::map<std::string, std::string> params = {{"distance", "euclidean"}, {"linkage", "ward"}};
    return params;
}
