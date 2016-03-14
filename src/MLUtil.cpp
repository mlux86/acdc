#include <unordered_map>
#include <set>

#include "MLUtil.h"
#include "nanoflann.hpp"

Eigen::MatrixXd MLUtil::knnAffinityMatrix(const Eigen::MatrixXd & data, const unsigned k_, bool mutual)
{
    using RowMajorMatrixXd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    using EigenKdTree = nanoflann::KDTreeEigenMatrixAdaptor<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>;

    RowMajorMatrixXd r = data; //convert to row-major for easy access to data points

    const unsigned k = k_ + 1; // identical point is always included

    unsigned n = r.rows();
    unsigned dim = r.cols();

    if (k > n)
    {
        std::stringstream ss;
        ss << "Cannot find " << k << " neighbors in " << n << " data points!";
        throw std::runtime_error(ss.str());
    }

    EigenKdTree tree(dim, r);
    tree.index->buildIndex();

    std::vector<size_t> indexes(k);
    std::vector<double> dists(k);
    nanoflann::KNNResultSet<double> resultSet(k);
    resultSet.init(&indexes[0], &dists[0]);

    Eigen::MatrixXd affinities = Eigen::MatrixXd::Zero(n, n);

    for (unsigned i = 0; i < n; i++)
    {
        resultSet.init(&indexes[0], &dists[0]);
        tree.index->findNeighbors(resultSet, r.data() + i*dim, nanoflann::SearchParams());
        
        for (unsigned j = 0; j < k; j++)
        {
            unsigned idx = indexes[j];
            affinities(i, idx) = sqrt(dists[j]);
        }
    }

    if (mutual)
    {
        for (unsigned i = 0; i < n; i++)
        {
            for (unsigned j = 0; j < n; j++)
            {
                auto m = std::min(affinities(i, j), affinities(j, i));
                affinities(i, j) = m;
                affinities(j, i) = m;
            }
        }    
    }

    return affinities;
}

Eigen::MatrixXd MLUtil::pca(const Eigen::MatrixXd & data_, const unsigned ndims)
{
    Eigen::MatrixXd data = data_;

    unsigned n = data.rows();

    // subtract mean
    Eigen::VectorXd mean = data.colwise().sum() / n;
    data.rowwise() -= mean.transpose();

    // covariance matrix
    Eigen::MatrixXd cov = data.transpose() * data / (double) n;

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> saes(cov);
    Eigen::MatrixXd eigvecs = saes.eigenvectors();

    Eigen::MatrixXd proj = eigvecs.rightCols(ndims);

    return data * proj;
}

Eigen::MatrixXd MLUtil::pdist(const Eigen::MatrixXd & data)
{
    unsigned n = data.rows();

    Eigen::MatrixXd distances(n, n);

    for (unsigned i = 0; i < n; i++)
    {
        distances.row(i) = (data.rowwise() - data.row(i)).rowwise().norm();
    }

    return distances;
}

Eigen::VectorXd MLUtil::condensedPdist(const Eigen::MatrixXd & data)
{
    unsigned n = data.rows();

    Eigen::VectorXd distances(n*(n-1)/2);

    unsigned c = 0;
    for (unsigned i = 0; i < n; i++)
    {
        for (unsigned j = i+1; j < n; j++)
        {
            distances(c++) = (data.row(i) - data.row(j)).norm();
        }
    }

    return distances;
}

std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::VectorXd, Eigen::MatrixXd, Eigen::MatrixXd> MLUtil::canonicalCorrelation(const Eigen::MatrixXd & x_, const Eigen::MatrixXd & y_)
{
    unsigned n = x_.rows();
    unsigned dimX = x_.cols();
    unsigned dimY = y_.cols();

    if (n != y_.rows())
    {
        throw std::runtime_error("Number of rows of x must equal number of rows of y.");
    }

    Eigen::MatrixXd x = x_;
    Eigen::MatrixXd y = y_;

    // subtract means
    Eigen::VectorXd meanX = x.colwise().sum() / n;
    x.rowwise() -= meanX.transpose();
    Eigen::VectorXd meanY = y.colwise().sum() / n;
    y.rowwise() -= meanY.transpose();

    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qrX(x);
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qrY(y);
    unsigned rankX = qrX.rank();
    unsigned rankY = qrY.rank();
    Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic> permX = qrX.colsPermutation();
    Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic> permY = qrY.colsPermutation();    
    Eigen::MatrixXd qX_ = qrX.householderQ();
    Eigen::MatrixXd qY_ = qrY.householderQ();
    Eigen::MatrixXd rX_ = qrX.matrixR();
    Eigen::MatrixXd rY_ = qrY.matrixR();

    if (rankX == 0 || rankY == 0)
    {
        throw std::runtime_error("Rank of input matrices must not be zero.");
    }

    Eigen::MatrixXd qX = qX_.leftCols(rankX);
    Eigen::MatrixXd rX = rX_.topLeftCorner(rankX, rankX).triangularView<Eigen::Upper>();

    Eigen::MatrixXd qY = qY_.leftCols(rankY);
    Eigen::MatrixXd rY = rY_.topLeftCorner(rankY, rankY).triangularView<Eigen::Upper>();

    unsigned d = std::min(dimX, dimY);

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(qX.transpose() * qY, Eigen::ComputeThinU | Eigen::ComputeThinV);

    Eigen::VectorXd r = svd.singularValues();

    Eigen::MatrixXd a_ = rX.inverse() * svd.matrixU().leftCols(d) * sqrt((double)n - 1.0);
    Eigen::MatrixXd b_ = rY.inverse() * svd.matrixV().leftCols(d) * sqrt((double)n - 1.0);

    Eigen::MatrixXd a = permX * a_;
    Eigen::MatrixXd b = permY * b_;

    Eigen::MatrixXd u = x * a;
    Eigen::MatrixXd v = y * b;

    return std::make_tuple(a, b, r, u, v);
}

Eigen::MatrixXd MLUtil::alignBootstrap(const Eigen::MatrixXd & reference, const Eigen::MatrixXd & bootstrap, const std::vector<unsigned> & bootstrapIndexes)
{
    unsigned n = bootstrap.rows();
    unsigned dim = bootstrap.cols();

    Eigen::MatrixXd x = Eigen::MatrixXd::Zero(n, dim);
    Eigen::MatrixXd y = Eigen::MatrixXd::Zero(n, dim);

    for (unsigned i = 0; i < n; ++i)
    {
        x.row(i) = reference.row(bootstrapIndexes[i]);
        y.row(i) = bootstrap.row(i);
    }    

    auto ccaRes = MLUtil::canonicalCorrelation(x, y);

    Eigen::MatrixXd a = std::get<0>(ccaRes);
    Eigen::MatrixXd v = std::get<4>(ccaRes);

    Eigen::MatrixXd result = v * a.inverse();

    return result;
}

std::vector<unsigned> MLUtil::alignBootstrapLabels(const std::vector<unsigned> & referenceLabels, const std::vector<unsigned> & bootstrapLabels, const std::vector<unsigned> & bootstrapIndexes)
{
    unsigned n = bootstrapLabels.size();

    std::unordered_map<unsigned, unsigned> mp;
    std::set<unsigned> assigned;

    unsigned k = *std::max_element(referenceLabels.begin(), referenceLabels.end()) + 1;

    for (unsigned i = 0; i < n; i++)
    {
        unsigned lblFrom = bootstrapLabels.at(i);
        unsigned lblTo = referenceLabels.at(bootstrapIndexes.at(i));
        if (mp.count(lblFrom) == 0) 
        {   
            if (assigned.find(lblTo) == assigned.end())
            {
                mp[lblFrom] = lblTo;
                assigned.insert(lblTo);
            } else
            {
                mp[lblFrom] = k++;
            }
        }
    }

    std::vector<unsigned> result(n, 0);

    for (unsigned i = 0; i < n; i++)
    {
        result[i] = mp.at(bootstrapLabels.at(i));   
    }

    return result;
}
