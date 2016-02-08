#include "Util.h"

#include <iostream>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <math.h>  
#include "nanoflann.hpp"
#include <numeric>
#include <unordered_map>

#include <seqan/alignment_free.h> 
#include <seqan/sequence.h> 
#include <seqan/seq_io.h>

Util::Util()
{
}

Util::~Util()
{
}

void Util::allPermsRepetition(std::vector<std::string> & perms, const std::vector<char> & alphabet, std::string elem, const unsigned length)
{
    if (length == 0)
    {
        perms.push_back(elem);
        return;
    }

    for (auto c : alphabet)
    {
        Util::allPermsRepetition(perms, alphabet, elem + c, length - 1);
    }
}

std::vector<std::string> Util::allKmers(const unsigned kmerLength) 
{
    std::vector<char> alphabet = {'A', 'C', 'G', 'T'};
    std::vector<std::string> perms;
    Util::allPermsRepetition(perms, alphabet, "", kmerLength);
    return perms;
}

std::string Util::reverseComplement(const std::string & seq) 
{
    std::string revCompl(seq);
    std::reverse(revCompl.begin(), revCompl.end());    

    for (auto & c : revCompl)
    {
        switch(c)
        {
            case 'A':
                c = 'T';
                break;
            case 'C':
                c = 'G';
                break;
            case 'G':
                c = 'C';
                break;
            case 'T':
                c = 'A';
                break;                                                
        }
    }

    return revCompl;
}

void Util::matrixRemoveRow(Eigen::MatrixXd & matrix, unsigned int rowToRemove)
{
    unsigned int numRows = matrix.rows()-1;
    unsigned int numCols = matrix.cols();

    if( rowToRemove < numRows )
        matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.block(rowToRemove+1,0,numRows-rowToRemove,numCols);

    matrix.conservativeResize(numRows,numCols);
}

void Util::matrixRemoveColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove)
{
    unsigned int numRows = matrix.rows();
    unsigned int numCols = matrix.cols()-1;

    if( colToRemove < numCols )
        matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.block(0,colToRemove+1,numRows,numCols-colToRemove);

    matrix.conservativeResize(numRows,numCols);
}

Eigen::MatrixXd Util::loadMatrix(std::string filename, char delimiter)
{
    std::ifstream ifs(filename, std::ifstream::in);
    if (!ifs)
    {
        throw std::runtime_error("Cannot read file '" + filename + "'!");
    }

    std::string delim(1, delimiter);

    std::string line;
    std::getline(ifs, line);

    std::vector<std::string> strs;
    boost::split(strs, line, boost::is_any_of(delim));    

    auto rows = boost::lexical_cast<unsigned>(strs.front());
    auto cols = boost::lexical_cast<unsigned>(strs.back());

    Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(rows, cols);

    for (unsigned i = 0; i < rows; i++)
    {
        std::getline(ifs, line);

        strs.clear();
        boost::split(strs, line, boost::is_any_of(delim));

        for (unsigned j = 0; j < cols; j++)
        {
            auto val = boost::lexical_cast<float>(strs[j]);
            mat(i, j) = val;
        }
    }

    ifs.close();

    return mat;
}

void Util::saveMatrix(const Eigen::MatrixXd & mat, std::string filename, char delimiter)
{
    std::ofstream ofs(filename, std::ifstream::out);

    if (!ofs)
    {
        throw std::runtime_error("Cannot open file '" + filename + "' for writing!");
    }

    unsigned rows = mat.rows();
    unsigned cols = mat.cols();

    ofs << rows << delimiter << cols << '\n';

    for (unsigned i = 0; i < rows; i++)
    {
        for (unsigned j = 0; j < cols; j++)
        {
            ofs << mat(i, j);
            if (j < cols - 1)
            {
                ofs << delimiter;
            }
        }
        ofs << '\n';
    }

    ofs.close();
}

unsigned Util::estimateTsnePerplexity(const Eigen::MatrixXd & mat)
{
    return (unsigned) (pow(log(mat.rows()), 2)); 
}

unsigned long long Util::getFileSizeBytes(const std::string & filename)
{
    std::ifstream in(filename, std::ifstream::ate | std::ifstream::binary);
    return in.tellg(); 
}

Eigen::MatrixXd Util::knnAffinityMatrix(const Eigen::MatrixXd & data, const unsigned k_, bool mutual)
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
            if (!mutual)
            {
                affinities(idx, i) = sqrt(dists[j]); // for normal k-nn
            }
        }
    }

    if (mutual)
    {
        for (unsigned i = 0; i < n; i++)
        {
            for (unsigned j = 0; j <= i; j++)
            {
                auto m = std::min(affinities(i, j), affinities(j, i));
                affinities(i, j) = m;
                affinities(j, i) = m;
            }
        }    
    }

    return affinities;
}

Eigen::MatrixXd Util::pca(const Eigen::MatrixXd & data_, const unsigned ndims)
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

Eigen::MatrixXd Util::pdist(const Eigen::MatrixXd & data)
{
    unsigned n = data.rows();

    Eigen::MatrixXd distances(n, n);

    for (unsigned i = 0; i < n; i++)
    {
        distances.row(i) = (data.rowwise() - data.row(i)).rowwise().norm();
    }

    return distances;
}

std::vector< std::vector<unsigned> > Util::stratifiedSubsamplingIndices(const unsigned n, const unsigned k, const double ratio)
{
    std::vector<unsigned> indices(n);
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), std::mt19937{std::random_device{}()});

    unsigned m = (unsigned) ((double)n * ratio);
    unsigned step = (n - m) / k;

    std::vector< std::vector<unsigned> > result;

    for (unsigned i = 0; i < k; i++)
    {
        unsigned from = i*step;
        unsigned to = i*step + m;
        if (i == k-1)
        {
            to = std::max(to, n);
        }

        std::vector<unsigned> e(indices.begin()+from, indices.begin()+to);
        result.push_back(e); 
    }

    return result;
}

Json::Value Util::clusteringToJson(const Eigen::MatrixXd & mat, const std::vector<unsigned> & labels, const std::vector<std::string> & tooltips)
{
    unsigned n = mat.rows();

    std::vector<std::string> colors = {"#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999"};
    auto jsonMat = Json::Value(Json::arrayValue);
    
    for(unsigned i = 0; i < n; i++)
    {
        auto jsonMatPt = Json::Value();
        jsonMatPt["x"] = mat(i, 0);
        jsonMatPt["y"] = mat(i, 1);
        unsigned colorIdx = labels.at(i) % colors.size();
        jsonMatPt["color"] = colors[colorIdx];
        jsonMatPt["toolTipContent"] = tooltips[i];
        jsonMatPt["label"] = labels.at(i);
        jsonMat.append(jsonMatPt);
    }

    return jsonMat;
}

std::vector<std::string> Util::fileLinesToVec(const std::string & filename)
{
    std::ifstream ifs(filename);
    std::string str;
    std::vector<std::string> result;
    while (std::getline(ifs, str))
    {
        if(!str.empty())
        {
            result.push_back(str);
        }
    }     
    ifs.close();
    return result;
}

std::tuple<Eigen::MatrixXd, Eigen::MatrixXd, Eigen::VectorXd, Eigen::MatrixXd, Eigen::MatrixXd> Util::canonicalCorrelation(const Eigen::MatrixXd & x_, const Eigen::MatrixXd & y_)
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

Eigen::MatrixXd Util::alignBootstrap(const Eigen::MatrixXd & reference, const Eigen::MatrixXd & bootstrap, const std::vector<unsigned> & bootstrapIndexes)
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

    auto ccaRes = canonicalCorrelation(x, y);

    Eigen::MatrixXd a = std::get<0>(ccaRes);
    Eigen::MatrixXd v = std::get<4>(ccaRes);

    Eigen::MatrixXd result = v * a.inverse();

    return result;
}

std::vector<unsigned> Util::alignBootstrapLabels(const std::vector<unsigned> & referenceLabels, const std::vector<unsigned> & bootstrapLabels, const std::vector<unsigned> & bootstrapIndexes)
{
    unsigned n = bootstrapLabels.size();

    std::unordered_map<unsigned, unsigned> mp;

    for (unsigned i = 0; i < n; i++)
    {
        unsigned lbl = bootstrapLabels.at(i);
        if (mp.count(lbl) == 0)
        {            
            mp[lbl] = referenceLabels.at(bootstrapIndexes[i]);
        }
    }

    std::vector<unsigned> result(n, 0);

    for (unsigned i = 0; i < n; i++)
    {
        result[i] = mp.at(bootstrapLabels.at(i));   
    }

    return result;
}

std::unique_ptr<std::string> Util::filterFasta(const std::string & fasta, const std::vector<std::string> contigs)
{

    std::stringstream ss;
    
    seqan::SeqFileIn seqFileIn(fasta.c_str());
    seqan::StringSet<seqan::CharString> ids;
    seqan::StringSet<seqan::String<seqan::Iupac> > seqs;

    seqan::readRecords(ids, seqs, seqFileIn);

    unsigned n = seqan::length(ids);

    for (unsigned i = 0; i < n; i++)
    {
        std::string id;
        move(id, ids[i]);
        if (std::find(contigs.begin(), contigs.end(), id) != contigs.end())
        {
            std::string seq;
            move(seq, seqs[i]);
            ss << '>' << id << '\n' << seq << '\n';
        }
    }

    std::unique_ptr<std::string> result(new std::string(ss.str()));
    return result;
}