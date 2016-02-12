#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include "MatrixUtil.h"

void MatrixUtil::matrixRemoveRow(Eigen::MatrixXd & matrix, unsigned int rowToRemove)
{
    unsigned int numRows = matrix.rows()-1;
    unsigned int numCols = matrix.cols();

    if( rowToRemove < numRows )
        matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.block(rowToRemove+1,0,numRows-rowToRemove,numCols);

    matrix.conservativeResize(numRows,numCols);
}

void MatrixUtil::matrixRemoveColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove)
{
    unsigned int numRows = matrix.rows();
    unsigned int numCols = matrix.cols()-1;

    if( colToRemove < numCols )
        matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.block(0,colToRemove+1,numRows,numCols-colToRemove);

    matrix.conservativeResize(numRows,numCols);
}

Eigen::MatrixXd MatrixUtil::loadMatrix(std::string filename, char delimiter)
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

void MatrixUtil::saveMatrix(const Eigen::MatrixXd & mat, std::string filename, char delimiter)
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

Json::Value MatrixUtil::matrixToJSON(const Eigen::MatrixXd & mat)
{
    unsigned n = mat.rows();
    unsigned dim = mat.cols();

    auto arr = Json::Value(Json::arrayValue);

    for(unsigned i = 0; i < n; i++)
    {
        auto entry = Json::Value(Json::arrayValue);
        for(unsigned j = 0; j < dim; j++)
        {
            entry.append(mat(i, j));
        }
        arr.append(entry);
    }

    return arr;
}
