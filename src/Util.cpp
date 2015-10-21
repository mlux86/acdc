#include "Util.h"

#include <algorithm>
#include <fstream>
#include <cstdlib>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
// #include <boost/algorithm/string/classification.hpp>

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




