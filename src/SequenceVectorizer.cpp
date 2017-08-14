#include "Logger.h"
#include "SequenceVectorizer.h"
#include "SequenceUtil.h"
#include "Opts.h"

#include <string>
#include <vector>
#include <unordered_set>
#include <math.h>

#include <seqan/alignment_free.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

SequenceVectorizer::SequenceVectorizer(const unsigned kmerLength_, const unsigned windowWidth_, const unsigned windowStep_) : kmerLength(kmerLength_), windowWidth(windowWidth_), windowStep(windowStep_)
{
	buildFeatureKmers();
}


SequenceVectorizer::SequenceVectorizer(const std::string & fasta)
{
	inputFasta = fasta;
	buildParams();
	buildFeatureKmers();
}

SequenceVectorizer::~SequenceVectorizer() {}

void SequenceVectorizer::buildParams()
{
	minContigLength = Opts::minContigLength();
	kmerLength = Opts::windowKmerLength();
	windowWidth = Opts::windowWidth();
	windowStep = Opts::windowStep();
	targetNumPoints = Opts::targetNumPoints();

	loadFasta();
	estimateWindowParams();
}

void SequenceVectorizer::loadFasta()
{
	seqan::SeqFileIn seqFileIn(inputFasta.c_str());
	seqan::StringSet<seqan::CharString> ids_;
	seqan::StringSet<seqan::String<seqan::Iupac> > seqs_;

	seqan::readRecords(ids_, seqs_, seqFileIn);

	unsigned n = seqan::length(ids_);

	ids.clear();
	sequences.clear();

	// Use only contigs which are larger than the minimum length

	for (unsigned i = 0; i < n; i++)
	{
		std::string id;
		move(id, ids_[i]);
		seqan::Dna5String seq = seqs_[i];
		unsigned len = seqan::length(seq);

		if (len < minContigLength)
		{
			continue;
		}

		ids.push_back(id);
		sequences.push_back(seq);
	}
}

void SequenceVectorizer::estimateWindowParams()
{
	if(windowWidth == 0)
	{
		unsigned long long numNucleotides = 0;
		for (const auto & seq : sequences)
		{
			numNucleotides += seqan::length(seq);
		}
        DLOG << "numNucleotides = " << numNucleotides << std::endl;

        if (numNucleotides <= targetNumPoints)
        {
            throw std::runtime_error("Input data is too small for requested number of data points.");
        }

		// Window step is estimated using the number of nucleotides of all usable contigs

        windowStep = (unsigned) ceil((double)numNucleotides / (double)targetNumPoints);
        windowWidth = 2 * windowStep;
        DLOG << "k=" << kmerLength << "   "
                << "windowWidth=" << windowWidth << "   "
                << "windowStep=" << windowStep << "\n";
	}

	if(windowStep == 0)
	{
		windowStep = windowWidth / 2;
	}
}

void SequenceVectorizer::buildFeatureKmers()
{
	unsigned idx = 0;
	for (auto kmer : SequenceUtil::allKmers(kmerLength))
	{
		std::string revCompl = SequenceUtil::reverseComplement(kmer);
		std::string & final = revCompl < kmer ? revCompl : kmer;

		if(kmerIndexes.find(final) != kmerIndexes.end())
		{
			kmerIndexes[kmer] = kmerIndexes[final];
		} else
		{
			kmerIndexes[kmer] = idx++;
		}

		uniqueKmers.insert(final);
	}
}

std::pair<Eigen::MatrixXd, std::vector<WindowRange>> SequenceVectorizer::vectorize(seqan::Dna5String & sequence) const
{

	unsigned len = seqan::length(sequence);
	unsigned n = (unsigned) (((int)len - (int)windowWidth) / (int)windowStep) + 1;

	if (n == 0) // include contigs smaller than window width
	{
		n++;
	}

	// Perform actual vectorization using sliding windows

	Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(n, getDim());
	std::vector<WindowRange> windows(n);

	for (unsigned i = 0; i < n; i++)
	{
		WindowRange wr;
		wr.from = i*windowStep;
		wr.to = std::min(len, i*windowStep+windowWidth);
		seqan::Infix<seqan::Dna5String>::Type window = infix(sequence, wr.from, wr.to);

		windows[i] = wr;

		seqan::String<unsigned> counts;
		seqan::countKmers(counts, window, kmerLength);

		auto kmerIter = kmerIndexes.cbegin();
		for(const auto & cnt : counts)
		{
			std::string kmer = kmerIter->first;
			unsigned idx = kmerIter->second;
			mat(i, idx) += cnt;
			kmerIter++;
		}
	}

	// Normalize if needed

	if (normalize)
	{
		Eigen::MatrixXd m = mat.rowwise().sum().asDiagonal();
		for (unsigned i = 0; i < n; i++)
		{
			if (abs(m(i, i)) < 1e-10)
			{
				m(i, i) = 1;
			}
		}
		mat = m.inverse() * mat;
	}

	return std::make_pair(mat, windows);
}

SequenceVectorizationResult SequenceVectorizer::vectorize() const
{
	unsigned n = ids.size();

	SequenceVectorizationResult svr;

	Eigen::MatrixXd data(0, getDim());

	// Vectorize each contig individually and build up full data matrix

	for (unsigned i = 0; i < n; i++)
	{
		std::string id = ids.at(i);
		seqan::Dna5String seq = sequences.at(i);

		try
		{
			std::pair<Eigen::MatrixXd, std::vector<WindowRange>> r = vectorize(seq);
			Eigen::MatrixXd & mat = r.first;
			unsigned rows = mat.rows();

			auto tmp = data;
			data = Eigen::MatrixXd(tmp.rows() + mat.rows(), getDim());
			data << tmp, mat;

			for (unsigned j = 0; j < rows; j++)
			{
				svr.contigs.push_back(id);
			}

			svr.windows.insert(svr.windows.end(), r.second.begin(), r.second.end());
		}
		catch(const std::exception& e)
		{
			DLOG << e.what() << '\n';
		}

		svr.contigSizes[id] = seqan::length(seq);

	}

	svr.data = data;

	return svr;
}

unsigned SequenceVectorizer::getDim() const { return uniqueKmers.size(); }
std::set<std::string> SequenceVectorizer::getFeatures() const { return uniqueKmers; }
unsigned SequenceVectorizer::getKmerLength() const { return kmerLength; }
unsigned SequenceVectorizer::getWindowWidth() const { return windowWidth; }
unsigned SequenceVectorizer::getWindowStep() const { return windowStep; }
void SequenceVectorizer::setNormalize(const bool normalize_) { normalize = normalize_; }
bool SequenceVectorizer::getNormalize() const { return normalize; }
