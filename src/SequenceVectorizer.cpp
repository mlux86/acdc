#include "Logger.h"
#include "SequenceVectorizer.h"
#include "SequenceUtil.h"
#include "IOUtil.h"

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


SequenceVectorizer::SequenceVectorizer(const std::string & fasta, const Opts & opts)
{
	inputFasta = fasta;
	buildParams(opts);
	buildFeatureKmers();
}

SequenceVectorizer::~SequenceVectorizer() {}

void SequenceVectorizer::buildParams(const Opts & opts)
{
	minContigLength = opts.minContigLength();
	kmerLength = opts.windowKmerLength();
	windowWidth = opts.windowWidth();
	windowStep = opts.windowStep();
	targetNumPoints = opts.targetNumPoints();

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

Eigen::MatrixXd SequenceVectorizer::vectorize(seqan::Dna5String & sequence) const
{

	unsigned len = seqan::length(sequence);
	unsigned n = (unsigned) (((int)len - (int)windowWidth) / (int)windowStep) + 1;	

	if (n == 0) // include contigs smaller than window width
	{
		n++;
	}

	Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(n, getDim());

	for (unsigned i = 0; i < n; i++)
	{
		unsigned from = i*windowStep;
		unsigned to = std::min(len, i*windowStep+windowWidth);
		seqan::Infix<seqan::Dna5String>::Type window = infix(sequence, from, to);

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

	return mat;
}

SequenceVectorizationResult SequenceVectorizer::vectorize() const
{
	unsigned n = ids.size();

	SequenceVectorizationResult svr;

	Eigen::MatrixXd data(0, getDim());

	for (unsigned i = 0; i < n; i++)
	{
		std::string id = ids.at(i);
		seqan::Dna5String seq = sequences.at(i);

		try 
		{
			Eigen::MatrixXd mat = vectorize(seq);
			unsigned rows = mat.rows();

			// concatenate matrices TODO that's quite expensive
			auto tmp = data;
			data = Eigen::MatrixXd(tmp.rows() + mat.rows(), getDim());
			data << tmp, mat;
			
			for (unsigned j = 0; j < rows; j++)
			{
				svr.contigs.push_back(id);				
			}
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
