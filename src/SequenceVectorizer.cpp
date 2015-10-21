#include "SequenceVectorizer.h"
#include "Util.h"

#include <seqan/alignment_free.h> 
#include <seqan/sequence.h> 
#include <seqan/seq_io.h>
#include <string>
#include <vector>
#include <unordered_set>

SequenceVectorizer::SequenceVectorizer(const unsigned kmerLength_, const unsigned windowSize_, const unsigned windowStep_) : kmerLength(kmerLength_), windowSize(windowSize_), windowStep(windowStep_)
{
	unsigned idx = 0;
	for (auto kmer : Util::allKmers(kmerLength))
	{
		std::string revCompl = Util::reverseComplement(kmer);
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

SequenceVectorizer::~SequenceVectorizer() {}

unsigned SequenceVectorizer::getDim() const { return uniqueKmers.size(); }
std::set<std::string> SequenceVectorizer::getFeatures() const { return uniqueKmers; }
unsigned SequenceVectorizer::getKmerLength() const { return kmerLength; }
unsigned SequenceVectorizer::getWindowSize() const { return windowSize; }
unsigned SequenceVectorizer::getWindowStep() const { return windowStep; }

void SequenceVectorizer::setNormalize(const bool normalize_) { normalize = normalize_; }
bool SequenceVectorizer::getNormalize() const { return normalize; }

Eigen::MatrixXd SequenceVectorizer::vectorize(seqan::Dna5String & sequence) const
{
	unsigned n = ((int)(length(sequence) - (int)windowSize) / (int)windowStep) + 1;
	std::cout << n << std::endl;
	Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(n, getDim());

	for (unsigned i = 0; i < n; i++)
	{
		seqan::Infix<seqan::Dna5String>::Type window = infix(sequence, i*windowStep, i*windowStep+windowSize);

		seqan::String<unsigned> counts;
		countKmers(counts, window, kmerLength);

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
		Eigen::MatrixXd m = mat.rowwise().sum().asDiagonal().inverse();
		mat = m * mat; 
	}

	return mat;
}

std::pair< Eigen::MatrixXd, std::vector<std::string> > SequenceVectorizer::vectorize(const std::string & fasta) const
{
	seqan::SeqFileIn seqFileIn(fasta.c_str());
	seqan::StringSet<seqan::CharString> ids;
	seqan::StringSet<seqan::Dna5String> seqs;

	readRecords(ids, seqs, seqFileIn);

	unsigned n = length(ids);

	Eigen::MatrixXd result(0, getDim());
	std::vector<std::string> labels;

	for (unsigned i = 0; i < n; i++)
	{
		std::string id;
		move(id, ids[i]);
		auto & seq = seqs[i];

		auto mat = vectorize(seq);
		unsigned rows = mat.rows();

		// concatenate matrices TODO that's quite expensive
		auto tmp = result;
		result = Eigen::MatrixXd(tmp.rows() + mat.rows(), getDim());
		result << tmp, mat;
		
		for (unsigned j = 0; j < rows; j++)
		{
			labels.push_back(id);
		}
	}
	
	return std::make_pair(result, labels);
}