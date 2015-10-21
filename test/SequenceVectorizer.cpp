#include "catch.h"

#include "SequenceVectorizer.h"

#include <iostream>

TEST_CASE("Correct features", "[SequenceVectorizer]")
{

	SequenceVectorizer sv2(2, 1000, 500);
	REQUIRE(sv2.getFeatures().size() == 10);
	SequenceVectorizer sv3(3, 1000, 500);
	REQUIRE(sv3.getFeatures().size() == 32);
	SequenceVectorizer sv4(4, 1000, 500);
	REQUIRE(sv4.getFeatures().size() == 136);

	std::vector<std::string> target2 = {
		"AA", "AC", "AG", "AT", "CA", "CC", "CG", "GA", "GC", "TA"
	};
	auto feat2 = sv2.getFeatures();
	REQUIRE(std::equal(feat2.begin(), feat2.end(), target2.begin()));

	std::vector<std::string> target3 = {
		"AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "ATA", "ATC", "ATG", "CAA", "CAC",
		"CAG", "CCA", "CCC", "CCG", "CGA", "CGC", "CTA", "CTC", "GAA", "GAC", "GCA", "GCC", "GGA", "GTA", "TAA", "TCA"
	};
	auto feat3 = sv3.getFeatures();
	REQUIRE(std::equal(feat3.begin(), feat3.end(), target3.begin()));

	std::vector<std::string> target4 = {
		"AAAA", "AAAC", "AAAG", "AAAT", "AACA", "AACC", "AACG", "AACT", "AAGA", "AAGC", "AAGG", "AAGT", "AATA", "AATC",
		"AATG", "AATT", "ACAA", "ACAC", "ACAG", "ACAT", "ACCA", "ACCC", "ACCG", "ACCT", "ACGA", "ACGC", "ACGG", "ACGT",
		"ACTA", "ACTC", "ACTG", "AGAA", "AGAC", "AGAG", "AGAT", "AGCA", "AGCC", "AGCG", "AGCT", "AGGA", "AGGC", "AGGG",
		"AGTA", "AGTC", "AGTG", "ATAA", "ATAC", "ATAG", "ATAT", "ATCA", "ATCC", "ATCG", "ATGA", "ATGC", "ATGG", "ATTA",
		"ATTC", "ATTG", "CAAA", "CAAC", "CAAG", "CACA", "CACC", "CACG", "CAGA", "CAGC", "CAGG", "CATA", "CATC", "CATG",
		"CCAA", "CCAC", "CCAG", "CCCA", "CCCC", "CCCG", "CCGA", "CCGC", "CCGG", "CCTA", "CCTC", "CGAA", "CGAC", "CGAG",
		"CGCA", "CGCC", "CGCG", "CGGA", "CGGC", "CGTA", "CGTC", "CTAA", "CTAC", "CTAG", "CTCA", "CTCC", "CTGA", "CTGC",
		"CTTA", "CTTC", "GAAA", "GAAC", "GACA", "GACC", "GAGA", "GAGC", "GATA", "GATC", "GCAA", "GCAC", "GCCA", "GCCC",
		"GCGA", "GCGC", "GCTA", "GGAA", "GGAC", "GGCA", "GGCC", "GGGA", "GGTA", "GTAA", "GTAC", "GTCA", "GTGA", "GTTA",
		"TAAA", "TACA", "TAGA", "TATA", "TCAA", "TCCA", "TCGA", "TGAA", "TGCA", "TTAA"		
	};
	auto feat4 = sv4.getFeatures();
	REQUIRE(std::equal(feat4.begin(), feat4.end(), target4.begin()));

}

TEST_CASE("Vectorize different sequences", "[SequenceVectorizer]")
{

	seqan::Dna5String seq1 = "CTGGA";
	seqan::Dna5String seq2 = "CGCGGCATTTACCTGACAGCGACGTAGTCT";
	seqan::Dna5String seq3 = "GTGCGAGCTGAACGGGCAATGCAGGAAGAGTTCTACCTGGAACTGAAAGAA";


	{	
		SequenceVectorizer sv(2, 5, 1); 
		auto mat = sv.vectorize(seq1);
		Eigen::MatrixXd target(1,10);
		target << 0, 0, 1, 0, 1, 1, 0, 1, 0, 0;
		REQUIRE(mat.isApprox(target));
	}
	{
		SequenceVectorizer sv(2, 8, 5); 	
		auto mat = sv.vectorize(seq2);
		Eigen::MatrixXd target(5, 10);
		target << 0, 0, 0, 1, 1, 1, 2, 0, 2, 0,
		          2, 1, 0, 1, 1, 1, 0, 0, 0, 1,
		          0, 2, 1, 0, 2, 1, 0, 1, 0, 0,
		          0, 2, 1, 0, 1, 0, 1, 1, 1, 0,
		          0, 3, 1, 0, 0, 0, 1, 1, 0, 1;
		REQUIRE(mat.isApprox(target));
	}
	{
		SequenceVectorizer sv(3, 15, 15); 	
		auto mat = sv.vectorize(seq2);
		Eigen::MatrixXd target(2, 32);
		target << 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 2, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0,
                  0, 0, 0, 0, 1, 0, 2, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 2, 0, 0, 0, 1, 0, 0;
		REQUIRE(mat.isApprox(target));
	}
	{
		SequenceVectorizer sv(4, 10, 5); 	
		auto mat = sv.vectorize(seq3);
		Eigen::MatrixXd target(9, 136);
		target << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
		          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		          0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0,
		          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
                  0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0,
                  1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
                  1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
                  0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;
		REQUIRE(mat.isApprox(target));
	}	

}

TEST_CASE("Normalization", "[SequenceVectorizer]")
{

	seqan::Dna5String seq = "CTGAACGG";

	{	
		SequenceVectorizer sv(2, 5, 2); 
		sv.setNormalize(true);		
		auto mat = sv.vectorize(seq);
		std::cout << mat << std::endl;
		Eigen::MatrixXd target(2,10);
		target << 0.25, 0, 0.25, 0, 0.25, 0, 0, 0.25, 0, 0, 
                  0.25, 0.25, 0, 0, 0, 0, 0.25, 0.25, 0, 0;
		REQUIRE(mat.isApprox(target));
	}

}