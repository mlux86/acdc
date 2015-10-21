#include <iostream>

#include "SequenceVectorizer.h"
#include "BarnesHutSNEBridge.h"
#include "Util.h"

int main()
{

	SequenceVectorizer sv(4, 4000, 2000);

	auto dat = sv.vectorize("/home/mlux/downloads/contigs.fa");

	Eigen::MatrixXd m = dat.first.rowwise().maxCoeff().asDiagonal().inverse();
	dat.first = m * dat.first; 	

	auto reduced = BarnesHutSNEBridge::runBarnesHutSNE(dat.first, 2, .5, 50);

	Util::saveMatrix(reduced, "/tmp/reduced", ' ');

	return 0;
}