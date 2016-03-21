#pragma once

#include <Eigen/Dense>
#include "Opts.h"
#include "tsne.h"

class BarnesHutSNEAdapter
{

private:
	BarnesHutSNEAdapter();
	~BarnesHutSNEAdapter();

	static unsigned estimateTsnePerplexity(const Eigen::MatrixXd & mat);

	static double * loadData(const Eigen::MatrixXd & eigendata);
	static Eigen::MatrixXd saveData(double* data, int n, int d);

public:
	static Eigen::MatrixXd runBarnesHutSNE(const Eigen::MatrixXd & eigendata, const Opts & opts);
	
};
