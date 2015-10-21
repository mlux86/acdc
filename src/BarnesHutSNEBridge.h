#ifndef __BarnesHutSNEBridge__
#define __BarnesHutSNEBridge__

#include <eigen3/Eigen/Dense>
#include "tsne.h"

class BarnesHutSNEBridge
{

private:
	BarnesHutSNEBridge();
	~BarnesHutSNEBridge();

	static double * loadData(const Eigen::MatrixXd & eigendata);
	static Eigen::MatrixXd saveData(double* data, int n, int d);

public:
	static Eigen::MatrixXd runBarnesHutSNE(const Eigen::MatrixXd & eigendata, const unsigned targetDim, const double theta, const double perplexity);
	
};

#endif