#ifndef __BarnesHutSNEBridge__
#define __BarnesHutSNEBridge__

#include <eigen3/Eigen/Dense>
#include "Opts.h"
#include "tsne.h"

class BarnesHutSNEBridge
{

private:
	BarnesHutSNEBridge();
	~BarnesHutSNEBridge();

	static double * loadData(const Eigen::MatrixXd & eigendata);
	static Eigen::MatrixXd saveData(double* data, int n, int d);

public:
	static Eigen::MatrixXd runBarnesHutSNE(const Eigen::MatrixXd & eigendata, const Opts & opts);
	
};

#endif