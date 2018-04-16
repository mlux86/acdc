#pragma once

#include <Eigen/Dense>

class BarnesHutSNEAdapter
{

private:
	BarnesHutSNEAdapter();
	~BarnesHutSNEAdapter();

	// Estimates t-Sne perplexity by a rule of thumb
	static unsigned estimateTsnePerplexity(const Eigen::MatrixXd & mat);

	// Loads t-Sne matrix from an Eigen object
	static double * loadData(const Eigen::MatrixXd & eigendata);

	// Copies t-Sne matrix to an Eigen object
	static Eigen::MatrixXd saveData(double* data, int n, int d);

public:
	// Run BH-Sne on the eigendata object
	// Copied and modified from the original BH-SNE implementation
	static Eigen::MatrixXd runBarnesHutSNE(const Eigen::MatrixXd & eigendata, unsigned seed);
	
};
