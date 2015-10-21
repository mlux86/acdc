#include "BarnesHutSNEBridge.h"

BarnesHutSNEBridge::BarnesHutSNEBridge()
{
}

BarnesHutSNEBridge::~BarnesHutSNEBridge()
{
}

double * BarnesHutSNEBridge::loadData(const Eigen::MatrixXd & eigendata)
{
	int N = eigendata.rows();
	int D = eigendata.cols();
	double * data = (double *) malloc(D * N * sizeof(double));

	int nD = 0;
    for(int n = 0; n < N; n++) 
    {
        for(int d = 0; d < D; d++) 
        {
            data[nD + d] = eigendata(n, d);
        }
        nD += D;
    }

    return data;
}

Eigen::MatrixXd BarnesHutSNEBridge::saveData(double* data, int N, int D)
{
	Eigen::MatrixXd eigendata(N, D);

	int nD = 0;
    for(int n = 0; n < N; n++) 
    {
        for(int d = 0; d < D; d++) 
        {
            eigendata(n, d) = data[nD + d];
        }
        nD += D;
    }

	return eigendata;	
}

// copied and modified from the original BH-SNE implementation
Eigen::MatrixXd BarnesHutSNEBridge::runBarnesHutSNE(const Eigen::MatrixXd & eigendata, const unsigned targetDim, const double theta, const double perplexity)
{

    TSNE* tsne = new TSNE();

    int N = eigendata.rows();
    int D = eigendata.cols();
    
    // Read the parameters and the dataset
	double * data = BarnesHutSNEBridge::loadData(eigendata);
	
	// Fire up the SNE implementation
	double * Y = (double*) malloc(N * targetDim * sizeof(double));
    if (Y == NULL)
	{ 
		printf("Memory allocation failed!\n"); 
		exit(1); 
	}
    srand(time(NULL));
	tsne->run(data, N, D, Y, targetDim, perplexity, theta);
	
	// Save the results
	Eigen::MatrixXd result = BarnesHutSNEBridge::saveData(Y, N, targetDim);
    
    // Clean up
	free(data); 
	free(Y);
    
    delete tsne;

    return result;

}
