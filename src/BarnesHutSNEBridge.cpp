#include "BarnesHutSNEBridge.h"

#include "Util.h"
#include <iostream>

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
Eigen::MatrixXd BarnesHutSNEBridge::runBarnesHutSNE(const Eigen::MatrixXd & eigendata, const Opts & opts)
{

    TSNE* tsne = new TSNE();

    int N = eigendata.rows();
    int D = eigendata.cols();
    
    // estimate perplexity if necessary
    unsigned perplexity = opts.tsnePerplexity();
    if (perplexity == 0)
    {
        perplexity = Util::estimateTsnePerplexity(eigendata);
        std::cout << "Estimated perplexity: " << perplexity << std::endl;
    }

    // read the parameters and the dataset
	double * data = BarnesHutSNEBridge::loadData(eigendata);
	
	// fire up the SNE implementation
	double * Y = (double*) malloc(N * opts.tsneDim() * sizeof(double));
    if (Y == NULL)
	{ 
		printf("Memory allocation failed!\n"); 
		exit(1); 
	}
    srand(time(NULL));
	tsne->run(data, N, D, Y, opts.tsneDim(), perplexity, opts.tsneTheta());
	
	// save the results
	Eigen::MatrixXd result = BarnesHutSNEBridge::saveData(Y, N, opts.tsneDim());
    
    // clean up
	free(data); 
	free(Y);
    
    delete tsne;

    return result;

}
