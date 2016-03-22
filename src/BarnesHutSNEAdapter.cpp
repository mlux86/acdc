#include "Logger.h"
#include "BarnesHutSNEAdapter.h"
#include "MLUtil.h"
#include "Opts.h"
#include "tsne.h"

BarnesHutSNEAdapter::BarnesHutSNEAdapter()
{
}

BarnesHutSNEAdapter::~BarnesHutSNEAdapter()
{
}

unsigned BarnesHutSNEAdapter::estimateTsnePerplexity(const Eigen::MatrixXd & mat)
{
    return (unsigned) (pow(log(mat.rows()), 2)); 
}

double * BarnesHutSNEAdapter::loadData(const Eigen::MatrixXd & eigendata)
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

Eigen::MatrixXd BarnesHutSNEAdapter::saveData(double* data, int N, int D)
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

Eigen::MatrixXd BarnesHutSNEAdapter::runBarnesHutSNE(const Eigen::MatrixXd & eigendata)
{
    Eigen::MatrixXd tsneData;
    if (eigendata.cols() > Opts::tsnePcaDim())
    {
        DLOG << "Initially reducing dimension to " << Opts::tsnePcaDim() << " using PCA...\n";
        tsneData = MLUtil::pca(eigendata, Opts::tsnePcaDim());
    } else
    {
        tsneData = eigendata;
    }

    TSNE* tsne = new TSNE();

    int N = tsneData.rows();
    int D = tsneData.cols();
    
    // estimate perplexity if necessary
    unsigned perplexity = Opts::tsnePerplexity();
    if (perplexity == 0)
    {
        perplexity = BarnesHutSNEAdapter::estimateTsnePerplexity(tsneData);
    }

    // read the parameters and the dataset
	double * data = BarnesHutSNEAdapter::loadData(tsneData);
	
	// fire up the SNE implementation
	double * Y = (double*) malloc(N * Opts::tsneDim() * sizeof(double));
    if (Y == NULL)
	{ 
		throw std::runtime_error("Memory allocation failed.");
	}
    srand(time(NULL));
    
    DLOG << "n=" << N << "   " 
            << "dim=" << D << "   " 
            << "targetDim=" << Opts::tsneDim() << "   "
            << "perplexity=" << perplexity << "   "
            << "theta=" << Opts::tsneTheta() << "\n";
	tsne->run(data, N, D, Y, Opts::tsneDim(), perplexity, Opts::tsneTheta());
	
	// save the results
	Eigen::MatrixXd result = BarnesHutSNEAdapter::saveData(Y, N, Opts::tsneDim());
    
    // clean up
	free(data); 
	free(Y);
    
    delete tsne;

    return result;

}
