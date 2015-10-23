#ifndef __TarjansAlgorithm__
#define __TarjansAlgorithm__

#include <stack>
#include <vector>

#include "Clustering.h"

class TarjansAlgorithm
{
private:
	struct Node
	{
		unsigned adjacencyIdx = 0;
		int index = -1;
		int lowLink = -1;
		bool onStack = false;		
	};

	unsigned index = 0;
	std::vector<Node> nodes;
	std::stack<Node*> stack;

	const Eigen::MatrixXd * adjacencies;
	unsigned n;

	std::vector<unsigned> strongConnect(Node & v);

public:
	TarjansAlgorithm();
	~TarjansAlgorithm();
	
	ClusteringResult run(const Eigen::MatrixXd & adjacencies);

};

#endif