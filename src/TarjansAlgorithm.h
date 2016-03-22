#pragma once

#include <stack>
#include <vector>
#include <Eigen/Dense>

// forward declaration
struct ClusteringResult; 

// Finding strongly connected components of a graph
// Algorithm resembles Wikipedia Pseudo code:
// https://en.wikipedia.org/wiki/Tarjan%27s_strongly_connected_components_algorithm
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
