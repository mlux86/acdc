#include "TarjansAlgorithm.h"

TarjansAlgorithm::TarjansAlgorithm()
{
}

TarjansAlgorithm::~TarjansAlgorithm()
{
}

std::vector<unsigned> TarjansAlgorithm::strongConnect(Node & v)
{
	v.index = index;
	v.lowLink = index;
	index++;
	stack.push(&v);
	v.onStack = true;

	// consider successors of v
	for (unsigned e = 0; e < n; e++)
	{
		if ((*adjacencies)(v.adjacencyIdx, e) > 0) // edge
		{
			Node & w = nodes[e];
			if (w.index == -1)
			{
				strongConnect(w);
				v.lowLink = std::min(v.lowLink, w.lowLink);
			} else if (w.onStack)
			{
				v.lowLink = std::min(v.lowLink, w.index);
			}
		}
	}

	// If v is a root node, pop the stack and generate a strongly connected component
	std::vector<unsigned> componentIndexes;
	if (v.lowLink == v.index)
	{
		Node * w;
		do
		{
			w = stack.top();
			stack.pop();
			w->onStack = false;
			componentIndexes.push_back(w->adjacencyIdx);
		} while (w != &v);

	}

	return componentIndexes;
}

ClusteringResult TarjansAlgorithm::run(const Eigen::MatrixXd & adjacencies_)
{
	adjacencies = &adjacencies_;

	n = adjacencies->rows();
	if (n != adjacencies->cols())
	{
		throw std::runtime_error("Adjacency matrix must be quadratic!");
	}

	index = 0;

	nodes.clear();
	for (unsigned i = 0; i < n; i++)
	{
		nodes.emplace_back();
		nodes.back().adjacencyIdx = i;
	}

	stack = std::stack<Node*>();

	ClusteringResult res;
	res.labels.reserve(n);

	for (auto & node : nodes)
	{
		if (node.index == -1)
		{
			std::vector<unsigned> componentIndexes = strongConnect(node);
			if (componentIndexes.size() > 0)
			{
				res.numClusters++;				
				for (unsigned idx : componentIndexes)
				{
					res.labels[idx] = res.numClusters;
				}
			}
		}
	}

	return res; 	
}
