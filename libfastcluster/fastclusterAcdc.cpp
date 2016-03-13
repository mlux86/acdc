#define fc_isnan(X) ((X)!=(X))

#include <cstddef>
#include <limits>
#include <algorithm>
#include <new>
#include <exception>
#include <memory>
#include <iostream>

#include "fastcluster.cpp"

#include "fastclusterAcdc.h"

/*
  Convenience class for the output array: automatic counter.
*/
class linkage_output {
private:
  t_float * Z;

public:
  linkage_output(t_float * const Z_)
    : Z(Z_)
  {}

  void append(const t_index node1, const t_index node2, const t_float dist,
              const t_float size) {
    if (node1<node2) {
      *(Z++) = static_cast<t_float>(node1);
      *(Z++) = static_cast<t_float>(node2);
    }
    else {
      *(Z++) = static_cast<t_float>(node2);
      *(Z++) = static_cast<t_float>(node1);
    }
    *(Z++) = dist;
    *(Z++) = size;
  }
};

// The size of a node is either 1 (a single point) or is looked up from
// one of the clusters.
#define size_(r_) ( ((r_<N) ? 1 : Z_(r_-N,3)) )

template <const bool sorted>
static void generate_dendrogram(t_float * const Z, cluster_result & Z2, const t_index N) {
  // The array "nodes" is a union-find data structure for the cluster
  // identities (only needed for unsorted cluster_result input).
  union_find nodes(sorted ? 0 : N);
  if (!sorted) {
    std::stable_sort(Z2[0], Z2[N-1]);
  }

  linkage_output output(Z);
  t_index node1, node2;

  for (node const * NN=Z2[0]; NN!=Z2[N-1]; ++NN) {
    // Get two data points whose clusters are merged in step i.
    if (sorted) {
      node1 = NN->node1;
      node2 = NN->node2;
    }
    else {
      // Find the cluster identifiers for these points.
      node1 = nodes.Find(NN->node1);
      node2 = nodes.Find(NN->node2);
      // Merge the nodes in the union-find data structure by making them
      // children of a new node.
      nodes.Union(node1, node2);
    }
    output.append(node1, node2, NN->dist, size_(node1)+size_(node2));
  }
}

void fc_linkage(double * const D_, double * const Z_, long int N_, unsigned char method) 
{
    
    if (N_ < 1 ) 
    {
        throw std::runtime_error("At least one element is needed for clustering.");
    }

    if (N_ > MAX_INDEX/4 || static_cast<int64_t>(N_-1)>>(T_FLOAT_MANT_DIG-1) > 0) 
    {
        throw std::runtime_error("Data is too big, index overflow.");
    }

    t_index N = static_cast<t_index>(N_);

    cluster_result Z2(N-1);
    auto_array_ptr<t_index> members;

    // For these methods, the distance update formula needs the number of
    // data points in a cluster.
    if (method==METHOD_METR_AVERAGE || method==METHOD_METR_WARD || method==METHOD_METR_CENTROID) 
    {
        members.init(N, 1);
    }

    // Operate on squared distances for these methods.
    if (method==METHOD_METR_WARD || method==METHOD_METR_CENTROID || method==METHOD_METR_MEDIAN) 
    {
        for (t_float * DD = D_; DD!=D_+static_cast<std::ptrdiff_t>(N)*(N-1)/2; ++DD)
        {
            *DD *= *DD;
        }
    }

    switch (method) 
    {
        case METHOD_METR_SINGLE:
        MST_linkage_core(N, D_, Z2);
        break;
        case METHOD_METR_COMPLETE:
        NN_chain_core<METHOD_METR_COMPLETE, t_index>(N, D_, NULL, Z2);
        break;
        case METHOD_METR_AVERAGE:
        NN_chain_core<METHOD_METR_AVERAGE, t_index>(N, D_, members, Z2);
        break;
        case METHOD_METR_WEIGHTED:
        NN_chain_core<METHOD_METR_WEIGHTED, t_index>(N, D_, NULL, Z2);
        break;
        case METHOD_METR_WARD:
        NN_chain_core<METHOD_METR_WARD, t_index>(N, D_, members, Z2);
        break;
        case METHOD_METR_CENTROID:
        generic_linkage<METHOD_METR_CENTROID, t_index>(N, D_, members, Z2);
        break;
        case METHOD_METR_MEDIAN:
        generic_linkage<METHOD_METR_MEDIAN, t_index>(N, D_, NULL, Z2);
        break;
        default:
        throw std::runtime_error(std::string("Invalid method index."));
    }

    if (method==METHOD_METR_WARD || method==METHOD_METR_CENTROID || method==METHOD_METR_MEDIAN) 
    {
        Z2.sqrt();
    }    

    if (method==METHOD_METR_CENTROID || method==METHOD_METR_MEDIAN) 
    {
        generate_dendrogram<true>(Z_, Z2, N);
    }
    else 
    {
        generate_dendrogram<false>(Z_, Z2, N);
    }

}

