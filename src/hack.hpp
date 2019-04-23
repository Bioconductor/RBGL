#ifndef __HACK__
#define __HACK__

#include <boost/graph/graph_utility.hpp>
#include <boost/graph/betweenness_centrality.hpp>

using namespace std;
using namespace boost;

template <class _ForwardIterator, class _Compare>
_ForwardIterator
max_kludge_element(_ForwardIterator __first, _ForwardIterator __last, _Compare __comp)
{
    if (__first != __last)
    {
        _ForwardIterator __i = __first;
        while (++__i != __last)
            if (__comp(*__first, *__i))
                __first = __i;
    }
    return __first;
}


template <class _ForwardIterator>
_ForwardIterator
max_kludge_element(_ForwardIterator __first, _ForwardIterator __last)
{
    return _VSTD::max_element(__first, __last,
              __less<typename iterator_traits<_ForwardIterator>::value_type>());
}

/** Threshold termination function for the betweenness centrality
 * clustering algorithm.
 */
template<typename T>
struct bc_clustering_threshold
{
  typedef T centrality_type;

  /// Terminate clustering when maximum absolute edge centrality is
  /// below the given threshold.
  explicit bc_clustering_threshold(T threshold) 
    : threshold(threshold), dividend(1.0) {}
  
  /**
   * Terminate clustering when the maximum edge centrality is below
   * the given threshold.
   *
   * @param threshold the threshold value
   *
   * @param g the graph on which the threshold will be calculated
   *
   * @param normalize when true, the threshold is compared against the
   * normalized edge centrality based on the input graph; otherwise,
   * the threshold is compared against the absolute edge centrality.
   */
  template<typename Graph>
  bc_clustering_threshold(T threshold, const Graph& g, bool normalize = true)
    : threshold(threshold), dividend(1.0)
  {
    if (normalize) {
      typename graph_traits<Graph>::vertices_size_type n = num_vertices(g);
      dividend = T((n - 1) * (n - 2)) / T(2);
    }
  }

  /** Returns true when the given maximum edge centrality (potentially
   * normalized) falls below the threshold.
   */
  template<typename Graph, typename Edge>
  bool operator()(T max_centrality, Edge, const Graph&)
  {
    return (max_centrality / dividend) < threshold;
  }

 protected:
  T threshold;
  T dividend;
};

/** Graph clustering based on edge betweenness centrality.
 * 
 * This algorithm implements graph clustering based on edge
 * betweenness centrality. It is an iterative algorithm, where in each
 * step it compute the edge betweenness centrality (via @ref
 * brandes_betweenness_centrality) and removes the edge with the
 * maximum betweenness centrality. The @p done function object
 * determines when the algorithm terminates (the edge found when the
 * algorithm terminates will not be removed).
 *
 * @param g The graph on which clustering will be performed. The type
 * of this parameter (@c MutableGraph) must be a model of the
 * VertexListGraph, IncidenceGraph, EdgeListGraph, and Mutable Graph
 * concepts.
 *
 * @param done The function object that indicates termination of the
 * algorithm. It must be a ternary function object thats accepts the
 * maximum centrality, the descriptor of the edge that will be
 * removed, and the graph @p g.
 *
 * @param edge_centrality (UTIL/OUT) The property map that will store
 * the betweenness centrality for each edge. When the algorithm
 * terminates, it will contain the edge centralities for the
 * graph. The type of this property map must model the
 * ReadWritePropertyMap concept. Defaults to an @c
 * iterator_property_map whose value type is 
 * @c Done::centrality_type and using @c get(edge_index, g) for the 
 * index map.
 *
 * @param vertex_index (IN) The property map that maps vertices to
 * indices in the range @c [0, num_vertices(g)). This type of this
 * property map must model the ReadablePropertyMap concept and its
 * value type must be an integral type. Defaults to 
 * @c get(vertex_index, g).
 */
template<typename MutableGraph, typename Done, typename EdgeCentralityMap,
         typename VertexIndexMap>
void 
betweenness_centrality_clustering(MutableGraph& g, Done done,
                                  EdgeCentralityMap edge_centrality,
                                  VertexIndexMap vertex_index)
{
  typedef typename property_traits<EdgeCentralityMap>::value_type
    centrality_type;
  typedef typename graph_traits<MutableGraph>::edge_iterator edge_iterator;
  typedef typename graph_traits<MutableGraph>::edge_descriptor edge_descriptor;

  if (has_no_edges(g)) return;

  // Function object that compares the centrality of edges
  indirect_cmp<EdgeCentralityMap, std::less<centrality_type> > 
    cmp(edge_centrality);

  bool is_done;
  do {
    brandes_betweenness_centrality(g, 
                                   edge_centrality_map(edge_centrality)
                                   .vertex_index_map(vertex_index));
    std::pair<edge_iterator, edge_iterator> edges_iters = edges(g);
    edge_descriptor e = *max_kludge_element(edges_iters.first, edges_iters.second, cmp);
    is_done = done(get(edge_centrality, e), e, g);
    if (!is_done) remove_edge(e, g);
  } while (!is_done && !has_no_edges(g));
}

/**
 * \overload
 */ 
template<typename MutableGraph, typename Done, typename EdgeCentralityMap>
void 
betweenness_centrality_clustering(MutableGraph& g, Done done,
                                  EdgeCentralityMap edge_centrality)
{
  betweenness_centrality_clustering(g, done, edge_centrality,
                                    get(vertex_index, g));
}

/**
 * \overload
 */ 
template<typename MutableGraph, typename Done>
void
betweenness_centrality_clustering(MutableGraph& g, Done done)
{
  typedef typename Done::centrality_type centrality_type;
  std::vector<centrality_type> edge_centrality(num_edges(g));
  betweenness_centrality_clustering(g, done, 
    make_iterator_property_map(edge_centrality.begin(), get(edge_index, g)),
    get(vertex_index, g));
}

#endif 
