/* RBGL2.h -- R interface to Boost Graph Library header */
/* assumes -IboostIncl for RBGL/src */

#include <boost/config.hpp>

#include <iostream>
#include <cstdio>
#include <vector>
#include <iterator>
#include <algorithm>
#include <time.h>

#include <boost/utility.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/pending/integer_range.hpp>
#include <boost/pending/indirect_cmp.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/edmunds_karp_max_flow.hpp>
#include <boost/graph/push_relabel_max_flow.hpp>
#include <boost/generator_iterator.hpp>

extern "C" {
#include <Rinternals.h>
}

template <class DirectedS = boost::directedS,
          typename WeightT = double>
class R_adjacency_list
    : public boost::adjacency_list<boost::vecS, boost::vecS, DirectedS,
                                   boost::property<boost::vertex_color_t,
                                                   boost::default_color_type>,
                                   boost::property<boost::edge_weight_t,
                                                   WeightT> >
{
    typedef boost::adjacency_list<boost::vecS, boost::vecS, DirectedS,
                                  boost::property<boost::vertex_color_t,
                                                  boost::default_color_type>,
                                  boost::property<boost::edge_weight_t,
                                                  WeightT> > Base;
    typedef WeightT R_weight_type;
    BOOST_STATIC_ASSERT(boost::is_arithmetic<R_weight_type>::value);
public:
    typedef typename Base::graph_property_type graph_property_type;
    typedef typename Base::vertices_size_type vertices_size_type;
    typedef typename Base::edges_size_type edges_size_type;

    inline R_adjacency_list()
        : Base() { }
    inline R_adjacency_list(const graph_property_type& p)
        : Base(p) { }
    inline R_adjacency_list(const Base& x)
        : Base(x) { }
    inline R_adjacency_list(vertices_size_type num_vertices)
        : Base(num_vertices) { }
    inline R_adjacency_list(vertices_size_type num_vertices, 
                            const graph_property_type& p)
        : Base(num_vertices, p) { }
#if !defined(BOOST_MSVC) || BOOST_MSVC > 1300
    // Required by Iterator Constructible Graph
    template <class EdgeIterator>
    inline R_adjacency_list(EdgeIterator first, EdgeIterator last,
                            vertices_size_type n,
                            edges_size_type m = 0,
                            const graph_property_type& p = graph_property_type())
        : Base(first, last, n, m, p) { }

    template <class EdgeIterator, class EdgePropertyIterator>
    inline R_adjacency_list(EdgeIterator first, EdgeIterator last,
                            EdgePropertyIterator ep_iter,
                            vertices_size_type n,
                            edges_size_type m = 0,
                            const graph_property_type& p = graph_property_type())
        : Base(first, last, ep_iter, n, m, p) { }
#endif

    inline R_adjacency_list(SEXP num_verts_in,
                            SEXP num_edges_in,
                            SEXP R_edges_in,
                            SEXP R_weights_in)
        : Base(asInteger(num_verts_in))
    {
        if (!isNumeric(R_weights_in)) error("R_weights_in should be Numeric");
        if (!isInteger(R_edges_in)) error("R_edges_in should be integer");
        int num_nodes = asInteger(num_verts_in);
        int NE = asInteger(num_edges_in);
        int* edges_in = INTEGER(R_edges_in);
        if (isReal(R_weights_in)) {
            if (boost::is_integral<R_weight_type>::value)
                error("R_weights_in should be integer");
            else {
                double* weights_in = REAL(R_weights_in);
                for (int i = 0; i < NE ; i++, edges_in += 2, weights_in++) {
                    boost::add_edge(*edges_in, *(edges_in+1),
                                    *weights_in, *this);
                }
            }
        } else {
            int* weights_in = INTEGER(R_weights_in);
            for (int i = 0; i < NE ; i++, edges_in += 2, weights_in++) {
                boost::add_edge(*edges_in, *(edges_in+1), *weights_in, *this);
            }
        }
    }

    inline R_adjacency_list(SEXP num_verts_in,
                            SEXP num_edges_in,
                            SEXP R_edges_in)
        : Base(asInteger(num_verts_in))
    {
        if (!isInteger(R_edges_in)) error("R_edges_in should be integer");
        int num_nodes = asInteger(num_verts_in);
        int NE = asInteger(num_edges_in);
        int* edges_in = INTEGER(R_edges_in);
        for (int i = 0; i < NE ; i++, edges_in += 2) {
            boost::add_edge(*edges_in, *(edges_in+1), 1, *this);
        }
    }
};

/* Graph_di is directed with integer weights
/* Graph_ui is undirected with integer weights
/* Graph_dd is directed with double weights
/* Graph_ud is undirected with double weights */

typedef R_adjacency_list<boost::directedS, int> Graph_di;
typedef R_adjacency_list<boost::undirectedS, int> Graph_ui;
typedef R_adjacency_list<boost::directedS, double> Graph_dd;
typedef R_adjacency_list<boost::undirectedS, double> Graph_ud;


namespace boost
{
  template < typename Graph >
    std::pair < typename graph_traits < Graph >::vertex_descriptor,
    typename graph_traits < Graph >::degree_size_type >
    min_degree_vertex(Graph & g)
  {
    typename graph_traits < Graph >::vertex_descriptor p;
    typedef typename graph_traits < Graph >::degree_size_type size_type;
    size_type delta = std::numeric_limits < size_type >::max();
    typename graph_traits < Graph >::vertex_iterator i, iend;
    for (tie(i, iend) = vertices(g); i != iend; ++i)
      if (degree(*i, g) < delta)
      {
        delta = degree(*i, g);
        p = *i;
      }
    return std::make_pair(p, delta);
  }

  template < typename Graph, typename OutputIterator >
    void neighbors(const Graph & g,
                   typename graph_traits < Graph >::vertex_descriptor u,
                   OutputIterator result)
  {
    typename graph_traits < Graph >::adjacency_iterator ai, aend;
    for (tie(ai, aend) = adjacent_vertices(u, g); ai != aend; ++ai)
      *result++ = *ai;
  }
  template < typename Graph, typename VertexIterator,
    typename OutputIterator > void neighbors(const Graph & g,
                                             VertexIterator first,
                                             VertexIterator last,
                                             OutputIterator result)
  {
    for (; first != last; ++first)
      neighbors(g, *first, result);
  }

  template < typename VertexListGraph, typename OutputIterator >
  typename graph_traits < VertexListGraph >::degree_size_type
  edge_connectivity(VertexListGraph & g, OutputIterator disconnecting_set)
  {
    typedef typename graph_traits <
      VertexListGraph >::vertex_descriptor vertex_descriptor;
    typedef typename graph_traits <
      VertexListGraph >::degree_size_type degree_size_type;
    typedef color_traits < default_color_type > Color;
    typedef typename adjacency_list_traits < vecS, vecS,
      directedS >::edge_descriptor edge_descriptor;
    typedef adjacency_list < vecS, vecS, directedS, no_property,
      property < edge_capacity_t, degree_size_type,
      property < edge_residual_capacity_t, degree_size_type,
      property < edge_reverse_t, edge_descriptor > > > > FlowGraph;

    vertex_descriptor u, v, p, k;
    edge_descriptor e1, e2;
    bool inserted;
    typename graph_traits < VertexListGraph >::vertex_iterator vi, vi_end;
    degree_size_type delta, alpha_star, alpha_S_k;
    std::set < vertex_descriptor > S, neighbor_S;
    std::vector < vertex_descriptor > S_star, nonneighbor_S;
    std::vector < default_color_type > color(num_vertices(g));
    std::vector < edge_descriptor > pred(num_vertices(g));

    FlowGraph flow_g(num_vertices(g));
    typename property_map < FlowGraph, edge_capacity_t >::type
      cap = get(edge_capacity, flow_g);
    typename property_map < FlowGraph, edge_residual_capacity_t >::type
      res_cap = get(edge_residual_capacity, flow_g);
    typename property_map < FlowGraph, edge_reverse_t >::type
      rev_edge = get(edge_reverse, flow_g);

    typename graph_traits < VertexListGraph >::edge_iterator ei, ei_end;
    for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
      u = source(*ei, g), v = target(*ei, g);
      tie(e1, inserted) = add_edge(u, v, flow_g);
      cap[e1] = 1;
      tie(e2, inserted) = add_edge(v, u, flow_g);
      cap[e2] = 1;
      rev_edge[e1] = e2;
      rev_edge[e2] = e1;
    }

    tie(p, delta) = min_degree_vertex(g);
    S_star.push_back(p);
    alpha_star = delta;
    S.insert(p);
    neighbor_S.insert(p);
    neighbors(g, S.begin(), S.end(),
              std::inserter(neighbor_S, neighbor_S.begin()));
    std::set_difference(vertices(g).first, vertices(g).second,
                        neighbor_S.begin(), neighbor_S.end(),
                        std::back_inserter(nonneighbor_S));

    while (!nonneighbor_S.empty()) {
      k = nonneighbor_S.front();
      alpha_S_k = edmunds_karp_max_flow
        (flow_g, p, k, cap, res_cap, rev_edge, &color[0], &pred[0]);
      if (alpha_S_k < alpha_star) {
        alpha_star = alpha_S_k;
        S_star.clear();
        for (tie(vi, vi_end) = vertices(flow_g); vi != vi_end; ++vi)
          if (color[*vi] != Color::white())
            S_star.push_back(*vi);
      }
      S.insert(k);
      neighbor_S.insert(k);
      neighbors(g, k, std::inserter(neighbor_S, neighbor_S.begin()));
      nonneighbor_S.clear();
      std::set_difference(vertices(g).first, vertices(g).second,
                          neighbor_S.begin(), neighbor_S.end(),
                          std::back_inserter(nonneighbor_S));
    }

    std::vector < bool > in_S_star(num_vertices(g), false);
    typename std::vector < vertex_descriptor >::iterator si;
    for (si = S_star.begin(); si != S_star.end(); ++si)
      in_S_star[*si] = true;
    degree_size_type c = 0;
    for (si = S_star.begin(); si != S_star.end(); ++si) {
      typename graph_traits < VertexListGraph >::out_edge_iterator ei, ei_end;
      for (tie(ei, ei_end) = out_edges(*si, g); ei != ei_end; ++ei)
        if (!in_S_star[target(*ei, g)]) {
          *disconnecting_set++ = *ei;
          ++c;
        }
    }

    return c;
  }

}

