/* RBGL2.h -- R interface to Boost Graph Library header */
/* assumes -IboostIncl for RBGL/src */

#include <boost/config.hpp>

#include <iostream>
#include <stdio.h>
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

using namespace std;
using namespace boost;

/* Graph_di is directed with integer weights
/* Graph_ui is undirected with integer weights
/* Graph_dd is directed with double weights
/* Graph_ud is undirected with double weights */

typedef adjacency_list<vecS, vecS, directedS,
                       property<vertex_color_t, default_color_type>,
                       property<edge_weight_t, int> > Graph_di;
                                      
typedef adjacency_list<vecS, vecS, undirectedS,
                       property<vertex_color_t, default_color_type>,
                       property<edge_weight_t, int> > Graph_ui;
typedef adjacency_list<vecS, vecS, directedS,
                              property<vertex_color_t, default_color_type>,
                              property<edge_weight_t, double> > Graph_dd;
typedef adjacency_list<vecS, vecS, undirectedS,
                       property<vertex_color_t, default_color_type>,
                       property<edge_weight_t, double> > Graph_ud;

template <class GTYPE>
inline GTYPE
BGL_createGraphIntWeighted(SEXP num_verts_in,
                           SEXP num_edges_in,
                           SEXP R_edges_in,
                           SEXP R_weights_in)
{
    typedef std::pair<int, int> E;
    if (!isInteger(R_weights_in)) error("R_weights_in should be integer");
    if (!isInteger(R_edges_in)) error("R_edges_in should be integer");
    int num_nodes = asInteger(num_verts_in);
    int NE = asInteger(num_edges_in);
    std::vector<E> edge_vector;
    edge_vector.reserve(NE);
    std::vector<int> weights(INTEGER(R_weights_in),
                             INTEGER(R_weights_in)+NE);
    for (int i = 0, cur = 0; i < NE ; i++, cur += 2) {
        edge_vector.push_back(std::make_pair((int) INTEGER(R_edges_in)[cur], 
                                             (int) INTEGER(R_edges_in)[cur+1]));
    }

    return GTYPE(edge_vector.begin(), edge_vector.end(),
                 weights.begin(), num_nodes);
}

template <class GTYPE>
inline GTYPE
BGL_createGraphDoubleWeighted(SEXP num_verts_in,
                              SEXP num_edges_in,
                              SEXP R_edges_in,
                              SEXP R_weights_in)
{
    typedef std::pair<int, int> E;
    if (!isNumeric(R_weights_in)) error("R_weights_in should be double");
    if (!isInteger(R_edges_in)) error("R_edges_in should be integer");
    int num_nodes = asInteger(num_verts_in);
    int NE = asInteger(num_edges_in);
    std::vector<E> edge_vector;
    std::vector<double> weights;
    edge_vector.reserve(NE);
    weights.reserve(NE);
    if (isReal(R_weights_in)) {
        weights.insert(weights.end(),
                       REAL(R_weights_in), REAL(R_weights_in)+NE);
    } else {
        weights.insert(weights.end(),
                       INTEGER(R_weights_in), INTEGER(R_weights_in)+NE);
    }
    for (int i = 0, cur = 0; i < NE ; i++, cur += 2) {
        edge_vector.push_back(std::make_pair((int) INTEGER(R_edges_in)[cur], 
                                             (int) INTEGER(R_edges_in)[cur+1]));
    }

    return GTYPE(edge_vector.begin(), edge_vector.end(),
                 weights.begin(), num_nodes);
}

struct unity_generator {
    typedef double result_type;
    double operator()() { return 1.0; }
};

template <class GTYPE>
inline GTYPE
BGL_createGraphUnweighted(SEXP num_verts_in,
                          SEXP num_edges_in,
                          SEXP R_edges_in)
{
    typedef std::pair<int, int> E;
    if (!isInteger(R_edges_in)) error("R_edges_in should be integer");
    int num_nodes = asInteger(num_verts_in);
    int NE = asInteger(num_edges_in);
    std::vector<E> edge_vector;
    edge_vector.reserve(NE);
    for (int i = 0, cur = 0; i < NE ; i++, cur += 2) {
        edge_vector.push_back(std::make_pair((int) INTEGER(R_edges_in)[cur], 
                                             (int) INTEGER(R_edges_in)[cur+1]));
    }

    unity_generator gen;

    return GTYPE(edge_vector.begin(), edge_vector.end(),
                 make_generator_iterator(gen), num_nodes);
}


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

