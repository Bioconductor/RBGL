/* RBGL2.h -- R interface to Boost Graph Library header */
/* assumes -IboostIncl for RBGL/src */

#ifndef RBGL_RBGL_H
#define RBGL_RBGL_H

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
#include <boost/pending/integer_range.hpp>
#include <boost/pending/indirect_cmp.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/graph/johnson_all_pairs_shortest.hpp>
#include <boost/graph/edge_connectivity.hpp>

/* added vc 7 aug 04 */
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/transitive_closure.hpp>
#include <boost/property_map.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/pending/ct_if.hpp>
#include <boost/type_traits/same_traits.hpp>
/* END added vc 7 aug 04 */


extern "C" {
#include <Rinternals.h>
}

template <class DirectedS = boost::directedS, typename WeightT = double>
class R_adjacency_list
            : public boost::adjacency_list<boost::vecS, boost::vecS, DirectedS,
            boost::property<boost::vertex_color_t, boost::default_color_type>,
            boost::property<boost::edge_weight_t, WeightT> >
{
    typedef boost::adjacency_list<boost::vecS, boost::vecS, DirectedS,
    boost::property<boost::vertex_color_t, boost::default_color_type>,
    boost::property<boost::edge_weight_t, WeightT> > Base;
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
        int NE = asInteger(num_edges_in);
        int* edges_in = INTEGER(R_edges_in);
        for (int i = 0; i < NE ; i++, edges_in += 2) {
            boost::add_edge(*edges_in, *(edges_in+1), 1, *this);
        }
    }
};

/* Graph_di is directed with integer weights
   Graph_ui is undirected with integer weights
   Graph_dd is directed with double weights
   Graph_ud is undirected with double weights */

typedef R_adjacency_list<boost::directedS, int> Graph_di;
typedef R_adjacency_list<boost::undirectedS, int> Graph_ui;
typedef R_adjacency_list<boost::directedS, double> Graph_dd;
typedef R_adjacency_list<boost::undirectedS, double> Graph_ud;

#endif // RBGL_RBGL_H

