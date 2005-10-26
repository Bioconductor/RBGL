#ifndef RBGL_MINCUT_H
#define RBGL_MINCUT_H

#include <boost/graph/edge_connectivity.hpp>

namespace boost 
{

template < typename VertexListGraph, typename OutputIterator >
typename graph_traits<VertexListGraph>::edges_size_type
min_cut(VertexListGraph & g, OutputIterator s_set, OutputIterator vs_set)
{
    typedef graph_traits<VertexListGraph> Traits;
    typedef typename Traits::vertex_descriptor vertex_descriptor;
    typedef typename Traits::vertex_iterator   vertex_iterator;
    typedef typename Traits::degree_size_type  degree_size_type;
    typedef typename Traits::edge_iterator     edge_iterator;
    typedef typename Traits::edges_size_type   edges_size_type;

    typedef color_traits<default_color_type> Color;

    typedef adjacency_list_traits<vecS, vecS, directedS> Tr;
    typedef typename Tr::edge_descriptor Tr_edge_desc;
    typedef adjacency_list<vecS, vecS, directedS, no_property,
	          property<edge_capacity_t, degree_size_type,
	          property<edge_residual_capacity_t, degree_size_type,
	          property<edge_reverse_t, Tr_edge_desc> > > >
	    FlowGraph;
    typedef typename graph_traits<FlowGraph>::edge_descriptor edge_descriptor;

    vertex_descriptor u, v, p, k;
    vertex_iterator vi, vi_end;
    edge_descriptor e1, e2;
    edge_iterator ei, ei_end;
    edges_size_type delta, alpha_star, alpha_S_k;
    bool inserted;

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

    // This loop is to convert undirected graph to directed graph, 
    // each edge in undirected graph has unit capacity 
    // each arc in directed graph also has unit capacity 
    for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei) {
        u = source(*ei, g), v = target(*ei, g);
        tie(e1, inserted) = add_edge(u, v, flow_g);
        cap[e1] = 1;
        tie(e2, inserted) = add_edge(v, u, flow_g);
        cap[e2] = 1;
        rev_edge[e1] = e2;
        rev_edge[e2] = e1;
    }

    tie(p, delta) = detail::min_degree_vertex(g);
    S_star.push_back(p);
    alpha_star = delta;
    S.insert(p);
    neighbor_S.insert(p);
    detail::neighbors(g, S.begin(), S.end(),
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
	detail::neighbors(g, k, std::inserter(neighbor_S, neighbor_S.begin()));
        nonneighbor_S.clear();
        std::set_difference(vertices(g).first, vertices(g).second,
                            neighbor_S.begin(), neighbor_S.end(),
                            std::back_inserter(nonneighbor_S));
    }

    /* To extract mincut sets: S and (V - S) */
    std::vector < bool > in_S_star(num_vertices(g), false);
    typename std::vector < vertex_descriptor >::iterator si;
    for (si = S_star.begin(); si != S_star.end(); ++si)
          in_S_star[*si] = true;
    
    for (tie(vi, vi_end) = vertices(flow_g); vi != vi_end; ++vi)
	  if ( in_S_star[*vi] ) *s_set++ = *vi; 
          else *vs_set++ = *vi; 

    return alpha_star;
}

}
#endif  // RBGL_MINCUT_H

