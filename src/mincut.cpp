#include "RBGL.hpp"
#include "mincut.hpp"
#include <boost/graph/edmonds_karp_max_flow.hpp>
#include <boost/graph/push_relabel_max_flow.hpp>
#include <boost/graph/boykov_kolmogorov_max_flow.hpp> 

typedef enum { E_MF_Push_Relabel, E_MF_Edmonds_Karp, E_MF_Kolmogorov } E_MF_METHOD;

static SEXP BGL_max_flow_internal(SEXP num_verts_in, SEXP num_edges_in,
                                  SEXP R_edges_in, SEXP R_capacity_in,
                                  SEXP src, SEXP sink, E_MF_METHOD method )
{
    using namespace boost;

    typedef adjacency_list_traits<vecS, vecS, directedS> Tr;
    typedef Tr::edge_descriptor Tr_edge_desc;

    typedef adjacency_list<vecS, vecS, directedS, no_property,
    		property<edge_capacity_t, double,
    		property<edge_residual_capacity_t, double,
    		property<edge_reverse_t, Tr_edge_desc> > > >
		FlowGraph;

    typedef graph_traits<FlowGraph>::edge_iterator   edge_iterator;
    typedef graph_traits<FlowGraph>::vertex_descriptor vertex_descriptor;

    FlowGraph flow_g;

    property_map < FlowGraph, edge_capacity_t >::type
    cap = get(edge_capacity, flow_g);
    property_map < FlowGraph, edge_residual_capacity_t >::type
    res_cap = get(edge_residual_capacity, flow_g);
    property_map < FlowGraph, edge_reverse_t >::type
    rev_edge = get(edge_reverse, flow_g);

    Tr::edge_descriptor e1, e2;
    bool in1, in2;

    if (!isInteger(R_edges_in)) error("R_edges_in should be integer");

    int NV = INTEGER(num_verts_in)[0];
    int NE = asInteger(num_edges_in);
    int* edges_in = INTEGER(R_edges_in);
    int*    capacity_i = (isReal(R_capacity_in)) ? 0 : INTEGER(R_capacity_in);
    double* capacity_d = (isReal(R_capacity_in)) ? REAL(R_capacity_in) : 0;

    for (int i = 0; i < NE ; i++, edges_in += 2)
    {
        tie(e1, in1) = boost::add_edge(*edges_in, *(edges_in+1), flow_g);
        tie(e2, in2) = boost::add_edge(*(edges_in+1), *edges_in, flow_g);
        if ( !in1 || !in2 )
            error("unable to add edge: (%d, %d)", *edges_in, *(edges_in+1));

        // fill in capacity_map
        cap[e1] = capacity_i ? (*capacity_i++) : (*capacity_d++);
        cap[e2] = 0;

        // fill in reverse_edge_map
        rev_edge[e1] = e2;
        rev_edge[e2] = e1;
    }

    double maxflow = 0;

    int vsrc = INTEGER(src)[0];
    int vsink = INTEGER(sink)[0];

    if ( 0 <= vsrc && vsrc < NV && 0 <= vsink && vsink < NV )
    {
    	vertex_descriptor s = vertex(vsrc, flow_g);
    	vertex_descriptor t = vertex(vsink, flow_g);

    	if ( method == E_MF_Push_Relabel ) 
	     maxflow = push_relabel_max_flow(flow_g, s, t);
	else if ( method == E_MF_Edmonds_Karp )
             maxflow = edmonds_karp_max_flow(flow_g, s, t);

/*
needs different parameters as of 14 june 2012
	else if ( method == E_MF_Kolmogorov )
	{
//	     error("kolmogorov_max_flow from BGL doesn't work");
	     maxflow = boykov_kolmogorov_max_flow(flow_g, s, t);
	}
*/

	else
	     error("unknown method for max_flow");
    }

    SEXP ansList, conn, eList, fList;
    PROTECT(ansList = allocVector(VECSXP,3));
    PROTECT(conn = NEW_NUMERIC(1));
    PROTECT(eList = allocMatrix(INTSXP, 2, asInteger(num_edges_in)));
    PROTECT(fList = allocMatrix(REALSXP, 1, asInteger(num_edges_in)));

    REAL(conn)[0] = maxflow;

    edge_iterator ei, e_end;
    int i = 0, j = 0;
    for (tie(ei, e_end) = edges(flow_g); ei != e_end; ++ei)
        if (cap[*ei] > 0) {
            INTEGER(eList)[i++] = source(*ei, flow_g);
            INTEGER(eList)[i++] = target(*ei, flow_g);
            REAL(fList)[j++] = cap[*ei] - res_cap[*ei];
        }

    SET_VECTOR_ELT(ansList,0,conn);
    SET_VECTOR_ELT(ansList,1,eList);
    SET_VECTOR_ELT(ansList,2,fList);
    UNPROTECT(4);

    return(ansList);
}

extern "C"
{
    SEXP BGL_min_cut_U (SEXP num_verts_in, SEXP num_edges_in,
                        SEXP R_edges_in,   SEXP R_weights_in )
    {
        using namespace boost;

        Graph_ud g(num_verts_in, num_edges_in, R_edges_in, R_weights_in);

        typedef graph_traits < Graph_ud >::vertex_descriptor Vertex;
        typedef graph_traits < Graph_ud >::edges_size_type dst;
        std::vector<Vertex> s_set, vs_set;
        dst cut_capacity = min_cut(g, std::back_inserter(s_set), std::back_inserter(vs_set) );

        SEXP ansList, conn;
        SEXP sList, vsList;  // for subsets of nodes in mincut: S, V - S

        PROTECT(ansList = allocVector(VECSXP,3));
        PROTECT(conn = NEW_NUMERIC(1));
        PROTECT(sList = allocVector(INTSXP, s_set.size()));
        PROTECT(vsList = allocVector(INTSXP,vs_set.size()));

        REAL(conn)[0] = (double)cut_capacity;

        std::vector<Vertex>::iterator vi;
        int i = 0;
        for ( i = 0, vi = s_set.begin(); vi != s_set.end(); i++, ++vi)
            INTEGER(sList)[i] = *vi;

        for ( i = 0, vi = vs_set.begin(); vi != vs_set.end(); i++, ++vi)
            INTEGER(vsList)[i] = *vi;

        SET_VECTOR_ELT(ansList,0,conn);
        SET_VECTOR_ELT(ansList,1,sList);
        SET_VECTOR_ELT(ansList,2,vsList);
        UNPROTECT(4);

        return(ansList);
    }


    SEXP BGL_edmonds_karp_max_flow(SEXP num_verts_in, SEXP num_edges_in,
                                   SEXP R_edges_in, SEXP R_capacity_in, 
				   SEXP src, SEXP sink )
    {
        SEXP ansList = BGL_max_flow_internal(num_verts_in, num_edges_in,
                                             R_edges_in, R_capacity_in, 
					     src, sink, E_MF_Edmonds_Karp);
        return(ansList);
    }

    SEXP BGL_push_relabel_max_flow(SEXP num_verts_in, SEXP num_edges_in,
                                   SEXP R_edges_in, SEXP R_capacity_in, 
				   SEXP src, SEXP sink )
    {
        SEXP ansList = BGL_max_flow_internal(num_verts_in, num_edges_in,
                                             R_edges_in, R_capacity_in, 
					     src, sink, E_MF_Push_Relabel);
        return(ansList);
    }

    SEXP BGL_kolmogorov_max_flow(SEXP num_verts_in, SEXP num_edges_in,
                                   SEXP R_edges_in, SEXP R_capacity_in, 
				   SEXP src, SEXP sink )
    {
        SEXP ansList = BGL_max_flow_internal(num_verts_in, num_edges_in,
                                             R_edges_in, R_capacity_in, 
					     src, sink, E_MF_Kolmogorov);
        return(ansList);
    }

}

