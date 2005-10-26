#include "RBGL.hpp"
#include <boost/graph/circle_layout.hpp>
#include <boost/graph/kamada_kawai_spring_layout.hpp>

using namespace boost;

enum vertex_position_t { vertex_position };
namespace boost { BOOST_INSTALL_PROPERTY(vertex, position); }

struct point { double x, y; };

typedef enum { E_LAYOUT_CIRCLE, E_LAYOUT_KKSL } E_LAYOUT_METHOD;

extern "C"
{
    typedef adjacency_list<vecS, vecS, undirectedS,
    // vertex properties
    property<vertex_index_t, int,
    property<vertex_position_t, point> >,
    // edge properties
    property<edge_weight_t, double> >
    IndexGraph;

    SEXP BGL_layout_internal
    (E_LAYOUT_METHOD method, SEXP num_verts_in, SEXP num_edges_in, SEXP R_edges_in, SEXP radius,
     SEXP R_weights_in, SEXP edge_or_side, SEXP es_length)
    {
        IndexGraph g;

        if (!isInteger(R_edges_in)) error("R_edges_in should be integer");

        int NE = asInteger(num_edges_in);
        int* edges_in = INTEGER(R_edges_in);

        for (int i = 0; i < NE ; i++, edges_in += 2)
            boost::add_edge(*edges_in, *(edges_in+1), g);

        if ( method == E_LAYOUT_CIRCLE )
        {
            double r = REAL(radius)[0];
            circle_graph_layout(g, get(vertex_position, g), r);
        }
        else if ( method == E_LAYOUT_KKSL )
        {
            property_map<IndexGraph, vertex_position_t>::type 
		    p = get(vertex_position, g);
            property_map<IndexGraph, edge_weight_t>::type 
		    w = get(edge_weight, g);

            int* weight_i = (isReal(R_weights_in)) ? 0 : INTEGER(R_weights_in);
            double* weight_d = (isReal(R_weights_in)) ? REAL(R_weights_in) : 0;

            graph_traits< IndexGraph>::edge_iterator e, e_end;
            for ( tie(e, e_end) = edges(g); e != e_end; ++e )
                w[*e] = weight_i ? (*weight_i++) : (*weight_d++);

            double l = REAL(es_length)[0];// * 50;
            bool e_or_s = LOGICAL(edge_or_side)[0];

            // what "radius" for call to circle_graph_layout?
            circle_graph_layout(g, get(vertex_position, g), l);
            bool ok;
	    if ( e_or_s )
		    ok = kamada_kawai_spring_layout(g, p, w, edge_length(l));
	    else
		    ok = kamada_kawai_spring_layout(g, p, w, side_length(l));
        }

        SEXP anslst, poslst;
        PROTECT(anslst = allocVector(VECSXP,1));
        PROTECT(poslst = allocMatrix(REALSXP,2, num_vertices(g)));

        int i = 0;
        property_map < IndexGraph, vertex_position_t >:: type
        position = get(vertex_position, g);
        graph_traits<IndexGraph>::vertex_iterator vi, vi_end;
        for (tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi)
        {
            REAL(poslst)[i++] = (double)position[*vi].x;
            REAL(poslst)[i++] = (double)position[*vi].y;
        }

        SET_VECTOR_ELT(anslst,0,poslst);
        UNPROTECT(2);
        return(anslst);
    }

    SEXP BGL_circle_layout
    (SEXP num_verts_in, SEXP num_edges_in, SEXP R_edges_in, SEXP radius)
    {
        SEXP anslst, dummy1=0, dummy2=0, dummy3=0;
        anslst = BGL_layout_internal(E_LAYOUT_CIRCLE, num_verts_in, num_edges_in, R_edges_in, radius, dummy1, dummy2, dummy3);

        return(anslst);
    }

    SEXP BGL_kamada_kawai_spring_layout
    (SEXP num_verts_in, SEXP num_edges_in, SEXP R_edges_in,
     SEXP R_weights_in, SEXP edge_or_side, SEXP es_length)
    {
        SEXP anslst, dummy=0;
        anslst = BGL_layout_internal(E_LAYOUT_KKSL, num_verts_in, num_edges_in, R_edges_in, dummy, R_weights_in, edge_or_side, es_length);
        return(anslst);
    }
}

