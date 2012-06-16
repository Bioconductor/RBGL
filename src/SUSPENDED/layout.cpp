#include "RBGL.hpp"
#include <boost/graph/circle_layout.hpp>
#include <boost/graph/kamada_kawai_spring_layout.hpp>
#include <boost/graph/random_layout.hpp>

// Work-around:
// the original .hpp simply does k*k/d (and the like) without checking d==0,
// this is particularly troublesome when the given graph is disconnected
#include "fruchterman_reingold.hpp"
//#include <boost/graph/fruchterman_reingold.hpp>

#include <boost/graph/gursoy_atun_layout.hpp>

#include <boost/graph/simple_point.hpp>
#include <boost/random/linear_congruential.hpp>

using namespace boost;

enum vertex_position_t { vertex_position };
namespace boost { BOOST_INSTALL_PROPERTY(vertex, position); }

struct point { double x, y; };

typedef enum { E_LAYOUT_CIRCLE, 
               E_LAYOUT_KKSL, 
               E_LAYOUT_RANDOM, 
               E_LAYOUT_FRFD, 
               E_LAYOUT_GA } E_LAYOUT_METHOD;

extern "C"
{
    typedef adjacency_list<vecS, vecS, undirectedS,
    		// vertex properties
    		property<vertex_index_t, int,
    		property<vertex_position_t, simple_point<double> > >,
    		// edge properties
    		property<edge_weight_t, double> >
    		IndexGraph;

    SEXP BGL_layout_internal (E_LAYOUT_METHOD method, 
         SEXP num_verts_in, SEXP num_edges_in, SEXP R_edges_in, 
         SEXP radius,
         SEXP R_weights_in, SEXP edge_or_side, SEXP es_length,
         SEXP R_minX, SEXP R_maxX, SEXP R_minY, SEXP R_maxY,
         SEXP R_width, SEXP R_height
         )
    {
        if (!isInteger(R_edges_in)) error("R_edges_in should be integer");

        int NV = asInteger(num_verts_in);
        int NE = asInteger(num_edges_in);
        int* edges_in = INTEGER(R_edges_in);

        IndexGraph g(NV);

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
        else if ( method == E_LAYOUT_RANDOM)
        {
            double minX = REAL(R_minX)[0];
            double maxX = REAL(R_maxX)[0];
            double minY = REAL(R_minY)[0];
            double maxY = REAL(R_maxY)[0];
    
            property_map<IndexGraph, vertex_position_t>::type 
		    p = get(vertex_position, g);

            minstd_rand gen;
            random_graph_layout(g, p, minX, maxX, minY, maxY, gen);
        }
        else if ( method == E_LAYOUT_FRFD)
        {
            double w = REAL(R_width)[0];
            double h = REAL(R_height)[0];
            
            property_map<IndexGraph, vertex_position_t>::type 
		    p = get(vertex_position, g);

            minstd_rand gen;
            random_graph_layout(g, p, -w/2, w/2, -h/2, h/2, gen);
            fruchterman_reingold_force_directed_layout(g, p, w, h);
        }
        else if ( method == E_LAYOUT_GA )
        {
            double r = REAL(radius)[0];
            double w = REAL(R_width)[0];
            double h = REAL(R_height)[0];

            property_map<IndexGraph, vertex_position_t>::type 
		    p = get(vertex_position, g);

            // TODO: fill in correct codes
            //square_topology < > lst;
            //gursoy_atun_layout(g, lst, p);
        }
        else
        {
        }

        SEXP poslst;
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

        UNPROTECT(1);
        return(poslst);
    }

    SEXP BGL_circle_layout
    (SEXP num_verts_in, SEXP num_edges_in, SEXP R_edges_in, SEXP radius)
    {
        SEXP dummy_weights=0, dummy_e_or_s=0, dummy_es_length=0;
        SEXP dummy_minX=0, dummy_maxX=0, dummy_minY=0, dummy_maxY=0;
        SEXP dummy_width=0, dummy_height=0;

        SEXP anslst; 
        anslst = BGL_layout_internal(E_LAYOUT_CIRCLE, 
                 num_verts_in, num_edges_in, R_edges_in, 
                 radius, 
                 dummy_weights, dummy_e_or_s, dummy_es_length,
                 dummy_minX, dummy_maxX, dummy_minY, dummy_maxY,
                 dummy_width, dummy_height);

        return(anslst);
    }

    SEXP BGL_kamada_kawai_spring_layout
    (SEXP num_verts_in, SEXP num_edges_in, SEXP R_edges_in,
     SEXP R_weights_in, SEXP edge_or_side, SEXP es_length)
    {
        SEXP dummy_radius=0;
        SEXP dummy_minX=0, dummy_maxX=0, dummy_minY=0, dummy_maxY=0;
        SEXP dummy_width=0, dummy_height=0;

        SEXP anslst; 
        anslst = BGL_layout_internal(E_LAYOUT_KKSL, 
                 num_verts_in, num_edges_in, R_edges_in, 
                 dummy_radius, 
                 R_weights_in, edge_or_side, es_length,
                 dummy_minX, dummy_maxX, dummy_minY, dummy_maxY,
                 dummy_width, dummy_height);
        return(anslst);
    }

    SEXP BGL_random_layout
    (SEXP num_verts_in, SEXP num_edges_in, SEXP R_edges_in, 
     SEXP R_minX, SEXP R_maxX, SEXP R_minY, SEXP R_maxY)
    {
        SEXP dummy_radius=0; 
        SEXP dummy_weights=0, dummy_e_or_s=0, dummy_es_length=0;
        SEXP dummy_width=0, dummy_height=0;

        SEXP anslst; 
        anslst = BGL_layout_internal(E_LAYOUT_RANDOM, 
                 num_verts_in, num_edges_in, R_edges_in, 
                 dummy_radius, 
                 dummy_weights, dummy_e_or_s, dummy_es_length,
                 R_minX, R_maxX, R_minY, R_maxY,
                 dummy_width, dummy_height);
        return(anslst);
    }

    SEXP BGL_FRFD_layout
    (SEXP num_verts_in, SEXP num_edges_in, SEXP R_edges_in, 
     SEXP R_width, SEXP R_height)
    {
        SEXP dummy_radius=0;
        SEXP dummy_weights=0, dummy_e_or_s=0, dummy_es_length=0;
        SEXP dummy_minX=0, dummy_maxX=0, dummy_minY=0, dummy_maxY=0;

        SEXP anslst;
        anslst = BGL_layout_internal(E_LAYOUT_FRFD,
                 num_verts_in, num_edges_in, R_edges_in,
                 dummy_radius,
                 dummy_weights, dummy_e_or_s, dummy_es_length,
                 dummy_minX, dummy_maxX, dummy_minY, dummy_maxY,
                 R_width, R_height);
        return(anslst);
    }

    SEXP BGL_gursov_atun_layout
    (SEXP num_verts_in, SEXP num_edges_in, SEXP R_edges_in, 
     SEXP R_width, SEXP R_height, SEXP R_radius)
    {
        SEXP dummy_weights=0, dummy_e_or_s=0, dummy_es_length=0;
        SEXP dummy_minX=0, dummy_maxX=0, dummy_minY=0, dummy_maxY=0;

        SEXP anslst;
        anslst = BGL_layout_internal(E_LAYOUT_FRFD,
                 num_verts_in, num_edges_in, R_edges_in,
                 R_radius,
                 dummy_weights, dummy_e_or_s, dummy_es_length,
                 dummy_minX, dummy_maxX, dummy_minY, dummy_maxY,
                 R_width, R_height);
        return(anslst);
    }

}

