#include "RBGL.hpp"
#include "Basic2DMatrix.hpp"
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/edge_connectivity.hpp>
#include <boost/graph/transitive_closure.hpp>
#include <boost/graph/biconnected_components.hpp>
#include <boost/graph/sequential_vertex_coloring.hpp>

/* need a template with C++ linkage for BFS */
/* adapted from Siek's bfs-example.cpp */

template < typename TimeMap > class bfs_time_visitor
    : public boost::default_bfs_visitor {
    typedef typename boost::property_traits < TimeMap >::value_type T;
public:
    bfs_time_visitor(TimeMap tmap, T & t):m_timemap(tmap), m_time(t) { }
    template < typename Vertex, typename Graph >
    void discover_vertex(Vertex u, const Graph & g) const
    {
        put(m_timemap, u, m_time++);
    }
    TimeMap m_timemap;
    T & m_time;
};

template < typename TimeMap > class dfs_time_visitor
    : public boost::default_dfs_visitor {
    typedef typename boost::property_traits < TimeMap >::value_type T;
public:
    dfs_time_visitor(TimeMap dmap, TimeMap fmap, T & t)
            :  m_dtimemap(dmap), m_ftimemap(fmap), m_time(t) {
    }
    template < typename Vertex, typename Graph >
    void discover_vertex(Vertex u, const Graph & g) const
    {
        put(m_dtimemap, u, m_time++);
    }
    template < typename Vertex, typename Graph >
    void finish_vertex(Vertex u, const Graph & g) const
    {
        put(m_ftimemap, u, m_time++);
    }
    TimeMap m_dtimemap;
    TimeMap m_ftimemap;
    T & m_time;
};

extern "C"
{
    SEXP BGL_tsort_D(SEXP num_verts_in, SEXP num_edges_in, SEXP R_edges_in)
    {
        // tsortbCG -- for bioConductor graph objects

        using namespace boost;
        typedef graph_traits < Graph_dd >::edge_descriptor Edge;
        typedef graph_traits < Graph_dd >::vertex_descriptor Vertex;
        Graph_dd g(num_verts_in, num_edges_in, R_edges_in);

        typedef property_map<Graph_dd, vertex_color_t>::type Color;
        graph_traits<Graph_dd>::vertex_iterator viter, viter_end;

        typedef std::list<Vertex> tsOrder;
        tsOrder tsord;
        SEXP tsout;

        PROTECT(tsout = NEW_NUMERIC(INTEGER(num_verts_in)[0]));

        try {
            topological_sort(g, std::front_inserter(tsord));

            int j = 0;
            for (tsOrder::iterator i = tsord.begin();
                    i != tsord.end(); ++i)
            {
                REAL(tsout)[j] = (double) *i;
                j++;
            }
        }
        catch ( not_a_dag )
        {
            Rprintf("not a dag, returning zeroes\n");
            for (int j = 0 ; j < INTEGER(num_verts_in)[0]; j++)
                REAL(tsout)[j] = 0.0;
        }
        UNPROTECT(1);

        return(tsout);
    } // end BGL_tsort_D


    SEXP BGL_bfs_D(SEXP num_verts_in, SEXP num_edges_in,
                   SEXP R_edges_in, SEXP R_weights_in, SEXP init_ind)
    {
        using namespace boost;

        typedef graph_traits < Graph_dd >::edge_descriptor Edge;
        typedef graph_traits < Graph_dd >::vertex_descriptor Vertex;
        Graph_dd g(num_verts_in, num_edges_in, R_edges_in, R_weights_in);

        typedef graph_traits < Graph_dd >::vertices_size_type size_type;

        const int N = INTEGER(num_verts_in)[0];
        // Typedefs
        typedef size_type* Iiter;

        // discover time properties
        std::vector < size_type > dtime(num_vertices(g));

        size_type time = 0;
        bfs_time_visitor < size_type * >vis(&dtime[0], time);
        breadth_first_search(g, vertex((int)INTEGER(init_ind)[0], g), visitor(vis));


        // use std::sort to order the vertices by their discover time
        std::vector < size_type > discover_order(N);
        integer_range < size_type > r(0, N);
        std::copy(r.begin(), r.end(), discover_order.begin());
        std::sort(discover_order.begin(), discover_order.end(),
                  indirect_cmp < Iiter, std::less < size_type > >(&dtime[0]));

        SEXP disc;
        PROTECT(disc = allocVector(INTSXP,N));

        int i;
        for (i = 0; i < N; ++i)
        {
            INTEGER(disc)[i] = discover_order[i];
        }

        UNPROTECT(1);
        return(disc);
    }


    SEXP BGL_dfs_D(SEXP num_verts_in, SEXP num_edges_in, SEXP R_edges_in,
                   SEXP R_weights_in)
    {
        using namespace boost;

        typedef graph_traits < Graph_dd >::edge_descriptor Edge;
        typedef graph_traits < Graph_dd >::vertex_descriptor Vertex;
        Graph_dd g(num_verts_in, num_edges_in, R_edges_in, R_weights_in);

        typedef graph_traits < Graph_dd >::vertices_size_type size_type;

        const int N = INTEGER(num_verts_in)[0];
        // Typedefs
        typedef size_type* Iiter;

        // discover time and finish time properties
        std::vector < size_type > dtime(num_vertices(g));
        std::vector < size_type > ftime(num_vertices(g));
        size_type t = 0;
        dfs_time_visitor < size_type * >vis(&dtime[0], &ftime[0], t);

        depth_first_search(g, visitor(vis));

        // use std::sort to order the vertices by their discover time
        std::vector < size_type > discover_order(N);
        integer_range < size_type > r(0, N);
        std::copy(r.begin(), r.end(), discover_order.begin());
        std::sort(discover_order.begin(), discover_order.end(),
                  indirect_cmp < Iiter, std::less < size_type > >(&dtime[0]));
        std::vector < size_type > finish_order(N);
        std::copy(r.begin(), r.end(), finish_order.begin());
        std::sort(finish_order.begin(), finish_order.end(),
                  indirect_cmp < Iiter, std::less < size_type > >(&ftime[0]));

        SEXP ansList;
        PROTECT(ansList = allocVector(VECSXP,2));
        SEXP disc;
        PROTECT(disc = allocVector(INTSXP,N));
        SEXP fin;
        PROTECT(fin = allocVector(INTSXP,N));

        int i;
        for (i = 0; i < N; ++i)
        {
            INTEGER(disc)[i] = discover_order[i];
            INTEGER(fin)[i] = finish_order[i];
        }

        SET_VECTOR_ELT(ansList,0,disc);
        SET_VECTOR_ELT(ansList,1,fin);
        UNPROTECT(3);
        return(ansList);
    }

    SEXP BGL_connected_components_U (SEXP num_verts_in,
                                     SEXP num_edges_in, SEXP R_edges_in, 
				     SEXP R_weights_in )
    {
        using namespace boost;
        SEXP outvec;

        typedef graph_traits < Graph_ud >::edge_descriptor Edge;
        typedef graph_traits < Graph_ud >::vertex_descriptor Vertex;
        Graph_ud g(num_verts_in, num_edges_in, R_edges_in, R_weights_in);

        int nvert = INTEGER(num_verts_in)[0] ;

        std::vector<int> component(num_vertices(g));
        connected_components(g, &component[0]);

        std::vector<int>::size_type k;

        PROTECT(outvec = allocVector(REALSXP,nvert));

        for (k = 0; k < component.size(); k++ )
            REAL(outvec)[k] = component[k];

        UNPROTECT(1);
        return(outvec);
    }

    SEXP BGL_strong_components_D (SEXP num_verts_in,
                                  SEXP num_edges_in, SEXP R_edges_in,
                                  SEXP R_weights_in )
    {
        using namespace boost;
        SEXP outvec;

        typedef graph_traits < Graph_dd >::edge_descriptor Edge;
        typedef graph_traits < Graph_dd >::vertex_descriptor Vertex;
        Graph_dd g(num_verts_in, num_edges_in, R_edges_in, R_weights_in);

        int nvert = INTEGER(num_verts_in)[0] ;

        std::vector<int> component(num_vertices(g));
        strong_components(g, &component[0]);

        std::vector<int>::size_type k;

        PROTECT(outvec = allocVector(REALSXP,nvert));

        for (k = 0; k < component.size(); k++ )
            REAL(outvec)[k] = component[k];

        UNPROTECT(1);
        return(outvec);
    }

    SEXP BGL_biconnected_components_U (SEXP num_verts_in,
                                  SEXP num_edges_in, SEXP R_edges_in,
                                  SEXP R_weights_in )
    {
        using namespace boost;
        SEXP outvec;

        typedef graph_traits < Graph_ud >::edge_descriptor Edge;
        typedef graph_traits < Graph_ud >::vertex_descriptor Vertex;
        Graph_ud g(num_verts_in, num_edges_in, R_edges_in, R_weights_in);

        int ne = INTEGER(num_edges_in)[0] ;

        // this is a bit cheating: use "edge_weight_t" for "edge_component_t"
        property_map < Graph_ud, edge_weight_t >::type 
              component = get(edge_weight, g);
        graph_traits < Graph_ud >::edge_iterator ei, ei_end;
        for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
            component[*ei] = (int)-1;
        int num_comps = biconnected_components(g, component);

        SEXP ansList, nc;
        PROTECT(ansList = allocVector(VECSXP,2));
        PROTECT(nc = NEW_INTEGER(1));
        PROTECT(outvec = allocVector(INTSXP,ne));

        INTEGER(nc)[0] = num_comps;
        int k = 0;
        for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
            INTEGER(outvec)[k++] = (int)component[*ei];

        SET_VECTOR_ELT(ansList,0,nc);
        SET_VECTOR_ELT(ansList,1,outvec);
        UNPROTECT(3);
        return(ansList);
    }

    SEXP BGL_articulation_points_U (SEXP num_verts_in,
                                  SEXP num_edges_in, SEXP R_edges_in,
                                  SEXP R_weights_in )
    {
        using namespace boost;
        SEXP outvec;

        typedef graph_traits < Graph_ud >::edge_descriptor Edge;
        typedef graph_traits < Graph_ud >::vertex_descriptor Vertex;
        Graph_ud g(num_verts_in, num_edges_in, R_edges_in, R_weights_in);

        std::vector<Vertex> art_points;
        articulation_points(g, std::back_inserter(art_points));

        SEXP ansList, nc;
        PROTECT(ansList = allocVector(VECSXP,2));
        PROTECT(nc = NEW_INTEGER(1));
        PROTECT(outvec = allocVector(INTSXP,art_points.size()));

        INTEGER(nc)[0] = art_points.size();

        for (int k = 0; k < art_points.size(); k++ )
            INTEGER(outvec)[k] = art_points[k];

        SET_VECTOR_ELT(ansList,0,nc);
        SET_VECTOR_ELT(ansList,1,outvec);
        UNPROTECT(3);
        return(ansList);
    }

    SEXP BGL_edge_connectivity_U (SEXP num_verts_in,
                                  SEXP num_edges_in, SEXP R_edges_in,
                                  SEXP R_weights_in )
    {
        using namespace boost;
        SEXP ansList, conn, edTmp;

        typedef graph_traits < Graph_ud >::edge_descriptor Edge;
        typedef graph_traits < Graph_ud >::vertex_descriptor Vertex;
        Graph_ud g(num_verts_in, num_edges_in, R_edges_in, R_weights_in);

        typedef graph_traits<Graph_ud>::degree_size_type dst;
        std::vector<Edge> disconnecting_set;
        std::vector<Edge>::iterator ei;
        dst c = edge_connectivity( g, std::back_inserter(disconnecting_set) );

        PROTECT(conn = NEW_NUMERIC(1));
        REAL(conn)[0] = (double)c;

        SEXP eList;
        PROTECT(ansList = allocVector(VECSXP,2));

        PROTECT(eList = allocVector(VECSXP,(int)c));

        SET_VECTOR_ELT(ansList,0,conn);

        int sind = 0;
        for (ei = disconnecting_set.begin(); ei != disconnecting_set.end();
                ++ei)
        {
            PROTECT(edTmp = NEW_NUMERIC(2));
            REAL(edTmp)[0] = (double)source(*ei,g);
            REAL(edTmp)[1] = (double)target(*ei,g);
            SET_VECTOR_ELT(eList,sind,edTmp);
            sind=sind+1;
            UNPROTECT(1);
        }

        SET_VECTOR_ELT(ansList,1,eList);
        UNPROTECT(3);
        return(ansList);
    }

    SEXP BGL_transitive_closure_D (SEXP num_verts_in, 
    		SEXP num_edges_in, SEXP R_edges_in )
    {
    	using namespace boost;
    	
        Graph_dd g(num_verts_in, num_edges_in, R_edges_in );
        Graph_dd TC;
    	
    	transitive_closure(g, TC);

        SEXP ansList, eList, vList;
        PROTECT(ansList = allocVector(VECSXP,2));
        PROTECT(vList = allocMatrix(INTSXP, 1, num_vertices(TC)));
        PROTECT(eList = allocMatrix(INTSXP, 2, num_edges(TC)));

        Graph_dd::vertex_iterator vi, v_end;
        int i = 0;
        for (i = 0, tie(vi, v_end) = vertices(TC); vi != v_end; ++vi)
	{
		INTEGER(vList)[i++] = *vi;
	}

        Graph_dd::edge_iterator ei, e_end;
        for (i = 0, tie(ei, e_end) = edges(TC); ei != e_end; ++ei)
        {
            INTEGER(eList)[i++] = source(*ei, TC);
            INTEGER(eList)[i++] = target(*ei, TC);
        }
    
        SET_VECTOR_ELT(ansList,0,vList);
        SET_VECTOR_ELT(ansList,1,eList);
        UNPROTECT(3);
        return(ansList);
    }

    SEXP BGL_sequential_vertex_coloring (SEXP num_verts_in, 
    		SEXP num_edges_in, SEXP R_edges_in)
    {
        using namespace boost;

        typedef graph_traits < Graph_ud >::edge_descriptor Edge;
        typedef graph_traits < Graph_ud >::vertex_descriptor Vertex;
        typedef graph_traits < Graph_ud >::vertices_size_type Vertex_Size_Type;
        typedef property_map < Graph_ud, vertex_index_t >::const_type vertex_index_map;

        Graph_ud g(num_verts_in, num_edges_in, R_edges_in);

        std::vector<Vertex_Size_Type> color_vec(num_vertices(g));
        iterator_property_map < Vertex_Size_Type*, vertex_index_map >
             color(&color_vec.front(), get(vertex_index, g));
        Vertex_Size_Type n = sequential_vertex_coloring( g, color );

        SEXP ansList, nc, cList;
        PROTECT(ansList = allocVector(VECSXP,2));
        PROTECT(nc = NEW_INTEGER(1));
        PROTECT(cList = allocVector(INTSXP, num_vertices(g)));
        INTEGER(nc)[0] = (int)n;

        Graph_ud::vertex_iterator vi, v_end;
        int i = 0;
        for (tie(vi, v_end) = vertices(g); vi != v_end; ++vi)
        {
            INTEGER(cList)[i++] = color_vec[*vi];
        }
    
        SET_VECTOR_ELT(ansList, 0, nc);
        SET_VECTOR_ELT(ansList, 1, cList);
        UNPROTECT(3);
        return(ansList);
    }

    SEXP BGL_astar_search_D (SEXP num_verts_in, 
    		SEXP num_edges_in, SEXP R_edges_in )
    {
        // TODO: fill in
    	using namespace boost;

        SEXP ansList;
        return(ansList);
    }
}
