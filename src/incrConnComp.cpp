#include "RBGL.hpp"
#include <boost/foreach.hpp>
#include <boost/pending/disjoint_sets.hpp>
#include <boost/graph/incremental_components.hpp>

extern "C"
{
	using namespace boost;
	using namespace std;

        typedef graph_traits<Graph_ud>::vertex_descriptor Vertex;
        typedef graph_traits<Graph_ud>::vertices_size_type size_type;
  
  typedef graph_traits<Graph_ud>::vertices_size_type VertexIndex;

	typedef Vertex* Parent;
	typedef size_type* Rank;
	typedef VertexIndex* Rank;

	static vector<size_type> rrank(1);
	static vector<Vertex> parent(1);
	static disjoint_sets<Rank, Parent> ds( &rrank[0], &parent[0]);
	static bool initialized = false;

	typedef component_index<unsigned int> Components;

	typedef enum {
		E_IC_INIT_INCREMENTAL_COMPONENT, 
		E_IC_INCREMENTAL_COMPONENT 
	} E_IC_METHOD;

    static SEXP BGL_incr_comp_internal(SEXP num_verts_in,
                SEXP num_edges_in, SEXP R_edges_in, E_IC_METHOD method)
    {
	using namespace boost;
        int NV = INTEGER(num_verts_in)[0] ;

        Graph_ud g(num_verts_in, num_edges_in, R_edges_in);

	if ( !initialized || method == E_IC_INIT_INCREMENTAL_COMPONENT )
	{
	    rrank.clear();   rrank.resize(NV, 0);
	    parent.clear(); parent.resize(NV, 0);

  	    disjoint_sets<Rank, Parent>  ds1(&rrank[0], &parent[0]);
  	    ds = ds1;
	    initialize_incremental_components(g, ds);

	    initialized = true;
	}

	if ( method == E_IC_INCREMENTAL_COMPONENT )
	{
	    incremental_components(g, ds);
	}

	
        Components components(parent.begin(), parent.end());
	unsigned int NC = components.size();

        SEXP anslst, conn, outvec;
        PROTECT(anslst = allocVector(VECSXP,NC+1));
	PROTECT(conn = NEW_INTEGER(1));

	INTEGER(conn)[0] = NC;
	SET_VECTOR_ELT(anslst,0,conn);

	int l;

        int k = 0;
	BOOST_FOREACH(VertexIndex current_index, components) {
            l = 0;
	    // Iterate through the child vertex indices for [current_index]
	    BOOST_FOREACH(VertexIndex child_index, components[current_index]) {
                l++;  // determine size of component k
            }
            PROTECT(outvec = allocMatrix(INTSXP, 1, l));  // storage for names of component k
            l = 0;
	    BOOST_FOREACH(VertexIndex child_index, components[current_index]) {
                INTEGER(outvec)[l++] =(int) child_index;
            }
	    SET_VECTOR_ELT(anslst,k+1,outvec);   // use offset into anslst, top element is NC
            k = k+1;  // was l+1 seems to be bug
        }

        UNPROTECT(NC+2);
        return(anslst);
    }
    SEXP BGL_init_incremental_components(SEXP num_verts_in,
		SEXP num_edges_in, SEXP R_edges_in)
    {
        SEXP ans = BGL_incr_comp_internal(num_verts_in,
                num_edges_in, R_edges_in, 
		E_IC_INIT_INCREMENTAL_COMPONENT);
	return(ans);
    }

    SEXP BGL_incremental_components(SEXP num_verts_in,
             SEXP num_edges_in, SEXP R_edges_in)
    {
        SEXP ans = BGL_incr_comp_internal(num_verts_in,
                num_edges_in, R_edges_in, 
		E_IC_INCREMENTAL_COMPONENT);
	return(ans);
    }

    SEXP BGL_same_component(SEXP num_verts_in, 
		    SEXP num_edges_in, SEXP R_edges_in, 
		    SEXP vert_1, SEXP vert_2)
    {
        if ( !initialized )
           error("graph is not prepared to handle incremental components.");

        Graph_ud g(num_verts_in, num_edges_in, R_edges_in);

        int NV = INTEGER(num_verts_in)[0];
	int v1 = INTEGER(vert_1)[0];
	int v2 = INTEGER(vert_2)[0];
	
	bool r = FALSE;

	if ( 0 <= v1 && v1 < NV && 0 <= v2 && v2 < NV )
	    r = same_component(vertex(v1, g), vertex(v2, g), ds);
	
	SEXP conn;
	PROTECT(conn = NEW_LOGICAL(1));
	LOGICAL(conn)[0] = r;
	UNPROTECT(1);
	return(conn);
    }

}

