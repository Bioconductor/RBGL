#include "RBGL.hpp"
#include <boost/pending/disjoint_sets.hpp>
#include <boost/graph/incremental_components.hpp>

extern "C" {
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
}

extern "C"
{
	using namespace boost;

        typedef graph_traits<Graph_ud>::vertex_descriptor Vertex;
        typedef graph_traits<Graph_ud>::vertices_size_type size_type;
	typedef Vertex* Parent;
	typedef size_type* Rank;

	static std::vector<size_type> rank(1);
	static std::vector<Vertex> parent(1);
	static disjoint_sets<Rank, Parent> ds(&rank[0], &parent[0]);
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
	    rank.clear();   rank.resize(NV, 0);
	    parent.clear(); parent.resize(NV, 0);

  	    disjoint_sets<Rank, Parent>  ds1(&rank[0], &parent[0]);
  	    ds = ds1;
	    initialize_incremental_components(g, ds);

	    initialized = true;
	}

	if ( method == E_IC_INCREMENTAL_COMPONENT )
	{
	    incremental_components(g, ds);
	}

	Components components(&parent[0], &parent[0] + parent.size());
	
	int NC = components.size();

        SEXP anslst, conn, outvec;
        PROTECT(anslst = allocVector(VECSXP,NC+1));
	PROTECT(conn = NEW_INTEGER(1));

	INTEGER(conn)[0] = NC;
	SET_VECTOR_ELT(anslst,0,conn);

	int l;
	Components::size_type k;
	Components::value_type::iterator j;

        for (k = 0; k < NC; k++ )
	{
		// count how many vertices in this component: no size()
		l = 0;
		for (j = components[k].begin(); j != components[k].end(); j++ ) 
			l++;
		PROTECT(outvec = allocMatrix(INTSXP, 1, l));

		l = 0;
		for (j = components[k].begin(); j != components[k].end(); j++ )
			INTEGER(outvec)[l++] = *j;

		SET_VECTOR_ELT(anslst,k+1,outvec);
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

