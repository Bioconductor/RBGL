// NOTE:   May 31, 2006
// gcc 4.1 has problem handling this function as it was part of interfaces.cpp;
// initial investigation points to some interplay among the header files;
// the current include-order works, while the following order generates lots
// of errors in compiling:
//       #include "RBGL.hpp"
//       #include <boost/graph/transitive_closure.hpp>
//       #include <boost/graph/graphviz.hpp>
// It's a work-around to separate this function from the rest and include header
// files differently.
// NOTE:   Jun 1, 2006
// The problem seems in 
//       #include <boost/graph/graphviz.hpp>
// since we don't need it, we remove it from RBGL.hpp, everything seems ok.

#include <boost/graph/transitive_closure.hpp>
#include "RBGL.hpp"

extern "C"
{
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

}

