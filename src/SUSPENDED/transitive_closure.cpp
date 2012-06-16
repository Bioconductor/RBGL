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

#include "RBGL.hpp"
#include <boost/graph/transitive_closure.hpp>
#include <boost/graph/graph_utility.hpp>


extern "C"
{
    SEXP BGL_transitive_closure_D (SEXP num_verts_in, 
    		SEXP num_edges_in, SEXP R_edges_in )
    {

/*
    	using namespace boost;
    	
        Graph_dd TC;
*/

using namespace boost;
  typedef property < vertex_name_t, char >Name;
  typedef property < vertex_index_t, std::size_t, Name > Index;
  typedef adjacency_list < listS, listS, directedS, Index > graph_t;
  typedef graph_traits < graph_t >::vertex_descriptor vertex_t;
  graph_t G;
  std::vector < vertex_t > verts(4);
  for (int i = 0; i < 4; ++i)
    verts[i] = add_vertex(Index(i, Name('a' + i)), G);
  Graph_dd g(num_verts_in, num_edges_in, R_edges_in );
  add_edge(verts[1], verts[2], G);
  add_edge(verts[1], verts[3], G);
  add_edge(verts[2], verts[1], G);
  add_edge(verts[3], verts[2], G);
  add_edge(verts[3], verts[0], G);

  adjacency_list <> TC;
  transitive_closure(G, TC);
/*
*/
    	
 //   	transitive_closure(g, TC);

/*
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

*/
    }

}

