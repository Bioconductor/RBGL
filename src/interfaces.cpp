#include "RBGL2.h"

extern "C" {
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
}

extern "C" 
{
	SEXP BGL_tsort_D(SEXP num_verts_in, SEXP num_edges_in, SEXP R_edges_in)
	{
	int NV = INTEGER(num_verts_in)[0];
	setupGraphTypes
	setTraits( Graph_dd )
	setUnweightedDoubleEdges( Graph_dd )
	   
	typedef list<Vertex> tsOrder;
	tsOrder tsord;
	SEXP tsout;
	
	PROTECT(tsout = NEW_NUMERIC(NV));
	
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
	      for (int j = 0 ; j < NV; j++)
	          REAL(tsout)[j] = 0.0;
	      }
	UNPROTECT(1);
	
	return(tsout);
	      
	} // end BGLtsortD


	SEXP BGL_KMST_D( SEXP num_verts_in, SEXP num_edges_in, 
	    SEXP R_edges_in, SEXP R_weights_in)
	{
	using namespace boost;
	
	Rprintf("warning: directed graph supplied; directions ignored.\n");
	
	setupGraphTypes
	setTraits( Graph_dd )
	setWeightedDoubleEdges( Graph_dd )
	
	std::vector < Edge > spanning_tree;
	
	kruskal_minimum_spanning_tree(g, std::back_inserter(spanning_tree));
	
	SEXP ansList;
	PROTECT(ansList = allocVector(VECSXP,2));
	SEXP ans;
	PROTECT(ans = allocMatrix(INTSXP,2,spanning_tree.size()));
	SEXP answt;
	PROTECT(answt = allocMatrix(REALSXP,1,spanning_tree.size()));
	int k = 0;
	int j = 0;
	
	for (std::vector < Edge >::iterator ei = spanning_tree.begin();
		ei != spanning_tree.end(); ++ei) 
		{
	   		INTEGER(ans)[k] = source(*ei,g);
		   	INTEGER(ans)[k+1] = target(*ei,g);
	   		REAL(answt)[j] = weight[*ei];
	   		k = k + 2;
	   		j = j + 1;
	  	}
	
	SET_VECTOR_ELT(ansList,0,ans);
	SET_VECTOR_ELT(ansList,1,answt);
	UNPROTECT(3);
	return(ansList);
	}
}
