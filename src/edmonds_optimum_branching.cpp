#include "RBGL.hpp"

#include "edmonds_optimum_branching.hpp"

extern "C"
{  

    SEXP edmondsOptimumBranching(SEXP num_verts_in,
                SEXP num_edges_in, SEXP R_edges_in, SEXP R_weights_in )
    {

        using namespace boost;

        Graph_dd g(num_verts_in, num_edges_in, R_edges_in, R_weights_in);

        property_map<Graph_dd, boost::edge_weight_t>::type 
		weights = get(edge_weight, g);

	property_map<Graph_dd, boost::vertex_index_t>::type 
		vertex_indices = get(vertex_index, g);

   typedef graph_traits<Graph_dd>::vertex_descriptor       Vertex;
   typedef graph_traits<Graph_dd>::edge_descriptor         Edge;

#ifdef RBGL_DEBUG
    std::cout << "This is the graph:\n";
    BOOST_FOREACH (Edge e, edges(g))
    {
        std::cout << "(" << source(e, g) << ", "
                  << target(e, g) << ")\t"
                  << get(weights, e) << "\n";
    }
#endif

	std::vector<Edge> branching;
	edmonds_optimum_branching<true, false, false>(g,
                                  vertex_indices,
                                  weights,
                                  static_cast<Vertex *>(0),
                                  static_cast<Vertex *>(0),
                                  std::back_inserter(branching));

#ifdef RBGL_DEBUG
    std::cout << "This is the maximum branching\n";
    BOOST_FOREACH (Edge e, branching)
    {
        std::cout << "(" << source(e, g) << ", "
                  << target(e, g) << ")\t"
                  << get(weights, e) << "\n";
    }
#endif

        SEXP ansList, ans, answt;
	PROTECT(ansList = Rf_allocVector(VECSXP,2));
        PROTECT(ans = Rf_allocMatrix(INTSXP, 2, branching.size()));
	PROTECT(answt = Rf_allocMatrix(REALSXP,1,branching.size()));

        int k = 0, j = 0;
	BOOST_FOREACH (Edge e, branching)
    	{
	    INTEGER(ans)[k++] = source(e, g);
	    INTEGER(ans)[k++] = target(e, g);
	    REAL(answt)[j++] = get(weights, e);
	}

	SET_VECTOR_ELT(ansList,0,ans);
	SET_VECTOR_ELT(ansList,1,answt);
        UNPROTECT(3);
        return(ansList);
    }

}

