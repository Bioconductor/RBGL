#include "RBGL.hpp"
#include <boost/graph/betweenness_centrality.hpp>
#include <boost/graph/bc_clustering.hpp>

extern "C" {
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
}

extern "C"
{
	using namespace boost;

   	typedef adjacency_list<vecS, vecS, undirectedS,
    		// vertex properties
    		property<vertex_index_t, int>, 
    		// edge properties
    		property<edge_weight_t, double, 
    		property<edge_centrality_t, double, 
    		property<edge_index_t, int> > > >
    		BCGraph;
	typedef graph_traits<BCGraph>::edge_descriptor Edge;

	SEXP BGL_brandes_betweenness_centrality(SEXP num_verts_in, 
		SEXP num_edges_in, SEXP R_edges_in, SEXP R_weights_in)
        {
    		BCGraph g;

		int NV = asInteger(num_verts_in);
		int NE = asInteger(num_edges_in);
		int* edges_in = INTEGER(R_edges_in);
		int* weights_in = INTEGER(R_weights_in);

		for (int i = 0; i < NE ; i++, edges_in += 2, weights_in++)
		    boost::add_edge(*edges_in, *(edges_in+1), *weights_in, g);

		SEXP anslst, bcvlst, bcelst, rbcvlst, dom;
		PROTECT(anslst = allocVector(VECSXP,4));
		PROTECT(bcvlst = allocMatrix(REALSXP, 1, NV));
		PROTECT(bcelst = allocMatrix(REALSXP, 1, NE));
		PROTECT(rbcvlst = allocMatrix(REALSXP, 1, NV));
		PROTECT(dom = NEW_NUMERIC(1));

		std::vector<double> c_map(NV);
		std::vector<double> e_map(NE);

		brandes_betweenness_centrality(g, 
			make_iterator_property_map(c_map.begin(), get(vertex_index, g), double()),
			make_iterator_property_map(e_map.begin(), get(edge_index, g), double()));

		for ( int v = 0; v < NV; v++ ) REAL(bcvlst)[v] = c_map[v];
		for ( int e = 0; e < NE; e++ ) REAL(bcelst)[e] = e_map[e];

		relative_betweenness_centrality(g, 
			make_iterator_property_map(c_map.begin(), get(vertex_index, g), double())); 

		for ( int v = 0; v < NV; v++ ) REAL(rbcvlst)[v] = c_map[v];
		
		double dominance = central_point_dominance(g,
			make_iterator_property_map(c_map.begin(), get(vertex_index, g), double())); 

		REAL(dom)[0] = dominance;

		SET_VECTOR_ELT(anslst,0,bcvlst);
		SET_VECTOR_ELT(anslst,1,bcelst);
		SET_VECTOR_ELT(anslst,2,rbcvlst);
		SET_VECTOR_ELT(anslst,3,dom);
		UNPROTECT(5);
		return(anslst);
	}

	class clustering_threshold : public bc_clustering_threshold<double>
	{
		typedef bc_clustering_threshold<double> inherited;
	public:
	clustering_threshold(double threshold, const BCGraph& g, bool normalize)
		: inherited(threshold, g, normalize), iter(1) { }

	bool operator()(double max_centrality, Edge e, const BCGraph& g)
	 {
		  ++iter;
		  return inherited::operator()(max_centrality, e, g);
	 }

	private:
		 unsigned int iter;
	};

	SEXP BGL_betweenness_centrality_clustering (SEXP num_verts_in, 
		SEXP num_edges_in, SEXP R_edges_in, SEXP R_weights_in,
		SEXP R_threshold,  SEXP R_normalize)
        {
    		BCGraph g;

		int NE = asInteger(num_edges_in);
		int* edges_in = INTEGER(R_edges_in);
		int* weights_in = INTEGER(R_weights_in);

		for (int i = 0; i < NE ; i++, edges_in += 2, weights_in++)
		    boost::add_edge(*edges_in, *(edges_in+1), *weights_in, g);

		double threshold = REAL(R_threshold)[0];
		bool normalize = LOGICAL(R_normalize)[0];

		betweenness_centrality_clustering(g,
			clustering_threshold(threshold, g, normalize),
			get(edge_centrality, g));

		SEXP anslst, cnt, bcvlst, bcelst;
		PROTECT(anslst = allocVector(VECSXP,3));
		PROTECT(cnt = NEW_INTEGER(1));
		PROTECT(bcvlst = allocMatrix(INTSXP, 2, num_edges(g)));
		PROTECT(bcelst = allocMatrix(REALSXP, 1, num_edges(g)));

		INTEGER(cnt)[0] = num_edges(g);

		property_map < BCGraph, edge_centrality_t >::type
		ec = get(edge_centrality, g);

		typedef graph_traits<BCGraph>::edge_iterator   edge_iterator;
		edge_iterator ei, e_end;
		int i = 0, j = 0;
		for ( tie(ei, e_end) = edges(g); ei != e_end; ++ei )
		{
			INTEGER(bcvlst)[i++] = source(*ei, g);
			INTEGER(bcvlst)[i++] = target(*ei, g);
			REAL(bcelst)[j++] = ec[*ei];
		}

		SET_VECTOR_ELT(anslst,0,cnt);
		SET_VECTOR_ELT(anslst,1,bcvlst);
		SET_VECTOR_ELT(anslst,2,bcelst);
		UNPROTECT(4);
		return(anslst);
	}
}

