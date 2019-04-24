#include "RBGL.hpp"
#include <boost/graph/graph_utility.hpp>


using namespace boost;

typedef adjacency_list<vecS, vecS, undirectedS,
	// vertex properties
        property<vertex_index_t, int, 
    	property<vertex_centrality_t, double> >, 
    	// edge properties
    	property<edge_weight_t, double, 
    	property<edge_centrality_t, double> > >
    	BCGraph;
typedef graph_traits<BCGraph>::edge_descriptor Edge;
typedef graph_traits<BCGraph>::edge_iterator EdgeIterator;
typedef boost::indirect_cmp<boost::adj_list_edge_property_map<boost::undirected_tag, double, double &, unsigned long,boost::property<boost::edge_weight_t, double, boost::property<boost::edge_centrality_t, double, boost::no_property> >, boost::edge_centrality_t>, std::__1::less<double> > EdgeCentralityCompare;
      
EdgeIterator
max_element(EdgeIterator __first, EdgeIterator __last, EdgeCentralityCompare __comp)
{
    if (__first != __last)
    {
        EdgeIterator __i = __first;
        while (++__i != __last)
            if (__comp(*__first, *__i))
                __first = __i;
    }
    return __first;
}

#include <boost/graph/bc_clustering.hpp>
#include <boost/graph/betweenness_centrality.hpp>



extern "C"
{


	SEXP BGL_brandes_betweenness_centrality(SEXP num_verts_in, 
		SEXP num_edges_in, SEXP R_edges_in, SEXP R_weights_in)
        {
    		BCGraph g;

		int NV = Rf_asInteger(num_verts_in);
		int NE = Rf_asInteger(num_edges_in);
		int* edges_in = INTEGER(R_edges_in);
		double* weights_in = REAL(R_weights_in);

		for (int i = 0; i < NE ; i++, edges_in += 2, weights_in++)
		    boost::add_edge(*edges_in, *(edges_in+1), *weights_in, g);

		SEXP anslst, bcvlst, enlst, bcelst, rbcvlst, dom;
		PROTECT(anslst = Rf_allocVector(VECSXP,5));
		PROTECT(bcvlst = Rf_allocMatrix(REALSXP, 1, NV));
		PROTECT(enlst = Rf_allocMatrix(INTSXP, 2, NE));
		PROTECT(bcelst = Rf_allocMatrix(REALSXP, 1, NE));
		PROTECT(rbcvlst = Rf_allocMatrix(REALSXP, 1, NV));
		PROTECT(dom = Rf_allocVector(REALSXP, 1));

		brandes_betweenness_centrality(g, 
                        centrality_map(get(vertex_centrality, g)).
			edge_centrality_map(get(edge_centrality, g)).
			weight_map(get(edge_weight, g)));

                property_map<BCGraph, vertex_centrality_t>::type 
			v_map = get(vertex_centrality, g);
                property_map<BCGraph, edge_centrality_t>::type 
			e_map = get(edge_centrality, g);

                graph_traits < BCGraph>::vertex_iterator vi, v_end;
                graph_traits < BCGraph>::edge_iterator ei, e_end;

                int v = 0, e = 0;

		for ( tie(vi, v_end) = vertices(g); vi != v_end; vi++ ) 
                    REAL(bcvlst)[v++] = v_map[*vi];
		for ( v = 0, tie(ei, e_end) = edges(g); ei != e_end ; ei++ ) 
		{
		    INTEGER(enlst)[v++] = source(*ei, g);
		    INTEGER(enlst)[v++] = target(*ei, g);
                    REAL(bcelst)[e++] = e_map[*ei];
                }

		relative_betweenness_centrality(g, get(vertex_centrality, g));
                v_map = get(vertex_centrality, g);

		for ( v = 0, tie(vi, v_end) = vertices(g); vi != v_end; vi++ ) 
                    REAL(rbcvlst)[v++] = v_map[*vi];
		
		double dominance = central_point_dominance(g,
		                  get(vertex_centrality, g));

		REAL(dom)[0] = dominance;

		SET_VECTOR_ELT(anslst,0,bcvlst);
		SET_VECTOR_ELT(anslst,1,bcelst);
		SET_VECTOR_ELT(anslst,2,rbcvlst);
		SET_VECTOR_ELT(anslst,3,dom);
		SET_VECTOR_ELT(anslst,4,enlst);
		UNPROTECT(6);
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
#if DEBUG
                  std::cout << "Iter: " << iter << " Max Centrality: "
                       << (max_centrality / dividend) << std::endl;
#endif
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

		int NE = Rf_asInteger(num_edges_in);
		int* edges_in = INTEGER(R_edges_in);
		double* weights_in = REAL(R_weights_in);

		for (int i = 0; i < NE ; i++, edges_in += 2, weights_in++)
		    boost::add_edge(*edges_in, *(edges_in+1), *weights_in, g);

		double threshold = REAL(R_threshold)[0];
		bool normalize = LOGICAL(R_normalize)[0];

		betweenness_centrality_clustering(g,
			clustering_threshold(threshold, g, normalize),
			get(edge_centrality, g));

		// betweenness_centrality_clustering(g,
		// 	clustering_threshold(threshold, g, normalize));

		SEXP anslst, cnt, bcvlst, bcelst;
		PROTECT(anslst = Rf_allocVector(VECSXP,3));
		PROTECT(cnt = Rf_allocVector(INTSXP, 1));
		PROTECT(bcvlst = Rf_allocMatrix(INTSXP, 2, num_edges(g)));
		PROTECT(bcelst = Rf_allocMatrix(REALSXP, 1, num_edges(g)));

		INTEGER(cnt)[0] = num_edges(g);

		property_map < BCGraph, edge_centrality_t >::type
		         ec = get(edge_centrality, g);

		typedef graph_traits<BCGraph>::edge_iterator   edge_iterator;
		edge_iterator ei, e_end;

#if DEBUG
                std::cout << " edge centralities: ";
		for ( tie(ei, e_end) = edges(g); ei != e_end; ++ei )
                        std::cout << " " << ec[*ei];
                std::cout << std::endl;
#endif

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

