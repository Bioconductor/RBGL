#include "RBGL.hpp"
#include <boost/graph/isomorphism.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/bandwidth.hpp>
#include <boost/graph/profile.hpp>
#include <boost/graph/wavefront.hpp>
#include <boost/graph/minimum_degree_ordering.hpp>
#include <boost/graph/sloan_ordering.hpp>

extern "C" {
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
}

extern "C"
{
	// These are more easily done in R alone, no point to use BGL
	// copy_graph
	// transpose_graph

	SEXP BGL_isomorphism(
		SEXP num_verts_in1, SEXP num_edges_in1, SEXP R_edges_in1,  
		SEXP num_verts_in2, SEXP num_edges_in2, SEXP R_edges_in2
			 )
	{
		using namespace boost;

		bool r = FALSE;

		const int NV1 = asInteger(num_verts_in1);
		const int NV2 = asInteger(num_verts_in2);
		const int NE1 = asInteger(num_edges_in1);
		const int NE2 = asInteger(num_edges_in2);

		if ( NV1 == NV2 )
		{
		typedef adjacency_list < vecS, listS, undirectedS,
			property < vertex_index_t, int> > VEGraph;

		VEGraph g1(NV1), g2(NV2);

		std::vector<graph_traits<VEGraph>::vertex_descriptor> v1(NV1), v2(NV2);
		property_map<VEGraph, vertex_index_t>::type
			v1_index_map = get(vertex_index, g1),
			v2_index_map = get(vertex_index, g2);

		graph_traits<VEGraph>::vertex_iterator i, end;
		int id = 0;
		for (tie(i, end) = vertices(g1); i != end; ++i, ++id) 
		{
		      put(v1_index_map, *i, id);
		      v1[id] = *i;
		}
		id = 0;
		for (tie(i, end) = vertices(g2); i != end; ++i, ++id) 
		{
		      put(v2_index_map, *i, id);
		      v2[id] = *i;
		}
		
		int* edges_in = INTEGER(R_edges_in1);
		for ( int i = 0; i < NE1; i++, edges_in += 2 )
			add_edge(v1[*edges_in], v1[*(edges_in+1)], g1);
		
		edges_in = INTEGER(R_edges_in2);
		for ( int i = 0; i < NE2; i++, edges_in += 2 )
			add_edge(v2[*edges_in], v2[*(edges_in+1)], g2);
		
		std::vector<graph_traits<VEGraph>::vertex_descriptor> f(NV1);

		r = isomorphism(g1, g2, isomorphism_map 
			(make_iterator_property_map(f.begin(), v1_index_map, f[0])));
		}

		SEXP ansList, conn;
		PROTECT(ansList = allocVector(VECSXP,1));
		PROTECT(conn = NEW_LOGICAL(1));

		LOGICAL(conn)[0] = r;

		SET_VECTOR_ELT(ansList,0,conn);
		UNPROTECT(2);
		return(ansList);
	}

	SEXP BGL_cuthill_mckee_ordering(SEXP num_verts_in, SEXP num_edges_in,
                         SEXP R_edges_in)
	{
		using namespace boost;

		typedef graph_traits<Graph_ud>::vertex_descriptor Vertex;
		typedef graph_traits<Graph_ud>::vertices_size_type size_type;

		const int N = asInteger(num_verts_in);
		std::vector<Vertex>    inv_perm(N);
		std::vector<size_type> perm(N);

		Graph_ud g(num_verts_in, num_edges_in, R_edges_in);

		cuthill_mckee_ordering(g, inv_perm.rbegin(), get(vertex_color, g), make_degree_map(g));

		SEXP ansList, invpermList, robw, rbw;
		PROTECT(ansList = allocVector(VECSXP,3));
		PROTECT(invpermList = allocVector(INTSXP,N));
		PROTECT(robw = NEW_INTEGER(1));
		PROTECT(rbw = NEW_INTEGER(1));

		property_map<Graph_ud, vertex_index_t>::type index_map = get(vertex_index, g);
		std::vector<Vertex>::const_iterator i;
		int j = 0;

		for ( i = inv_perm.begin(); i != inv_perm.end(); i++ )
			INTEGER(invpermList)[j++] = index_map[*i];
		for ( size_type c = 0; c != inv_perm.size(); ++c )
			perm[index_map[inv_perm[c]]] = c;

		INTEGER(robw)[0] = bandwidth(g);
		INTEGER(rbw)[0] = bandwidth(g, make_iterator_property_map(&perm[0], index_map, perm[0]));

		SET_VECTOR_ELT(ansList,0,invpermList);
		SET_VECTOR_ELT(ansList,1,robw);
		SET_VECTOR_ELT(ansList,2,rbw);

		UNPROTECT(4);
		return(ansList);
	}

	SEXP BGL_min_degree_ordering(SEXP num_verts_in, SEXP num_edges_in,
                         SEXP R_edges_in,   SEXP R_delta)
	{
		using namespace boost;

		int delta = asInteger(R_delta);
		const int NV = asInteger(num_verts_in);
		typedef graph_traits<Graph_dd>::vertex_descriptor Vertex;
		Graph_dd g(num_verts_in, num_edges_in, R_edges_in);
		std::vector<int> inverse_perm(NV, 0);
		std::vector<int> perm(NV, 0);
		std::vector<int> degree(NV, 0);
		std::vector<int> supernode_sizes(NV, 1);
		property_map<Graph_dd, vertex_index_t>::type id = get(vertex_index, g);

		minimum_degree_ordering(g, 
			make_iterator_property_map(&degree[0], id, degree[0]), 
			&inverse_perm[0], 
			&perm[0], 
			make_iterator_property_map(&supernode_sizes[0], id, supernode_sizes[0]), 
			delta, id);

		SEXP ansList, invpermList, permList;
		PROTECT(ansList = allocVector(VECSXP,2));
		PROTECT(invpermList = allocVector(INTSXP,NV));
		PROTECT(permList = allocVector(INTSXP,NV));

		std::vector<int>::const_iterator i;
		int j = 0;

		for ( i = inverse_perm.begin(); i != inverse_perm.end(); i++ )
			INTEGER(invpermList)[j++] = inverse_perm[*i];

		j = 0;
		for ( i = perm.begin(); i != perm.end(); i++ )
			INTEGER(permList)[j++] = perm[*i];

		SET_VECTOR_ELT(ansList,0,invpermList);
		SET_VECTOR_ELT(ansList,1,permList);
		UNPROTECT(3);
		return(ansList);
	}

	SEXP BGL_sloan_ordering(SEXP num_verts_in, SEXP num_edges_in,
                     SEXP R_edges_in, SEXP R_weights_in, SEXP R_W1, SEXP R_W2 )
	{
		using namespace boost;

		const int NV = asInteger(num_verts_in);
		const int NE = asInteger(num_edges_in);
		const int W1 = asInteger(R_W1);
		const int W2 = asInteger(R_W2);

		typedef adjacency_list< setS, vecS, undirectedS, 
		    property< vertex_color_t, default_color_type,
		    property< vertex_degree_t, int,
		    property< vertex_priority_t, double > > > > IGraph;

		typedef graph_traits<IGraph>::vertex_descriptor Vertex;
		typedef graph_traits<IGraph>::vertices_size_type size_type;
		typedef std::pair<std::size_t, std::size_t> Pair;

		IGraph g(NV);

		int* edges_in = INTEGER(R_edges_in);
		for ( int i = 0; i < NE; i++, edges_in += 2 )
			add_edge(*edges_in, *(edges_in+1), g);
		
		property_map<IGraph, vertex_index_t>::type index_map = get(vertex_index, g);

		std::vector<Vertex>    sloan_order(NV);
		std::vector<size_type> perm(NV);

		sloan_ordering(g, sloan_order.begin(), get(vertex_color, g), 
			make_degree_map(g), get(vertex_priority, g), W1, W2);

		SEXP ansList, sList, rbw, rpf, rmw, raw, rrw;
		PROTECT(ansList = allocVector(VECSXP,6));
		PROTECT(sList = allocVector(INTSXP,NV));
		PROTECT(rbw = NEW_INTEGER(1));
		PROTECT(rpf = NEW_INTEGER(1));
		PROTECT(rmw = NEW_INTEGER(1));
		PROTECT(raw = NEW_INTEGER(1));
		PROTECT(rrw = NEW_INTEGER(1));

		std::vector<Vertex>::const_iterator i;
		int j = 0;

		for ( i = sloan_order.begin(); i != sloan_order.end(); i++ )
			INTEGER(sList)[j++] = index_map[*i];

		for ( size_type c = 0; c != sloan_order.size(); ++c )
			perm[index_map[sloan_order[c]]] = c;

		INTEGER(rbw)[0] = bandwidth(g, make_iterator_property_map(&perm[0], index_map, perm[0]));
		INTEGER(rpf)[0] = profile(g, make_iterator_property_map(&perm[0], index_map, perm[0]));
		INTEGER(rmw)[0] = max_wavefront(g, make_iterator_property_map(&perm[0], index_map, perm[0]));
		INTEGER(raw)[0] = aver_wavefront(g, make_iterator_property_map(&perm[0], index_map, perm[0]));
		INTEGER(rrw)[0] = rms_wavefront(g, make_iterator_property_map(&perm[0], index_map, perm[0]));

		SET_VECTOR_ELT(ansList,0,sList);
		SET_VECTOR_ELT(ansList,1,rbw);
		SET_VECTOR_ELT(ansList,2,rpf);
		SET_VECTOR_ELT(ansList,3,rmw);
		SET_VECTOR_ELT(ansList,4,raw);
		SET_VECTOR_ELT(ansList,5,rrw);

		UNPROTECT(7);
		return(ansList);
	}

}

