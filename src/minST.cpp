#include "RBGL.hpp"
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>

extern "C"
{

    SEXP BGL_KMST_D( SEXP num_verts_in, SEXP num_edges_in,
                     SEXP R_edges_in, SEXP R_weights_in)
    {
        using namespace boost;

        typedef graph_traits < Graph_dd >::edge_descriptor Edge;
        typedef graph_traits < Graph_dd >::vertex_descriptor Vertex;
        Graph_dd g(num_verts_in, num_edges_in, R_edges_in, R_weights_in);
        property_map < Graph_dd, edge_weight_t >::type weight = get(edge_weight, g);

        std::vector < Edge > spanning_tree;

        kruskal_minimum_spanning_tree(g, std::back_inserter(spanning_tree));

        SEXP ansList, ans, answt;
        PROTECT(ansList = Rf_allocVector(VECSXP,2));
        PROTECT(ans = Rf_allocMatrix(INTSXP,2,spanning_tree.size()));
        PROTECT(answt = Rf_allocMatrix(REALSXP,1,spanning_tree.size()));
        int k = 0, j = 0;

        for (std::vector < Edge >::iterator ei = spanning_tree.begin();
                ei != spanning_tree.end(); ++ei)
        {
            INTEGER(ans)[k++] = source(*ei,g);
            INTEGER(ans)[k++] = target(*ei,g);
            REAL(answt)[j++] = weight[*ei];
        }

        SET_VECTOR_ELT(ansList,0,ans);
        SET_VECTOR_ELT(ansList,1,answt);
        UNPROTECT(3);
        return(ansList);
    } 

    SEXP BGL_KMST_U( SEXP num_verts_in, SEXP num_edges_in,
                     SEXP R_edges_in, SEXP R_weights_in)
    {
        using namespace boost;

        typedef graph_traits < Graph_ud >::edge_descriptor Edge;
        typedef graph_traits < Graph_ud >::vertex_descriptor Vertex;
        Graph_ud g(num_verts_in, num_edges_in, R_edges_in, R_weights_in);
        property_map < Graph_ud, edge_weight_t >::type weight = get(edge_weight, g);

        std::vector < Edge > spanning_tree;

        kruskal_minimum_spanning_tree(g, std::back_inserter(spanning_tree));

        SEXP ansList, ans, answt;
        PROTECT(ansList = Rf_allocVector(VECSXP,2));
        PROTECT(ans = Rf_allocMatrix(INTSXP,2,spanning_tree.size()));
        PROTECT(answt = Rf_allocMatrix(REALSXP,1,spanning_tree.size()));

        int k = 0, j = 0;
        for (std::vector < Edge >::iterator ei = spanning_tree.begin();
                ei != spanning_tree.end(); ++ei)
        {
            INTEGER(ans)[k++] = source(*ei,g);
            INTEGER(ans)[k++] = target(*ei,g);
            REAL(answt)[j++] = weight[*ei];
        }

        SET_VECTOR_ELT(ansList,0,ans);
        SET_VECTOR_ELT(ansList,1,answt);
        UNPROTECT(3);
        return(ansList);
    } 

    SEXP BGL_PRIM_U( SEXP num_verts_in, SEXP num_edges_in,
                     SEXP R_edges_in, SEXP R_weights_in)
    {
        using namespace boost;

        typedef graph_traits < Graph_ud >::edge_descriptor Edge;
        typedef graph_traits < Graph_ud >::vertex_descriptor Vertex;

        Graph_ud g(num_verts_in, num_edges_in, R_edges_in, R_weights_in);

	int NV = Rf_asInteger(num_verts_in);
        std::vector <Vertex> parent(NV);

        prim_minimum_spanning_tree(g, &parent[0]);

        property_map<Graph_ud, edge_weight_t>::type weight = get(edge_weight, g);

        SEXP ansList, ans, answt;
        PROTECT(ansList = Rf_allocVector(VECSXP,2));
        PROTECT(ans = Rf_allocMatrix(INTSXP,2,NV));
        PROTECT(answt = Rf_allocMatrix(REALSXP,1,NV));

        int k = 0, j = 0;
        for (unsigned int v = 0; v < num_vertices(g); ++v)
	{
            INTEGER(ans)[k++] = parent[v];
            INTEGER(ans)[k++] = v;
            REAL(answt)[j++] = ( parent[v] == v) ? 
			 0 : get(weight, edge(parent[v], v, g).first);
	}

        SET_VECTOR_ELT(ansList,0,ans);
        SET_VECTOR_ELT(ansList,1,answt);
        UNPROTECT(3);
        return(ansList);
    }

}

