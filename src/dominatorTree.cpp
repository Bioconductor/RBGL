#include "RBGL.hpp"
#include <boost/graph/dominator_tree.hpp>

extern "C"
{
    SEXP BGL_dominator_tree(SEXP num_verts_in, SEXP num_edges_in, 
         SEXP R_edges_in, SEXP start )
    {
    	using namespace std;
    	using namespace boost;

        const int NV = asInteger(num_verts_in);
        const int NE = asInteger(num_edges_in);
        int v = asInteger(start);

        typedef adjacency_list < vecS, listS, bidirectionalS,
                        property < vertex_index_t, int> > VEGraph;

         VEGraph g(NV);

         typedef property_map < VEGraph, vertex_index_t > ::type IndexMap;
         IndexMap v1_index_map = get(vertex_index, g);
         vector<graph_traits<VEGraph>::vertex_descriptor> v1(NV);
         graph_traits<VEGraph>::vertex_iterator vi, v_end;
         int i = 0;
         for (tie(vi, v_end) = vertices(g); vi != v_end; ++vi, ++i)
         {
               put(v1_index_map, *vi, i);
               v1[i] = *vi;
         }

         int* edges_in = INTEGER(R_edges_in);
         for ( i = 0; i < NE; i++, edges_in += 2 )
         {
             add_edge(v1[*edges_in], v1[*(edges_in+1)], g);
         }

        typedef graph_traits < VEGraph > ::vertex_descriptor Vertex;

        vector< Vertex > domTreePredVector =
            vector<Vertex>(num_vertices(g), graph_traits<VEGraph>::null_vertex());
        iterator_property_map< vector < Vertex > ::iterator, IndexMap > 
          domTreePredMap =
             make_iterator_property_map(domTreePredVector.begin(), v1_index_map);

    	lengauer_tarjan_dominator_tree(g, vertex(v, g), domTreePredMap);

        SEXP ansList;
        PROTECT(ansList = allocVector(INTSXP,num_vertices(g)));

        for (tie(vi, v_end) = vertices(g); vi != v_end; ++vi)
	{
            i = get(v1_index_map, *vi);
            if (get(domTreePredMap, *vi) != graph_traits<VEGraph>::null_vertex())
		INTEGER(ansList)[i] = get(v1_index_map, get(domTreePredMap, *vi));
            else
                INTEGER(ansList)[i] = i; // (numeric_limits<int>::max)();
	}

        UNPROTECT(1);
        return(ansList);

    }

}

