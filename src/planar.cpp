#include "RBGL.hpp"

#include <boost/graph/boyer_myrvold_planar_test.hpp>
#include <boost/graph/planar_face_traversal.hpp>
#include <boost/graph/planar_canonical_ordering.hpp>
#include <boost/graph/chrobak_payne_drawing.hpp>
#include <boost/graph/is_kuratowski_subgraph.hpp>
#include <boost/graph/make_connected.hpp>
#include <boost/graph/make_biconnected_planar.hpp>
#include <boost/graph/make_maximal_planar.hpp>
#include <boost/graph/max_cardinality_matching.hpp>

// NOTE: 
//    this "undef" is to avoid following conflict:
// in boost: boost/tuple/tuple.hpp:get(tuples::cons<HT, TT>& c)
// in R:     Rinternals.h: #define cons                  Rf_cons
#undef cons
#include <boost/graph/is_straight_line_drawing.hpp>

using namespace boost;

// Note: if these are inside planarFaceTraversal, it causes major 
//      hiccups in g++, this could be reproduced in the example codes
/////////////////////////////////////////////////////////////////////
template <typename Vertex, typename Edge>
struct my_output_visitor : public planar_face_traversal_visitor
{
  my_output_visitor(): v_vis(0), f_vis(0) { }
  ~my_output_visitor() 		{ }
  void begin_face() 		{ v_vis.clear(); }
  void end_face() 		{ f_vis.push_back(v_vis); }
  void next_vertex(Vertex v) 	{ v_vis.push_back(v); }

  typedef std::vector< Vertex > v_vis_t;
  v_vis_t v_vis;
  std::vector< v_vis_t > f_vis;
};

template<typename Graph, typename Vertex>
struct my_add_edge_visitor : public default_add_edge_visitor
{
   void visit_vertex_pair(Vertex u, Vertex v, Graph& g)
   {
        add_edge(u, v, g);
//        std::cout << " add edge: " << u << " " << v << std::endl;
        e_vis.push_back(std::make_pair(u, v));
   }

   typedef std::vector< std::pair< Vertex, Vertex> > e_vis_t;
   e_vis_t e_vis;
};

//////////////////////////////////////////////////////////////////////
//a class to hold the coordinates of the straight line embedding
struct coord_t
{
  std::size_t x;
  std::size_t y;
};

//////////////////////////////////////////////////////////////////////

extern "C"
{
    typedef adjacency_list
    	       < vecS, vecS, undirectedS,
      	         property<vertex_index_t, int>,
      	         property<edge_index_t, int>
    	       > planarGraph;

    typedef graph_traits<planarGraph>::vertex_descriptor Vertex;
    typedef graph_traits<planarGraph>::edge_descriptor Edge;

    typedef property_map<planarGraph, vertex_index_t>::type Vertex_Index_t;
    typedef property_map<planarGraph, edge_index_t>::type Edge_Index_t;

    typedef std::vector<Vertex> Vertex_Vec_t;

    typedef std::vector< Edge > Edge_Vec_t;
    typedef std::vector< Edge_Vec_t > embedding_storage_t;
    typedef boost::iterator_property_map
    	    < embedding_storage_t::iterator, Vertex_Index_t > embedding_t;

    typedef std::vector< coord_t > Coord_Vec_t;
    typedef boost::iterator_property_map
   	    < Coord_Vec_t::iterator, Vertex_Index_t > straight_line_drawing_t;

    Edge_Index_t 	e_index;
    Coord_Vec_t 	straight_line_drawing_storage(0);
    embedding_storage_t embedding_storage(0);

    graph_traits<planarGraph>::edges_size_type 	edge_count;
    graph_traits<planarGraph>::edge_iterator 	ei, ei_end;
    graph_traits<planarGraph>::vertex_iterator 	vi, vi_end;

     void initPlanarGraph(planarGraph* g,
			    SEXP num_verts_in,
                            SEXP num_edges_in,
                            SEXP R_edges_in)
     {
        if ( !isInteger(R_edges_in) ) error("R_edges_in should be integer");

        int NE = asInteger(num_edges_in);
        int* edges_in = INTEGER(R_edges_in);

        for (int i = 0; i < NE ; i++, edges_in += 2) 
	{
            boost::add_edge(*edges_in, *(edges_in+1), 1, *g);
        }
     }

     SEXP boyerMyrvoldPlanarityTest(SEXP num_verts_in, 
			SEXP num_edges_in, SEXP R_edges_in )
     {
	planarGraph g;
	initPlanarGraph(&g, num_verts_in, num_edges_in, R_edges_in);

	SEXP ans;
	PROTECT(ans = NEW_INTEGER(1));

	INTEGER(ans)[0] = boyer_myrvold_planarity_test(g);

	UNPROTECT(1);
	return ans;
     }

     SEXP planarFaceTraversal(SEXP num_verts_in, 
                        SEXP num_edges_in, SEXP R_edges_in )
     {
	planarGraph g;
	initPlanarGraph(&g, num_verts_in, num_edges_in, R_edges_in);

	// Initialize the interior edge index
	e_index = get(edge_index, g);
	edge_count = 0;
	for ( tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
	    put(e_index, *ei, edge_count++);

  	embedding_storage.clear(); 
	embedding_storage.resize(num_vertices(g));

  	if ( boyer_myrvold_planarity_test(
		   boyer_myrvold_params::graph = g,
                   boyer_myrvold_params::embedding = &embedding_storage[0]) )
  	{
  	   my_output_visitor<Vertex, Edge> v_vis;
  	   planar_face_traversal(g, &embedding_storage[0], v_vis);

//#if RBGL_DEBUG 
//  	   std::cout << "we get the following: " << std::endl;
//  	   for ( int i = 0; i < v_vis.f_vis.size(); i++ )
//  	   {
//   	      for ( int j = 0; j < v_vis.f_vis[i].size(); j++ )
//      	         std::cout << v_vis.f_vis[i][j] << " ";
//    	      std::cout << std::endl;
//  	   }
//#endif 
           SEXP ans, ansList;
           PROTECT(ansList = allocVector(VECSXP,v_vis.f_vis.size()));

           for ( int i = 0; i < v_vis.f_vis.size(); i++ )
           {
                PROTECT(ans = allocVector(INTSXP, v_vis.f_vis[i].size()));

                for ( int j = 0; j < v_vis.f_vis[i].size(); j++ )
                    INTEGER(ans)[j] = v_vis.f_vis[i][j];

                SET_VECTOR_ELT(ansList,i, ans);
           }

           UNPROTECT(1+v_vis.f_vis.size());
           return ansList;
  	}
  	else
	{
//       	   std::cout << "Input graph is not planar" << std::endl;

           SEXP ans;
           PROTECT(ans = NEW_INTEGER(1));

           INTEGER(ans)[0] = 0;

           UNPROTECT(1);
           return ans;
	}
     }

     SEXP planarCanonicalOrdering(SEXP num_verts_in,
                        SEXP num_edges_in, SEXP R_edges_in )
     {
	planarGraph g;
	initPlanarGraph(&g, num_verts_in, num_edges_in, R_edges_in);

        // Initialize the interior edge index
	e_index = get(edge_index, g);
        edge_count = 0;
        for ( tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
            put(e_index, *ei, edge_count++);

  	embedding_storage.clear(); 
	embedding_storage.resize(num_vertices(g));

  	embedding_t embedding(embedding_storage.begin(), get(vertex_index,g));

        if ( boyer_myrvold_planarity_test(
		   boyer_myrvold_params::graph = g,
                   boyer_myrvold_params::embedding = embedding) )
        {
           my_add_edge_visitor<planarGraph, Vertex> e1_vis;

           make_connected(g, get(vertex_index, g), e1_vis);

           make_biconnected_planar(g, &embedding_storage[0],
                        get(edge_index, g), e1_vis);
	
           my_add_edge_visitor<planarGraph, Vertex> e_vis;
           make_maximal_planar(g, &embedding_storage[0],
                        get(vertex_index, g), get(edge_index, g), e_vis);

	   // input has to be maximal planar graph with 2+ nodes
  	   Vertex_Vec_t ordering;
  	   planar_canonical_ordering(g, embedding, std::back_inserter(ordering));

//#if RBGL_DEBUG
  	   Vertex_Vec_t::iterator oi, oi_end = ordering.end();
//  	   std::cout << "The planar canonical ordering is: ";
//
//  	   for ( oi = ordering.begin(); oi != oi_end; ++oi)
//    	       std::cout << *oi << " ";
//  	   std::cout << std::endl;
//
  	   straight_line_drawing_storage.clear(); straight_line_drawing_storage.resize(num_vertices(g));
  	   straight_line_drawing_t straight_line_drawing
    		   (straight_line_drawing_storage.begin(), get(vertex_index,g));
//
  	   // Compute the straight line drawing
  	   chrobak_payne_straight_line_drawing(g,  
                                      embedding,
                                      ordering.begin(),
                                      ordering.end(),
                                      straight_line_drawing);
//
//  	   std::cout << "The straight line drawing is: " << std::endl;
//  	   for(tie(vi,vi_end) = vertices(g); vi != vi_end; ++vi)
//    	   {
//      	      std::cout << *vi << " -> (" 
//		   << straight_line_drawing[*vi].x << ", " 
//		   << straight_line_drawing[*vi].y << ")" 
//		   << std::endl;
//    	   }
//#endif

	   SEXP ans;
	   PROTECT(ans = allocVector(INTSXP, ordering.size()));

	   for ( int i = 0; i < ordering.size(); i++ )
	      INTEGER(ans)[i] = ordering[i];

	   UNPROTECT(1);
	   return ans;
	}
	else
	{
//	   std::cout << "Input graph is not planar" << std::endl;

	   SEXP ans;
	   PROTECT(ans = NEW_INTEGER(1));

	   INTEGER(ans)[0] = 0;

	   UNPROTECT(1);
	   return ans;
	}
     }

     SEXP chrobakPayneStraightLineDrawing(SEXP num_verts_in,
                        SEXP num_edges_in, SEXP R_edges_in )
     {
	planarGraph g;
	initPlanarGraph(&g, num_verts_in, num_edges_in, R_edges_in);

  	embedding_storage.clear(); embedding_storage.resize(num_vertices(g));
  	embedding_t embedding(embedding_storage.begin(), get(vertex_index,g));

	// input must be maximal planar graph with 3+ vertices
  	if ( boyer_myrvold_planarity_test(
		   boyer_myrvold_params::graph = g,
                   boyer_myrvold_params::embedding = embedding) )
	{
           my_add_edge_visitor<planarGraph, Vertex> e1_vis;

           make_connected(g, get(vertex_index, g), e1_vis);

           make_biconnected_planar(g, &embedding_storage[0],
                        get(edge_index, g), e1_vis);

           my_add_edge_visitor<planarGraph, Vertex> e_vis;
           make_maximal_planar(g, &embedding_storage[0],
                        get(vertex_index, g), get(edge_index, g), e_vis);

  	   std::vector< Vertex > ordering;
  	   planar_canonical_ordering(g, embedding, std::back_inserter(ordering));

  	   // property map to hold the mapping from vertices to coord_t's
  	   straight_line_drawing_storage.clear(); 
	   straight_line_drawing_storage.resize(num_vertices(g));
  	   straight_line_drawing_t straight_line_drawing
        	   (straight_line_drawing_storage.begin(), get(vertex_index,g));

  	   // Compute the straight line drawing
  	   chrobak_payne_straight_line_drawing(g, 
                                      embedding, 
                                      ordering.begin(),
                                      ordering.end(),
                                      straight_line_drawing);

//#if RBGL_DEBUG
//  	   std::cout << "The straight line drawing is: " << std::endl;
//  	   for(tie(vi,vi_end) = vertices(g); vi != vi_end; ++vi)
//    	   {
//      	      std::cout << *vi << " -> (" 
//		<< straight_line_drawing[*vi].x << ", " 
//		<< straight_line_drawing[*vi].y << ")" 
//                << std::endl;
//    	   }
//#endif

	   SEXP ans;
	   PROTECT(ans = allocMatrix(INTSXP, 2, num_vertices(g)));

	   int j = 0;
  	   for(tie(vi,vi_end) = vertices(g); vi != vi_end; ++vi)
	   {
	      INTEGER(ans)[j++] = straight_line_drawing[*vi].x;
	      INTEGER(ans)[j++] = straight_line_drawing[*vi].y;
	   }

	   UNPROTECT(1);
	   return ans;
	}
	else
	{
//	   std::cout << "Input graph is not planar" << std::endl;

	   SEXP ans;
	   PROTECT(ans = NEW_INTEGER(1));

	   INTEGER(ans)[0] = 0;

	   UNPROTECT(1);
	   return ans;
	}
     }

     SEXP isStraightLineDrawing(SEXP num_verts_in,
                        SEXP num_edges_in, SEXP R_edges_in,
			SEXP drawing_in )
     {
	planarGraph g;
	initPlanarGraph(&g, num_verts_in, num_edges_in, R_edges_in);
 
	bool rstl;

  	straight_line_drawing_storage.clear(); 
	straight_line_drawing_storage.resize(num_vertices(g)); 
  	straight_line_drawing_t straight_line_drawing 
		(straight_line_drawing_storage.begin(), get(vertex_index,g) );

	for ( int i = 0, j = 0; i < num_vertices(g); i++ )
	{
	   straight_line_drawing[i].x = INTEGER(drawing_in)[j++];
	   straight_line_drawing[i].y = INTEGER(drawing_in)[j++];
//      	    std::cout << i << " -> (" 
//		<< straight_line_drawing[i].x << ", " 
//		<< straight_line_drawing[i].y << ")"
//                << std::endl;
	}

  	// Verify that the drawing is actually a plane drawing
  	if (is_straight_line_drawing(g, straight_line_drawing))
	{
//    		std::cout << "Is a plane drawing." << std::endl;
		rstl = 1;
	}
  	else
	{
//    		std::cout << "Is not a plane drawing." << std::endl;
		rstl = 0;
	}

	SEXP ans;
	PROTECT(ans = NEW_INTEGER(1));

	INTEGER(ans)[0] = rstl;
	UNPROTECT(1);
	return ans;
     }

     SEXP isKuratowskiSubgraph(SEXP num_verts_in,
                        SEXP num_edges_in, SEXP R_edges_in )
     {
	planarGraph g;
	initPlanarGraph(&g, num_verts_in, num_edges_in, R_edges_in);

  	// Initialize the interior edge index
	e_index = get(edge_index, g);
  	edge_count = 0;
  	for(tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
    		put(e_index, *ei, edge_count++);
  
	bool is_planar = 0, is_kuratowski = 0;

  	Edge_Vec_t kuratowski_edges;
  	if ( boyer_myrvold_planarity_test(
		   boyer_myrvold_params::graph = g,
                   boyer_myrvold_params::kuratowski_subgraph = 
                                       std::back_inserter(kuratowski_edges)) )
	{
	   is_planar = 1;
//    	   std::cout << "Input graph is planar" << std::endl;
	}
  	else
    	{
	   is_planar = 0;
//      	   std::cout << "Input graph is not planar" << std::endl;

//      	   std::cout << "Edges in the Kuratowski subgraph: ";

      	   Edge_Vec_t::iterator ki, ki_end = kuratowski_edges.end();
//      	   for ( ki = kuratowski_edges.begin(); ki != ki_end; ++ki )
//              std::cout << *ki << " "; 
//      	   std::cout << std::endl;
//
//      	   std::cout << "Is a kuratowski subgraph? ";
      	   if (is_kuratowski_subgraph(g, kuratowski_edges.begin(), kuratowski_edges.end()))
	   {
 //             std::cout << "Yes." << std::endl;
	      is_kuratowski = 1;
	   }
      	   else
	   {
//              std::cout << "No." << std::endl;
	      is_kuratowski = 0;
	   }
    	}

	SEXP ans1, ans2, k_vec, ansList;

	PROTECT(ansList = allocVector(VECSXP,3));

	PROTECT(ans1 = NEW_INTEGER(1));
	INTEGER(ans1)[0] = is_planar;
	SET_VECTOR_ELT(ansList,0,ans1);

	PROTECT(ans2 = NEW_INTEGER(1));
	INTEGER(ans2)[0] = is_kuratowski;
	SET_VECTOR_ELT(ansList,1,ans2);

	PROTECT(k_vec = allocMatrix(INTSXP, 2, kuratowski_edges.size()));

	int i, j = 0;
	Edge_Vec_t::iterator ki, ki_end = kuratowski_edges.end();
  	for ( ki = kuratowski_edges.begin(); ki != ki_end; ++ki )
	{
	      INTEGER(k_vec)[j++] = source(*ki, g);
	      INTEGER(k_vec)[j++] = target(*ki, g);
	}

	SET_VECTOR_ELT(ansList,2, k_vec);

	UNPROTECT(4);
	return ansList;
     }

     SEXP makeConnected(SEXP num_verts_in,
                        SEXP num_edges_in, SEXP R_edges_in )
     {
	planarGraph g;
	initPlanarGraph(&g, num_verts_in, num_edges_in, R_edges_in);

	my_add_edge_visitor<planarGraph, Vertex> e_vis;

	make_connected(g, get(vertex_index, g), e_vis);

//  	for ( int i = 0; i < e_vis.e_vis.size(); i++ )
//  	{
//            std::cout << e_vis.e_vis[i].first << " "
//                  << e_vis.e_vis[i].second
//                  << std::endl;
//  	}

	// output is a graph
	SEXP ans;
	PROTECT(ans = allocMatrix(INTSXP, 2, num_edges(g)));

	int i;
  	for( i = 0, tie(ei, ei_end) = edges(g); ei != ei_end; ++ei )
	{
	   INTEGER(ans)[i++] = source(*ei, g);
	   INTEGER(ans)[i++] = target(*ei, g);
	}

	UNPROTECT(1);
	return ans;
     }

     SEXP makeBiconnectedPlanar(SEXP num_verts_in,
                        SEXP num_edges_in, SEXP R_edges_in )
     {
	planarGraph g;
	initPlanarGraph(&g, num_verts_in, num_edges_in, R_edges_in);

  	//Initialize the interior edge index
	e_index = get(edge_index, g);
  	edge_count = 0;
  	for(tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
    		put(e_index, *ei, edge_count++);
  
  	embedding_storage.clear(); 
	embedding_storage.resize(num_vertices(g));

	bool is_planar = 0;

  	if ( boyer_myrvold_planarity_test(
		   boyer_myrvold_params::graph = g,
                   boyer_myrvold_params::embedding = &embedding_storage[0]) )
	{
 //   	   std::cout << "Input graph is planar" << std::endl;

	   is_planar = 1;

	   // input must be connected planar graph
	   my_add_edge_visitor<planarGraph, Vertex> e_vis;

	   make_connected(g, get(vertex_index, g), e_vis);

 	   make_biconnected_planar(g, &embedding_storage[0], get(edge_index, g), e_vis);

// following drawn from printing code below
        }
  	   if (!boyer_myrvold_planarity_test(g))
	   	is_planar = 0;

  //	   for ( int i = 0; i < e_vis.e_vis.size(); i++ )
//  	   {
//              std::cout << e_vis.e_vis[i].first << " "
//                  	<< e_vis.e_vis[i].second
//                  	<< std::endl;
//  	   }

 // 	   if (boyer_myrvold_planarity_test(g))
//    	      std::cout << "Also, the graph is still planar." << std::endl;
//  	   else
//    	      std::cout << "But the graph is not still planar." << std::endl;
//	}
//  	else
//	{
//    	   std::cout << "Input graph is not planar" << std::endl;
//
//	   is_planar = 0;
//	}
  
	// output is a graph
        SEXP ans, b_vec, ansList;

        PROTECT(ansList = allocVector(VECSXP,2));

        PROTECT(ans = NEW_INTEGER(1));
        INTEGER(ans)[0] = is_planar;
        SET_VECTOR_ELT(ansList,0,ans);

	PROTECT(b_vec = allocMatrix(INTSXP, 2, num_edges(g)));

	int i;
  	for( i = 0, tie(ei, ei_end) = edges(g); ei != ei_end; ++ei )
	{
	   INTEGER(b_vec)[i++] = source(*ei, g);
	   INTEGER(b_vec)[i++] = target(*ei, g);
	}

        SET_VECTOR_ELT(ansList,1,b_vec);

	UNPROTECT(3);
	return ansList;
     }

     SEXP makeMaximalPlanar(SEXP num_verts_in,
                        SEXP num_edges_in, SEXP R_edges_in )
     {
	planarGraph g;
	initPlanarGraph(&g, num_verts_in, num_edges_in, R_edges_in);

  	//Initialize the interior edge index
	e_index = get(edge_index, g);
  	edge_count = 0;
  	for(tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
    		put(e_index, *ei, edge_count++);
  
	bool is_planar = 0;

  	embedding_storage.clear(); 
	embedding_storage.resize(num_vertices(g));
  	if ( boyer_myrvold_planarity_test(
		   boyer_myrvold_params::graph = g,
                   boyer_myrvold_params::embedding = &embedding_storage[0]) )
	{
//    	   std::cout << "Input graph is planar" << std::endl;
	   is_planar = 1;

  	   // input to make_maximal_planar must be biconnected
	   my_add_edge_visitor<planarGraph, Vertex> e1_vis;

	   make_connected(g, get(vertex_index, g), e1_vis);

  	   make_biconnected_planar(g, &embedding_storage[0], 
			get(edge_index, g), e1_vis);

  	   // Re-initialize the edge index, since we just added a few edges
  	   edge_count = 0;
  	   for(tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
    	   	put(e_index, *ei, edge_count++);

  	   embedding_storage.clear(); 
	   embedding_storage.resize(num_vertices(g));
  	   if ( boyer_myrvold_planarity_test(boyer_myrvold_params::graph = g,
                           boyer_myrvold_params::embedding = &embedding_storage[0]))
	   {
 //   	      std::cout << "After calling make_biconnected, the graph is still planar" 
//              		<< std::endl;
	   }
  	   else
	   {
	      // should NOT be here
//    	      std::cout << "After calling make_biconnected, the graph is not planar" 
//              		<< std::endl;
	   }

  	   // input must be biconnected w/ 3+ nodes
	   my_add_edge_visitor<planarGraph, Vertex> e_vis;
  	   make_maximal_planar(g, &embedding_storage[0], 
			get(vertex_index, g), get(edge_index, g), e_vis);

	}
  	else
	{
 //   	   std::cout << "Input graph is not planar" << std::endl;
	   is_planar = 0;
	}
  
	// output is a graph
        SEXP ans, b_vec, ansList;

        PROTECT(ansList = allocVector(VECSXP,2));

        PROTECT(ans = NEW_INTEGER(1));
        INTEGER(ans)[0] = is_planar;
        SET_VECTOR_ELT(ansList,0,ans);

        PROTECT(b_vec = allocMatrix(INTSXP, 2, num_edges(g)));

        int i;
        for( i = 0, tie(ei, ei_end) = edges(g); ei != ei_end; ++ei )
        {
           INTEGER(b_vec)[i++] = source(*ei, g);
           INTEGER(b_vec)[i++] = target(*ei, g);
        }

        SET_VECTOR_ELT(ansList,1,b_vec);

        UNPROTECT(3);
        return ansList; 

     }

    SEXP edmondsMaxCardinalityMatching(SEXP num_verts_in, 
				SEXP num_edges_in,SEXP R_edges_in)
    {
	typedef graph_traits<Graph_ui>::vertex_descriptor Vertex;

	Graph_ui g(num_verts_in, num_edges_in, R_edges_in);

	graph_traits<Graph_ui>::vertex_iterator vi, vi_end;

	std::vector<Vertex> mate(num_vertices(g));

	bool is_max = checked_edmonds_maximum_cardinality_matching(g, &mate[0]);

	SEXP ans, b_vec, ansList;

	PROTECT(ansList = allocVector(VECSXP,2));

        PROTECT(ans = NEW_INTEGER(1));
        INTEGER(ans)[0] = is_max;
        SET_VECTOR_ELT(ansList,0,ans);

        int m_cnt;
        for( m_cnt = 0, tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi )
        {
	   if ( mate[*vi] != graph_traits< Graph_ui >::null_vertex() &&
		*vi < mate[*vi] )
	   m_cnt++;
	}

        PROTECT(b_vec = allocMatrix(INTSXP, 2, m_cnt));

	int i;
        for( i = 0, tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi )
        {
	   if ( mate[*vi] != graph_traits< Graph_ui >::null_vertex() &&
		*vi < mate[*vi] )
	   {
           	INTEGER(b_vec)[i++] = *vi;
           	INTEGER(b_vec)[i++] = mate[*vi];
	   }
        }

        SET_VECTOR_ELT(ansList,1,b_vec);

	UNPROTECT(3);
	return ansList;

    }
}

