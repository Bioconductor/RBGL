#include "RBGL.hpp"

extern "C" {
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
}

/* need a template with C++ linkage for BFS */
/* adapted from Siek's bfs-example.cpp */

template < typename TimeMap > class bfs_time_visitor
    : public boost::default_bfs_visitor {
  typedef typename boost::property_traits < TimeMap >::value_type T;
public:
  bfs_time_visitor(TimeMap tmap, T & t):m_timemap(tmap), m_time(t) { }
  template < typename Vertex, typename Graph >
    void discover_vertex(Vertex u, const Graph & g) const
  {
    put(m_timemap, u, m_time++);
  }
  TimeMap m_timemap;
  T & m_time;
};

template < typename TimeMap > class dfs_time_visitor
    : public boost::default_dfs_visitor {
  typedef typename boost::property_traits < TimeMap >::value_type T;
public:
  dfs_time_visitor(TimeMap dmap, TimeMap fmap, T & t)
:  m_dtimemap(dmap), m_ftimemap(fmap), m_time(t) {
  }
  template < typename Vertex, typename Graph >
    void discover_vertex(Vertex u, const Graph & g) const
  {
    put(m_dtimemap, u, m_time++);
  }
  template < typename Vertex, typename Graph >
    void finish_vertex(Vertex u, const Graph & g) const
  {
    put(m_ftimemap, u, m_time++);
  }
  TimeMap m_dtimemap;
  TimeMap m_ftimemap;
  T & m_time;
};

extern "C" 
{
	SEXP BGL_tsort_D(SEXP num_verts_in, SEXP num_edges_in, SEXP R_edges_in)
	{
	// tsortbCG -- for bioConductor graph objects 
	  
    using namespace boost;
    typedef graph_traits < Graph_dd >::edge_descriptor Edge;
    typedef graph_traits < Graph_dd >::vertex_descriptor Vertex;
    Graph_dd g(num_verts_in, num_edges_in, R_edges_in);

    typedef property_map<Graph_dd, vertex_color_t>::type Color;
	graph_traits<Graph_dd>::vertex_iterator viter, viter_end;
	 
	typedef std::list<Vertex> tsOrder;
	tsOrder tsord;
	SEXP tsout;
	
	PROTECT(tsout = NEW_NUMERIC(INTEGER(num_verts_in)[0]));
	
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
	      for (int j = 0 ; j < INTEGER(num_verts_in)[0]; j++)
	          REAL(tsout)[j] = 0.0;
	      }
	UNPROTECT(1);

	return(tsout);
	} // end BGL_tsort_D
 

	SEXP BGL_KMST_D( SEXP num_verts_in, SEXP num_edges_in, 
	    SEXP R_edges_in, SEXP R_weights_in)
	{
	using namespace boost;
	
//	Rprintf("warning: directed graph supplied; directions ignored.\n");
	
    typedef graph_traits < Graph_dd >::edge_descriptor Edge;
    typedef graph_traits < Graph_dd >::vertex_descriptor Vertex;
    Graph_dd g(num_verts_in, num_edges_in, R_edges_in, R_weights_in);
    property_map < Graph_dd, edge_weight_t >::type weight = get(edge_weight, g);
	
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
	} //end BGL_KMST_D


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
	} //end BGL_KMST_U

/*
	SEXP BGL_PRIM_U( SEXP num_verts_in, SEXP num_edges_in, 
	    SEXP R_edges_in, SEXP R_weights_in)
	{
	using namespace boost;
	
    typedef graph_traits < Graph_ud >::edge_descriptor Edge;
    typedef graph_traits < Graph_ud >::vertex_descriptor Vertex;
    Graph_ud g(num_verts_in, num_edges_in, R_edges_in, R_weights_in);
    property_map < Graph_ud, edge_weight_t >::type weight = get(edge_weight, g);
	
	std::vector <Vertex> parent;
	
	prim_minimum_spanning_tree(g, &parent[0]);

        property_map<Graph_ud, edge_weight_t>::type weight = get(edge_weight, g);
        int total_wgt = 0;
        for (int v = 0; v < num_vertices(g); ++v)
         if (parent([v]) != v)
            total_wgt += get(weight, edge(parent[v], v, g).first);
        Rprintf("total is %d\n", total_wgt);
	
	return(R_NilValue);
	}   */

	SEXP BGL_bfs_D(SEXP num_verts_in, SEXP num_edges_in, 
		SEXP R_edges_in, SEXP R_weights_in, SEXP init_ind)
	{
	using namespace boost;
	
    typedef graph_traits < Graph_dd >::edge_descriptor Edge;
    typedef graph_traits < Graph_dd >::vertex_descriptor Vertex;
    Graph_dd g(num_verts_in, num_edges_in, R_edges_in, R_weights_in);
	
	typedef graph_traits < Graph_dd >::vertices_size_type size_type;
	
	const int N = INTEGER(num_verts_in)[0];
	// Typedefs
	typedef size_type* Iiter;
	
	// discover time properties
	std::vector < size_type > dtime(num_vertices(g));
	
	size_type time = 0;
	bfs_time_visitor < size_type * >vis(&dtime[0], time);
	breadth_first_search(g, vertex((int)INTEGER(init_ind)[0], g), visitor(vis));
	
	
	// use std::sort to order the vertices by their discover time
	std::vector < size_type > discover_order(N);
	integer_range < size_type > r(0, N);
	std::copy(r.begin(), r.end(), discover_order.begin());
	std::sort(discover_order.begin(), discover_order.end(),
	            indirect_cmp < Iiter, std::less < size_type > >(&dtime[0]));
	
	SEXP disc;
	PROTECT(disc = allocVector(INTSXP,N));
	
	int i;
	  for (i = 0; i < N; ++i)
	    {
	    INTEGER(disc)[i] = discover_order[i];
	    }
	
	UNPROTECT(1);
	return(disc);
	}


	SEXP BGL_dfs_D(SEXP num_verts_in, SEXP num_edges_in, SEXP R_edges_in,
	    SEXP R_weights_in)
	{
	using namespace boost;
	
    typedef graph_traits < Graph_dd >::edge_descriptor Edge;
    typedef graph_traits < Graph_dd >::vertex_descriptor Vertex;
    Graph_dd g(num_verts_in, num_edges_in, R_edges_in, R_weights_in);
	
	typedef graph_traits < Graph_dd >::vertices_size_type size_type;
	
	const int N = INTEGER(num_verts_in)[0];
	  // Typedefs
	typedef size_type* Iiter;
	
	  // discover time and finish time properties
	std::vector < size_type > dtime(num_vertices(g));
	std::vector < size_type > ftime(num_vertices(g));
	size_type t = 0;
	dfs_time_visitor < size_type * >vis(&dtime[0], &ftime[0], t);
	
	depth_first_search(g, visitor(vis));
	
	  // use std::sort to order the vertices by their discover time
	std::vector < size_type > discover_order(N);
	integer_range < size_type > r(0, N);
	std::copy(r.begin(), r.end(), discover_order.begin());
	std::sort(discover_order.begin(), discover_order.end(),
	            indirect_cmp < Iiter, std::less < size_type > >(&dtime[0]));
	std::vector < size_type > finish_order(N);
	std::copy(r.begin(), r.end(), finish_order.begin());
	std::sort(finish_order.begin(), finish_order.end(),
	            indirect_cmp < Iiter, std::less < size_type > >(&ftime[0]));
	
	SEXP ansList;
	PROTECT(ansList = allocVector(VECSXP,2));
	SEXP disc;
	PROTECT(disc = allocVector(INTSXP,N));
	SEXP fin;
	PROTECT(fin = allocVector(INTSXP,N));
	
	int i;
	for (i = 0; i < N; ++i)
	    {
	    INTEGER(disc)[i] = discover_order[i];
	    INTEGER(fin)[i] = finish_order[i];
	    }
	
	SET_VECTOR_ELT(ansList,0,disc);
	SET_VECTOR_ELT(ansList,1,fin);
	UNPROTECT(3);
	return(ansList);
	}

	SEXP BGL_dijkstra_shortest_paths_D (SEXP num_verts_in, 
		SEXP num_edges_in, SEXP R_edges_in,
		SEXP R_weights_in, SEXP init_ind)
	{
	using namespace boost;
	
    typedef graph_traits < Graph_dd >::edge_descriptor Edge;
    typedef graph_traits < Graph_dd >::vertex_descriptor Vertex;
    Graph_dd g(num_verts_in, num_edges_in, R_edges_in, R_weights_in);
	
	int N = num_vertices(g);
	std::vector<Vertex> p(N);
	std::vector<double> d(N);
	
	dijkstra_shortest_paths(g, vertex((int)INTEGER(init_ind)[0], g),
	         predecessor_map(&p[0]).distance_map(&d[0]));
	
	SEXP dists, pens, ansList;
	PROTECT(dists = allocVector(REALSXP,N));
	PROTECT(pens = allocVector(INTSXP,N));
	graph_traits < Graph_dd >::vertex_iterator vi, vend;
	for (tie(vi, vend) = vertices(g); vi != vend; ++vi) {
	    REAL(dists)[*vi] = d[*vi];
	    INTEGER(pens)[*vi] = p[*vi];
	  }
	PROTECT(ansList = allocVector(VECSXP,2));
	SET_VECTOR_ELT(ansList,0,dists);
	SET_VECTOR_ELT(ansList,1,pens);
	
	UNPROTECT(3);
	return(ansList);
	}

	SEXP BGL_connected_components_U (SEXP num_verts_in, 
		SEXP num_edges_in, SEXP R_edges_in,
		SEXP R_weights_in )
	{
	using namespace boost;
        SEXP outvec;
	
    typedef graph_traits < Graph_ud >::edge_descriptor Edge;
    typedef graph_traits < Graph_ud >::vertex_descriptor Vertex;
    Graph_ud g(num_verts_in, num_edges_in, R_edges_in, R_weights_in);

	int nvert = INTEGER(num_verts_in)[0] ;

    std::vector<int> component(num_vertices(g));
    connected_components(g, &component[0]);

    std::vector<int>::size_type k;
    
	PROTECT(outvec = allocVector(REALSXP,nvert));

	for (k = 0; k < component.size(); k++ )
	    REAL(outvec)[k] = component[k];
	
	UNPROTECT(1);
	return(outvec);
	}

	SEXP BGL_strong_components_D (SEXP num_verts_in, 
		SEXP num_edges_in, SEXP R_edges_in,
		SEXP R_weights_in )
	{
	using namespace boost;
        SEXP outvec;
	
    typedef graph_traits < Graph_dd >::edge_descriptor Edge;
    typedef graph_traits < Graph_dd >::vertex_descriptor Vertex;
    Graph_dd g(num_verts_in, num_edges_in, R_edges_in, R_weights_in);

	int nvert = INTEGER(num_verts_in)[0] ;

    std::vector<int> component(num_vertices(g));
    strong_components(g, &component[0]);

    std::vector<int>::size_type k;
    
	PROTECT(outvec = allocVector(REALSXP,nvert));

	for (k = 0; k < component.size(); k++ )
	    REAL(outvec)[k] = component[k];
	
	UNPROTECT(1);
	return(outvec);
	}

	SEXP BGL_edge_connectivity_U (SEXP num_verts_in, 
		SEXP num_edges_in, SEXP R_edges_in,
		SEXP R_weights_in )
	{
	using namespace boost;
    SEXP ansList, conn, edTmp;
	
    typedef graph_traits < Graph_ud >::edge_descriptor Edge;
    typedef graph_traits < Graph_ud >::vertex_descriptor Vertex;
    Graph_ud g(num_verts_in, num_edges_in, R_edges_in, R_weights_in);

    typedef graph_traits<Graph_ud>::degree_size_type dst;
    std::vector<Edge> disconnecting_set;
    std::vector<Edge>::iterator ei;
    dst c = edge_connectivity( g, std::back_inserter(disconnecting_set) );

    PROTECT(conn = NEW_NUMERIC(1));
    REAL(conn)[0] = (double)c;

    SEXP eList;
    PROTECT(ansList = allocVector(VECSXP,2));

    PROTECT(eList = allocVector(VECSXP,(int)c));

    SET_VECTOR_ELT(ansList,0,conn);

    int sind = 0;
        for (ei = disconnecting_set.begin(); ei != disconnecting_set.end();
			++ei)
          {
          PROTECT(edTmp = NEW_NUMERIC(2));
          REAL(edTmp)[0] = (double)source(*ei,g);
          REAL(edTmp)[1] = (double)target(*ei,g);
          SET_VECTOR_ELT(eList,sind,edTmp);
          sind=sind+1;
          UNPROTECT(1);
          }

        SET_VECTOR_ELT(ansList,1,eList);
	UNPROTECT(3);
	return(ansList);
	}



	SEXP BGL_johnson_all_pairs_shortest_paths_D(SEXP num_verts_in, 
		SEXP num_edges_in, SEXP R_edges_in,
		SEXP R_weights_in)
	{
  using namespace boost;
  typedef adjacency_list<vecS, vecS, directedS, no_property,
    property< edge_weight_t, double, property< edge_weight2_t, double > > > Graph;
  int nv = INTEGER(num_verts_in)[0]; 
  if (nv > 200) error("bug in BGL limits num nodes to fixed number, now set at 200; you can recompile with larger limit if you like\n");
  SEXP out;
  const int V = nv;
  typedef std::pair < int, int >Edge;
  Edge edge_array[] =
    { Edge(0, 1), Edge(0, 4), Edge(0, 2), Edge(1, 3), Edge(1, 4),
    Edge(2, 1), Edge(3, 2), Edge(3, 0), Edge(4, 3)
  };
  const std::size_t E = sizeof(edge_array) / sizeof(Edge);

  Graph g(edge_array, edge_array + E, V);

  property_map < Graph, edge_weight_t >::type w = get(edge_weight, g);
  double weights[] = { 3, -4, 8, 1, 7, 4, -5, 2, 6 };
  double *wp = weights;

  graph_traits < Graph >::edge_iterator e, e_end;
  for (boost::tie(e, e_end) = edges(g); e != e_end; ++e)
    w[*e] = *wp++;

  double D[200][200];
  johnson_all_pairs_shortest_paths(g, D); 
  PROTECT(out = NEW_NUMERIC(nv*nv));
  int k = 0;
  for (int i = 0 ; i < nv ; i++)
   for (int j = 0; j < nv; j++ )
      {
      REAL(out)[k] = D[i][j];
      k++;
      }
  UNPROTECT(1);
  return out;
}

/*
SEXP BGL_transitive_closure_D (SEXP num_verts_in, 
		SEXP num_edges_in, SEXP R_edges_in )
	{
	using namespace boost;
	
    typedef graph_traits < Graph_dd >::edge_descriptor Edge;
    typedef graph_traits < Graph_dd >::vertex_descriptor Vertex;
    Graph_dd g(num_verts_in, num_edges_in, R_edges_in );
	
	int N = num_vertices(g);
   
        Graph_dd TC(num_verts_in, num_edges_in, R_edges_in);
	
	transitive_closure(g, TC);

        Graph_dd::edge_iterator e, eend;

	for (boost::tie(e,eend) = edges(TC);
		e != eend; ++e) {
			std::cout << "1" << "\n";
		}
	return R_NilValue;
	}
*/

}
