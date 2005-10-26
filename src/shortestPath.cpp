#include "RBGL.hpp"
#include "Basic2DMatrix.hpp"
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/johnson_all_pairs_shortest.hpp>
#include <boost/graph/dag_shortest_paths.hpp>
#include <boost/graph/bellman_ford_shortest_paths.hpp>
#include <boost/graph/floyd_warshall_shortest.hpp>

extern "C"
{
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

    SEXP BGL_johnson_all_pairs_shortest_paths_D(SEXP num_verts_in,
            SEXP num_edges_in, SEXP R_edges_in,
            SEXP R_weights_in)
    {
        using namespace boost;
        typedef adjacency_list<vecS, vecS, directedS, no_property,
        property< edge_weight_t, double, property< edge_weight2_t, double > > > Graph;
        int nv = INTEGER(num_verts_in)[0];
        SEXP out;
        const int V = nv;
        typedef std::pair < int, int >Edge;

        Graph_dd g(num_verts_in, num_edges_in, R_edges_in, R_weights_in);

        Basic2DMatrix<double> D(nv, nv);

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

    SEXP BGL_bellman_ford_shortest_paths(SEXP num_verts_in,
                                         SEXP num_edges_in, SEXP R_edges_in,
                                         SEXP R_weights_in, SEXP init_ind)
    {
        using namespace boost;

        typedef adjacency_list < vecS, vecS, directedS,
        no_property, property < edge_weight_t, double> > EdgeGraph;

        int i, j;
        int NE = (int)INTEGER(num_edges_in)[0];
        int N = (int)INTEGER(num_verts_in)[0];
        int s = (int)INTEGER(init_ind)[0];

        EdgeGraph g(N);

        int* edges_in = INTEGER(R_edges_in);
        for ( i = 0; i < NE; i++, edges_in += 2 )
            add_edge(*edges_in, *(edges_in+1), g);

        std::vector<std::size_t> p(N);
        for ( i = 0; i < N; i++ ) p[i] = i;

        std::vector<double> d(N, std::numeric_limits<double>::max());
        d[s] = 0;

        property_map<EdgeGraph, edge_weight_t>::type w = get(edge_weight, g);

        int* weight_i = (isReal(R_weights_in)) ? 0 : INTEGER(R_weights_in);
        double* weight_d = (isReal(R_weights_in)) ? REAL(R_weights_in) : 0;

        graph_traits< EdgeGraph >::edge_iterator ei, ei_end;
        for ( tie(ei, ei_end) = edges(g); ei != ei_end; ++ei )
            w[*ei] = weight_i ? (*weight_i++) : (*weight_d++);

        bool r = bellman_ford_shortest_paths(g, N,
                 weight_map(w).predecessor_map(&p[0]).distance_map(&d[0]));

        SEXP conn, dList, pList, ansList;
        PROTECT(ansList = allocVector(VECSXP,3));
        PROTECT(conn = NEW_LOGICAL(1));
        PROTECT(dList = allocVector(REALSXP,N));
        PROTECT(pList = allocVector(INTSXP,N));

        LOGICAL(conn)[0] = r;

        for (i = 0; i < N; i++)
        {
            INTEGER(pList)[i] = p[i];
            REAL(dList)[i] = d[i];
        }

        SET_VECTOR_ELT(ansList,0,conn);
        SET_VECTOR_ELT(ansList,1,dList);
        SET_VECTOR_ELT(ansList,2,pList);

        UNPROTECT(4);
        return(ansList);
    }

    SEXP BGL_dag_shortest_paths(SEXP num_verts_in,
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

        dag_shortest_paths(g, vertex((int)INTEGER(init_ind)[0], g),
                           predecessor_map(&p[0]).distance_map(&d[0]));

        SEXP dists, pens, ansList;
        PROTECT(dists = allocVector(REALSXP,N));
        PROTECT(pens = allocVector(INTSXP,N));
        graph_traits < Graph_dd >::vertex_iterator vi, vend;
        for (tie(vi, vend) = vertices(g); vi != vend; ++vi) {
            if ( int(d[*vi]) == std::numeric_limits<int>::max() )
            {
                REAL(dists)[*vi] = R_NaN;
                INTEGER(pens)[*vi] = *vi;
            }
            else
            {
                REAL(dists)[*vi] = d[*vi];
                INTEGER(pens)[*vi] = p[*vi];
            }
        }
        PROTECT(ansList = allocVector(VECSXP,2));
        SET_VECTOR_ELT(ansList,0,dists);
        SET_VECTOR_ELT(ansList,1,pens);

        UNPROTECT(3);
        return(ansList);
    }

    SEXP BGL_floyd_warshall_all_pairs_shortest_paths_D(SEXP num_verts_in,
            SEXP num_edges_in, SEXP R_edges_in,
            SEXP R_weights_in)
    {
        using namespace boost;
        typedef adjacency_list<vecS, vecS, directedS, no_property,
        property< edge_weight_t, double, property< edge_weight2_t, double > > > Graph;
        int nv = INTEGER(num_verts_in)[0];
        SEXP out;
        const int V = nv;
        typedef std::pair < int, int >Edge;

        Graph_dd g(num_verts_in, num_edges_in, R_edges_in, R_weights_in);

        Basic2DMatrix<double> D(nv, nv);

        floyd_warshall_all_pairs_shortest_paths(g, D);

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
}

