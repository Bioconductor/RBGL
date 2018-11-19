#include "RBGL.hpp"
#include "Basic2DMatrix.hpp"
#include <boost/graph/johnson_all_pairs_shortest.hpp>
#include <boost/graph/edmonds_karp_max_flow.hpp>

extern "C"
{

#include <Rdefines.h>

    using namespace std;
    using namespace boost;

    typedef std::set<int> SubgraphAsSet;
    typedef std::vector< SubgraphAsSet > CliqueVector;
    typedef std::vector< CliqueVector > ResultCliquesType;

    static void addNewClique(CliqueVector& cliques, int i, int j)
    {
        SubgraphAsSet s;
        s.insert(i);
        s.insert(j);
        cliques.push_back(s);
    }

    static void findAllCliques(ResultCliquesType& rCliques,
                               Basic2DMatrix<double>& D)
    {
        CliqueVector cliques;
        SubgraphAsSet::iterator s;
        CliqueVector::iterator ci, cj;
        int i, j, k, N=0;
        const int nv = D.numrows();

        // N: max distance in given graph
        for ( i = 0; i < nv; i++ )
            for ( j = i+1; j < nv; j++ )
            {
                N = max(N, (int)D[i][j]);
                // each edge is 1-clique
                if ( D[i][j] == 1 ) addNewClique(cliques, i, j);
            }

        for ( k = 1; k <= N; k++ )
        {
            for ( i = 0; i < nv; i++ )
            {
                for ( ci = cliques.begin(); ci != cliques.end(); ci++ )
                {
                    // i is already in this clique
                    if ( (*ci).find(i) != (*ci).end() ) continue;

                    for ( s = (*ci).begin(); s != (*ci).end(); s++ )
                    {
                        if ( D[i][*s] > k || D[*s][i] > k ) break;
                    }

                    // add i to this clique
                    if ( s == (*ci).end() )
                    {
                        (*ci).insert(i);

                        // eliminate its subsequent subsets
                        for ( cj = ci + 1; cj != cliques.end(); )
                        {
                            if (includes((*ci).begin(), (*ci).end(),
                                         (*cj).begin(), (*cj).end()) )
                                cj = cliques.erase(cj);
                            else
                                cj++;
                        }
                    }
                }
            }
            rCliques.push_back(cliques);
        }

#if DEBUG
        cout << " Cliques: " << endl;
        for ( i = 0; i < rCliques.size(); i++ )
        {
            cout << i+1 << " cliques: " << endl;
            for ( ci=rCliques[i].begin(); ci!=rCliques[i].end(); ci++ )
            {
                cout << "    ";
                for ( s = (*ci).begin(); s != (*ci).end(); s++ )
                    cout << (*s)+1 << " ";
                cout << endl;
            }
        }
#endif
    }

    SEXP kCliques(SEXP num_verts_in, SEXP num_edges_in,
                  SEXP R_edges_in, SEXP R_weights_in)
    {
        // R_weights_in has to be INTEGER now
        int nv = INTEGER(num_verts_in)[0];

        Graph_ud g(num_verts_in, num_edges_in, R_edges_in, R_weights_in);

        Basic2DMatrix<double> D(nv, nv);

        // find out the shortest distance between any two nodes
        johnson_all_pairs_shortest_paths(g, D);

        // find k-cliques now
        ResultCliquesType rCliques;
        findAllCliques(rCliques, D);

        ResultCliquesType::iterator ci;
        CliqueVector::iterator vi;
        SubgraphAsSet::iterator si;
        int i, j, k;

        SEXP ansList, cList, sList;
        PROTECT(ansList = Rf_allocVector(VECSXP, (int)rCliques.size()));

        for ( i = 0, ci = rCliques.begin(); ci != rCliques.end(); i++, ci++)
        {
            PROTECT(cList = Rf_allocVector(VECSXP, (*ci).size()));
            for ( j = 0, vi = (*ci).begin(); vi != (*ci).end(); j++, vi++ )
            {
                PROTECT(sList = Rf_allocVector(INTSXP, (*vi).size()));
                for ( k = 0, si = (*vi).begin(); si != (*vi).end(); k++, si++ )
                {
                    INTEGER(sList)[k] = *si;
                }
                SET_VECTOR_ELT(cList,j,sList);
                UNPROTECT(1);
            }
            SET_VECTOR_ELT(ansList,i,cList);
            UNPROTECT(1);
        }
        UNPROTECT(1);
        return(ansList);
    }

    SEXP lambdaSets(SEXP num_verts_in, SEXP num_edges_in,
                    SEXP R_edges_in, SEXP R_capacity_in)
    {
        using namespace boost;

        typedef adjacency_list_traits<vecS, vecS, directedS> Tr;
        typedef Tr::edge_descriptor Tr_edge_desc;

        typedef adjacency_list<vecS, vecS, directedS, no_property,
        property<edge_capacity_t, double,
        property<edge_residual_capacity_t, double,
        property<edge_reverse_t, Tr_edge_desc> > > >
        FlowGraph;

        typedef graph_traits<FlowGraph>::vertex_descriptor vertex_descriptor;
	typedef graph_traits<FlowGraph>::edge_descriptor edge_descriptor;

        FlowGraph flow_g;

        property_map < FlowGraph, edge_capacity_t >::type
        cap = get(edge_capacity, flow_g);
        property_map < FlowGraph, edge_residual_capacity_t >::type
        res_cap = get(edge_residual_capacity, flow_g);
        property_map < FlowGraph, edge_reverse_t >::type
        rev_edge = get(edge_reverse, flow_g);

        edge_descriptor e1, e2;
        bool in1, in2;

        if (!Rf_isInteger(R_edges_in)) Rf_error("R_edges_in should be integer");

        int NV = INTEGER(num_verts_in)[0];
        int NE = Rf_asInteger(num_edges_in);
        int* edges_in = INTEGER(R_edges_in);
	int i, j, k, MaxC=0; 

        for (i = 0; i < NE ; i++, edges_in += 2)
        {
            tie(e1, in1) = boost::add_edge(*edges_in, *(edges_in+1), flow_g);
            tie(e2, in2) = boost::add_edge(*(edges_in+1), *edges_in, flow_g);
            if ( !in1 || !in2 )
                Rf_error("unable to add edge: (%d, %d)", *edges_in, *(edges_in+1));

            // fill in capacity_map
            cap[e1] = 1; 
            cap[e2] = 1; 

            // fill in reverse_edge_map
            rev_edge[e1] = e2;
            rev_edge[e2] = e1;
        }

        // ASSUMPTION: max_flow(u, v) = max_flow(v, u)
        // compute edge connectivities between u and v, (u < v)
	// we only need lower-left triangle of the matrix w/o diagonal
        Basic2DMatrix<int> CC(NV, NV);

        for ( i = 0; i < NV; i++ )
        {
            vertex_descriptor s = vertex(i, flow_g);

            for ( j = 0; j < i; j++ )
            {
                vertex_descriptor t = vertex(j, flow_g);

                CC[i][j] = (int)edmonds_karp_max_flow(flow_g, s, t);
                MaxC = max(MaxC, CC[i][j]);
            }
        }

        // calc lambda sets by successively partition V
        Basic2DMatrix<int> P(MaxC+1, NV);
        for ( k = 0; k <= MaxC; k++ )
        {
            for ( i = 0; i < NV; i++ ) P[k][i] = i;

            for ( i = 1; i < NV; i++ )
                for ( j = 0; j < i; j++ )
                    if ( CC[i][j] >= k )  P[k][i] = P[k][j];
        }

#if DEBUG
        cout << " edge-connectivity matrix: " << endl;
        for ( i = 0; i < NV; i++ )
        {
            for ( j = 0; j < NV; j++ ) cout << CC[i][j] << " ";
            cout << endl;
        }

        cout << " P matrix: " << endl;
        for ( k = 0; k <= MaxC; k++ )
        {
            cout << " k = " << k << ": ";
            for ( j = 0; j < NV; j++ ) cout << P[k][j] << " ";
            cout << endl;
        }
#endif

        SEXP ansList, conn, eList;
        PROTECT(ansList = Rf_allocVector(VECSXP,2));
        PROTECT(conn = Rf_allocVector(REALSXP, 1));
        PROTECT(eList = Rf_allocMatrix(INTSXP, MaxC+1, NV));

        REAL(conn)[0] = MaxC;

        for ( i = 0, j = 0; j < NV; j++ )
            for ( k = 0; k <= MaxC; k++ )
                INTEGER(eList)[i++] = P[k][j];

        SET_VECTOR_ELT(ansList,0,conn);
        SET_VECTOR_ELT(ansList,1,eList);
        UNPROTECT(3);

        return(ansList);
    }
}

