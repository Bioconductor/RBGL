#include "RBGL.hpp"
#include "Basic2DMatrix.hpp"
#include <boost/graph/johnson_all_pairs_shortest.hpp>


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
        PROTECT(ansList = allocVector(VECSXP, (int)rCliques.size()));

        for ( i = 0, ci = rCliques.begin(); ci != rCliques.end(); i++, ci++)
        {
            PROTECT(cList = allocVector(VECSXP, (*ci).size()));
            for ( j = 0, vi = (*ci).begin(); vi != (*ci).end(); j++, vi++ )
            {
                PROTECT(sList = allocVector(INTSXP, (*vi).size()));
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

}

