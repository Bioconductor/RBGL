#include "RBGL.hpp"

extern "C"
{
    using namespace boost;
    using namespace std;

    // Combinatorial Optimization: algorithms and complexity  (p. 403)
    // by C. H. Papadimitriou, K. Steiglitz
    //
    // A graph G = (V, E) is triangulated (i.e. chordal) if all cycles
    // [v1, v2, ..., vk] of length 4 or more have an chord, i.e., an edge
    // [vi, vj] with j != i +/- 1 (mod k)
    //
    // an equivalent definition of chordal graphs is:
    //
    // G is chordal iff either G is an empty graph, or
    // there is an v in V such that
    //     (i) the neighborhood of v (i.e., v and its adjacent nodes)
    //         forms a clique, and
    //     (ii) recursively, G-v is chordal

    static bool isClique(Graph_ud& g, graph_traits < Graph_ud >::vertex_descriptor v)
    {
        typedef graph_traits < Graph_ud >::edge_descriptor Edge;
        std::pair<Edge, bool> pp;

        Graph_ud::adjacency_iterator va, va_last, va1, va2;

        tie(va, va_last) = adjacent_vertices(v, g);


#if DEBUG
        cout << " no. nodes = " << num_vertices(g)
        << " no. edges = " << num_edges(g)
        << endl;
        cout << " connected nodes:   ";
        for ( va1 = va; va1 != va_last; va1++ )
        {
            cout << *va1+1 << "  ";
        }
        cout << endl;
#endif

        for ( va1 = va; va1 != va_last; va1++ )
        {
            va2 = va1;
            for ( va2 = va1+1; va2 != va_last; va2++ )
            {
#if DEBUG
                cout << "    " << *va1+1 << "--" << *va2+1 << endl;
#endif
                pp = edge(*va1, *va2, g);
                if ( !pp.second ) return FALSE;
            }
        }

        return TRUE;
    }

    static bool isTriangulatedInternal(Graph_ud& g)
    {
        bool r = TRUE;

        typedef graph_traits < Graph_ud >::vertex_iterator vertex_iterator;

        vertex_iterator vi, vi_end;
        Graph_ud::adjacency_iterator va, va_last;

        if ( num_edges(g) == 0 && num_vertices(g) == 0 )
            return TRUE;

        for ( tie(vi, vi_end) = vertices(g); vi != vi_end; vi++ )
        {
#if DEBUG
            cout << *vi+1 << endl;
#endif
            if ( isClique(g, *vi) )
            {
#if DEBUG
                cout << "         REMOVED " << endl;
#endif
                clear_vertex(*vi, g);
                remove_vertex(*vi, g);
                r = isTriangulatedInternal(g);
                return r;
            }
        }
        return FALSE;
    }

    SEXP isTriangulated(SEXP num_verts_in, SEXP num_edges_in, SEXP R_edges_in )
    {
        Graph_ud g(num_verts_in, num_edges_in, R_edges_in);
        SEXP ans;
        PROTECT(ans = NEW_INTEGER(1));

        INTEGER(ans)[0] = isTriangulatedInternal(g);

        UNPROTECT(1);
        return(ans);
    }

typedef vector<int>              oneCliqueType;
typedef vector<oneCliqueType>    allCliquesType;

#if DEBUG
    static void print_one_clique(oneCliqueType& clique)
    {
        cout << "     clique contains node(s): ";
        for ( unsigned int j = 0; j < clique.size(); j++ )
            cout << clique[j]+1 << "  "; 
        cout << endl;
    }

    static void print_all_cliques(allCliquesType& cliques, char* msg)
    {
        cout << msg << endl;
        for ( unsigned int i = 0; i < cliques.size(); i++ )
            print_one_clique(cliques[i]);
    }
#endif

    static bool inline isConnected(Graph_ud& g, int u, int v)
    {
        if ( u == v ) return TRUE;

        typedef graph_traits < Graph_ud >::edge_descriptor Edge;
        std::pair<Edge, bool> p;
        p = edge(u, v, g);
        return p.second;
    }

    static void extend_v2(Graph_ud& g, std::vector<int>& old_array, 
                          std::vector<int>& compsub,
                          int ne, int ce, int& c, 
                          allCliquesType& cliques)
    {
        std::vector<int> new_array(ce+1);    
        int nod, fixp;
        int newne, newce, i, j, count, pos, p, s, sel, minnod;
        minnod = ce;
        nod = 0;
        fixp = -1;
        pos = -1;
        s = -1;

        // determine each counter value and look for minimum
        for ( i = 1; i <= ce && minnod != 0; i++ )
        {
            p = old_array[i]; count = 0;

            // count disconnections
            for ( j = ne+1; j <= ce && count < minnod; j++ )
            {
                if ( !isConnected(g, p, old_array[j]) )
                {
                    count++;
                    pos = j;  // position of potential candidate
                }
            }

            // test new minimum
            if ( count < minnod )
            {
                fixp = p;
                minnod = count;
                if ( i <= ne ) { s = pos; }
                else { s = i; nod = 1; }
            }
        }

        // backtrack cycle
        if ( s <= ce )         // just to play it safe
        for ( nod = minnod + nod; nod >= 1; nod-- )
        {
            // interchange
            p = old_array[s];
            old_array[s] = old_array[ne+1];
            old_array[ne+1] = p;
            sel = p;

            newne = 0;

            // fill new set "not_set"
            for ( i = 1; i <= ne; i++ )
            {
                if ( isConnected(g, sel, old_array[i]) )
                    new_array[++newne] = old_array[i];
            }

            // fill new set "cand"
            newce = newne; 
            for ( i = ne+2; i <= ce; i++ )
            {
                if ( isConnected(g, sel, old_array[i]) )
                    new_array[++newce] = old_array[i];
            }

            c++;
            compsub[c] = sel;
            if ( newce == 0 )
            {
                oneCliqueType clique(c);
                for ( int loc = 1; loc <= c; loc++ )
                {
                    clique[loc-1] = compsub[loc];
                }
                cliques.push_back(clique);

#if DEBUG
                cout << "   found ";
                print_one_clique(clique);
#endif
            }
            else
            {
                if ( newne < newce )
                    extend_v2(g, new_array, compsub, newne, newce, c, cliques);
            }

            // remove from "compsub"
            c--;

            // add to "not_set"
            ne++;
            // select a condidate disconnected to the fixed point
            if ( nod > 1 )
            {
                // look for candidate
                for ( s = ne+1; s <= ce && isConnected(g, fixp, old_array[s]); s++ )
                    ;
            }

        }
    }

    static void bron_kerbosch_all_cliques(Graph_ud& g, allCliquesType& cliques)
    {
        int N = num_vertices(g);

        std::vector<int> ALL(N+1);
        std::vector<int> compsub(N+1, 0);

        for ( int i = 0; i <= N; i++ ) ALL[i] = i-1;  // node index for BGL graph

        int c = 0;
        extend_v2(g, ALL, compsub, 0, N, c, cliques);

#if DEBUG
        print_all_cliques(cliques, " all cliques: ");
#endif
    }

    SEXP maxClique(SEXP num_verts_in, SEXP num_edges_in, SEXP R_edges_in )
    {
        Graph_ud g(num_verts_in, num_edges_in, R_edges_in);

        allCliquesType cliques;

        bron_kerbosch_all_cliques(g, cliques);

        allCliquesType::iterator ci;

/*
        // keep only the max cliques
        int max_clique_size = 0;
        for ( ci = cliques.begin(); ci != cliques.end(); ci++ )
        {
            if ( (*ci).size() > max_clique_size )
               max_clique_size = (*ci).size();
        }

        for ( ci--; ci != cliques.begin(); ci-- )
            if ( (*ci).size() < max_clique_size )
               cliques.erase(ci);
*/

#if DEBUG
        print_all_cliques(cliques, " max clique(s): ");
#endif

        oneCliqueType::iterator cj;

        int i, j;
        SEXP ansList, cnodes;
        PROTECT(ansList = allocVector(VECSXP, cliques.size()));

        for ( i = 0, ci = cliques.begin(); ci != cliques.end(); i++, ci++ )
        {
            PROTECT(cnodes = allocVector(INTSXP, (*ci).size()));

            for ( j = 0, cj = (*ci).begin(); cj != (*ci).end(); j++, cj++ )
                INTEGER(cnodes)[j] = *cj+1;    // node index for R graph

            SET_VECTOR_ELT(ansList,i,cnodes);
            UNPROTECT(1);
        }

        UNPROTECT(1);
        return(ansList);
    }
}

