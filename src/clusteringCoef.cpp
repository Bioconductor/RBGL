#include "RBGL.hpp"

#include <stdlib.h>

#include <boost/graph/simple_point.hpp>

extern "C"
{

#include <Rdefines.h>
#include <R_ext/Random.h>
#include <Rmath.h>

    using namespace std;
    using namespace boost;

    static void delta_and_tau
     (const Graph_ud& g, vector<int>& v_delta, vector<int>& v_tau)
    {
        Graph_ud::vertex_iterator vi, v_end;
        Graph_ud::adjacency_iterator ui, u_end, wi, w_end;

        int dv = 0, tv = 0;

        v_delta.clear();
        v_tau.clear();

        for ( tie(vi, v_end) = vertices(g); vi != v_end; ++vi )
        {
            // delta(v)
            dv = 0;
            for ( tie(ui, u_end) = adjacent_vertices(*vi, g);
                    ui != u_end; ++ui )
            {
                wi = ui;
                for ( ++wi; wi != u_end; ++wi )
                    if ( edge(*ui, *wi, g).second ) dv++;
            }
            v_delta.push_back(dv);

            // tau(v)
            dv = degree(*vi, g);
            tv = dv * ( dv - 1 ) / 2;
            v_tau.push_back(tv);
        }
    }

    SEXP clusteringCoef(
        SEXP num_verts_in, SEXP num_edges_in, SEXP R_edges_in,
        SEXP weighted, SEXP R_v_weights_in)
    {
        int i;

        int NV = INTEGER(num_verts_in)[0];
        vector<double> v_weight(NV, 1);

        if ( INTEGER(weighted)[0] )
        {
            double* weights = REAL(R_v_weights_in);
            for ( i = 0; i < NV; i++ ) v_weight[i] = weights[i];
        }

        Graph_ud g(num_verts_in, num_edges_in, R_edges_in);
        vector<int> v_delta, v_tau;
        delta_and_tau(g, v_delta, v_tau);

        double nn = 0;		// count nodes w/ deg(v) >= 2
        double cG = 0;
	Graph_ud::vertex_descriptor v;
        for ( i = 0; i < NV; i++ )
        {
	    v = vertex(i, g);
            if ( out_degree(v, g) >= 2 && v_tau[i] > 0 )
            {
                cG += v_weight[i] * v_delta[i] / v_tau[i];
                nn += v_weight[i];
            }
        }

        if ( nn ) cG /= nn;

        SEXP ccoef;
        PROTECT(ccoef = NEW_NUMERIC(1));
        REAL(ccoef)[0] = cG;
        UNPROTECT(1);
        return(ccoef);
    }

    SEXP transitivity( SEXP num_verts_in, SEXP num_edges_in, SEXP R_edges_in)
    {
        int NV = INTEGER(num_verts_in)[0];

        Graph_ud g(num_verts_in, num_edges_in, R_edges_in);
        vector<int> v_delta, v_tau;
        delta_and_tau(g, v_delta, v_tau);

	double tG = 0;
        double sum_dv = 0, sum_tv = 0;
        for ( int i = 0; i < NV; i++ )
        {
            sum_dv += v_delta[i];
            sum_tv += v_tau[i];
        }
        if ( sum_tv ) tG = sum_dv / sum_tv;

#if DEBUG
        cout << " sum_dv = " << sum_dv
        << " sum_tv = " << sum_tv
        << " v_delta.size() = " << v_delta.size()
        << " v_tau.size() = " << v_tau.size()
	<< " tG = " << tG
        << endl;
#endif

        SEXP tcoef;
        PROTECT(tcoef = NEW_NUMERIC(1));
        REAL(tcoef)[0] = tG;
        UNPROTECT(1);
        return(tcoef);
    }

    // uniformly sample in [1, n]
    static inline int uniformRandomNumber(const int n)
    {
	int j = (int)(n * unif_rand()) + 1;   // unif_rand in [0, 1)
	return j;
    }

    static inline int findIndex(const int r, const vector<int>& W)
    {
	unsigned int i;
	for ( i = 1; i < W.size(); i++ ) if ( r <= W[i] ) break;
	return i;
    }

    // input: a node in graph g
    // output: one neighbor of node n chosen uniformly randomly
    static inline void uniformRandomAdjacentNode
	(const Graph_ud::vertex_descriptor& v, const Graph_ud& g, 
	 Graph_ud::vertex_descriptor& u,
	 Graph_ud::vertex_descriptor& w)
    {
	int nc = out_degree(v, g);

	Graph_ud::adjacency_iterator vi, v_end;
	tie(vi, v_end) = adjacent_vertices(v, g);

	switch (nc)
	{
	case 0: 
	case 1: u = w = *vi;
		break;
	case 2: 
		u = *vi; vi++;
		w = *vi; 
		break;
	default:
		{
		int r1 = uniformRandomNumber(nc);
		int r2 = uniformRandomNumber(nc);

		while ( r1 == r2 ) r2 = uniformRandomNumber(nc);

		for ( int i = 0; vi != v_end; vi++, i++ )
		{
		    if ( i == r1 ) u = *vi;
		    if ( i == r2 ) w = *vi;
		}

		break;
		}
	}
#if DEBUG
	cout << " uniformRandomAdjacentNode: " << endl;
	cout << " n = " << n << endl;
	cout << " nc = " << nc << endl;
	cout << " *vi = " << *vi << endl;
	cout << " u = " << u << endl;
	cout << " w = " << w << endl;
#endif
    }

    static inline void uniformRandomAdjacentNode_i
	(const int n, const Graph_ud& g, 
	 Graph_ud::vertex_descriptor& u,
	 Graph_ud::vertex_descriptor& w)
    {
	Graph_ud::vertex_descriptor v = vertex(n, g);
	uniformRandomAdjacentNode(v, g, u, w);
    }

    // Approximating Cw  
    //    Outline of the algorithm:
    //    Input: integer k; 
    //           array A[1..|V'|] of nodes V' = {v in V: d(v) >= 2}
    //           node weights w: V' -> N>0;
    //           adjacentcy array for each node
    //    Output: approximation of Cw
    //    Data: node variables: u, w;
    //           integer variables: r, l, j, W[0..|V'|]
    //    Algorithm:
    //    W[0] = 0
    //    for i = (1, ..., |V'|) do
    //       W[i] = W[i-1] + w(A[i])
    //    l = 0
    //    for i in (1, ..., k) do
    //    {
    //       r = UniformRandomNumber( {1,...,W[|V'|]} )
    //       j = FindIndex( j: W[j-1] < r <= W[j] )
    //       u = UniformRandomAdjacentNode(A[j])
    //       repeat
    //         w = UniformRandomAdjacentNode(A[j])
    //       until u != w
    //       if ( EdgeExists(u, w) then
    //          l = l + 1
    //    }
    //    return l/k

    SEXP clusteringCoefAppr(SEXP k_in,
        SEXP num_verts_in, SEXP num_edges_in, SEXP R_edges_in,
        SEXP weighted, SEXP R_v_weights_in)
    {
	// prepare for later unif_rand call
	GetRNGstate();

        Graph_ud g(num_verts_in, num_edges_in, R_edges_in);

        int i, l, r, j;

        int k = INTEGER(k_in)[0];
        int NV = INTEGER(num_verts_in)[0];
        vector<int> v_weight(NV, 1);
        vector<int> W(NV+1, 0);

#if DEBUG
	cout << " inside clusteringCoefAppr " 
		<< " k = " << k
		<< " NV = " << NV 
		<< endl;
#endif

        if ( INTEGER(weighted)[0] )
        {
            double* weights = REAL(R_v_weights_in);
            for ( i = 0; i < NV; i++ ) v_weight[i] = (int) weights[i];
        }

        Graph_ud::vertex_descriptor u=Graph_ud::null_vertex(), w=Graph_ud::null_vertex();

	W[0] = 0;
	for ( i = 1; i < NV+1; i++ ) W[i] = W[i-1] + v_weight[i-1];

	// TODO: limit nodes to those w/ degree >= 2
	//       pick a number within range uniformaly 
	//	 pick an adjacent node uniformaly randomly
	for ( l = 0, i = 0; i < k; i++ )
	{
	   r = uniformRandomNumber(W[NV]);
	   j = findIndex(r, W);
	   uniformRandomAdjacentNode_i(j-1, g, u, w);

	   if ( edge(u, w, g).second ) l++;

#if DEBUG
	   cout << " i = " << i;
	   cout << " r = " << r;
	   cout << " j = " << j;
	   cout << " l = " << l << endl;
#endif
	}

	double cG = double(l) / double(k);

        SEXP ccoef;
        PROTECT(ccoef = NEW_NUMERIC(1));
        REAL(ccoef)[0] = cG;
        UNPROTECT(1);
        return(ccoef);

    }

    inline bool prob_cmp(const simple_point<int>& p1, 
			 const simple_point<int>& p2)
	{ return p1.y > p2.y; }

    //  To find a random node w/ probability d(u) / sum(d(V))
    //  The following closely mirrors the codes on
    //  Unequal probability sampling; without-replacement case
    //
    //    /* Record element identities */
    //    for (i = 0; i < n; i++)
    //        perm[i] = i + 1;
    //
    //    //* Sort probabilities into descending order */
    //    //* Order element identities in parallel */
    //    revsort(p, perm, n);
    //
    //    //* Compute the sample */
    //    totalmass = 1;
    //    for (i = 0, n1 = n-1; i < nans; i++, n1--) {
    //        rT = totalmass * unif_rand();
    //        mass = 0;
    //        for (j = 0; j < n1; j++) {
    //            mass += p[j];
    //            if (rT <= mass)
    //                break;
    //        }
    //        ans[i] = perm[j];
    //        totalmass -= p[j];
    //        for(k = j; k < n1; k++) {
    //            p[k] = p[k + 1];
    //            perm[k] = perm[k + 1];
    //        }
    //    }
    static void ProbRandomNode
	(const Graph_ud::vertex_descriptor& v, const Graph_ud& g,
	 Graph_ud::vertex_descriptor& u)
    {

	typedef graph_traits<Graph_ud>::vertex_iterator vertex_iterator;

	vertex_iterator vi, v_end;

	int NV = num_vertices(g);
	std::vector < simple_point<int> > pp(num_vertices(g));

	int i = 0, totalmass = 0;
	for ( tie(vi, v_end) = vertices(g); vi != v_end; vi++, i++ )
	{
	   pp[i].x = i+1;
	   pp[i].y = out_degree(*vi, g);
	   totalmass += pp[i].y;
	}

	std::stable_sort(pp.begin(), pp.end(), prob_cmp);

	int j, k, n1, rT, mass;
	for ( i = 0, n1 = NV-1; i < NV; i++, n1-- )
	{
            rT = (int)(totalmass * unif_rand());
            mass = 0;
            for (j = 0; j < n1; j++) {
                mass += pp[j].y;
                if (rT <= mass) break;
            }
	    u = vertex(i, g);
	    if ( !edge(v, u, g).second ) break;

	    totalmass -= pp[j].y;
	    for ( k = j; k < n1; k++ ) pp[k] = pp[k+1];
	}
    }

    // Graph Generator:
    //    Outline of algorithm:
    //    Input: initial graph G: two connected nodes
    //           integer: n >= 3, d >= 2, o
    //    Output: graph G
    //    Algorithm:
    //    for ( i = (3, ..., n) do
    //    {
    //       v = NewNode()
    //       for 1, ..., Min(i-1, d) do
    //       {
    //          repeat 
    //            u = RandomNode( with prob du / sum(d(v) )
    //          until node EdgeExists(v, u)
    //          AddEdge(v, u)
    //       }
    //       for 1, ..., o do
    //       {
    //          u = RandomAdjacentNode(v)
    //          repeat
    //             w = RandomAdjacentNode(v)
    //          until w != u
    //          if ( node EdgeExists(u, w) then
    //             AddEdge(u, w)
    //       }
    //    }

    SEXP graphGenerator(SEXP n_in, SEXP d_in, SEXP o_in)
    {
	int i, j;
        int n = INTEGER(n_in)[0];
        int d = INTEGER(d_in)[0];
        int o = INTEGER(o_in)[0];

	GetRNGstate();	// get random number generator ready

	// initial graph with 2 connected nodes
	Graph_ud g(2);
	boost::add_edge(0, 1, g);

        Graph_ud::vertex_descriptor v, u, w=Graph_ud::null_vertex();
	
	for ( i = 3; i <= n; i++ )
	{
	   // generate a new node 
	   v = boost::add_vertex(g);

	   for ( j = 1; j <= min(i-1, d); j++ )
	   {
		ProbRandomNode(v, g, u); 
		boost::add_edge(v, u, g);
	   }

	   for ( j = 1; j <= o; j++ )
	   {
	        uniformRandomAdjacentNode(v, g, u, w);

		if ( !edge(u, w, g).second )
		   boost::add_edge(u, w, g);
	   }
	}

#if DEBUG
	typedef graph_traits<Graph_ud>::vertex_iterator vertex_iterator;
	vertex_iterator vi, v_end;
	cout << " no. of vertices: " << num_vertices(g)
	     << " no. of edges:    " << num_edges(g)
	     << endl;

	for ( tie(vi, v_end) = vertices(g); vi != v_end; vi++, i++ )
	{
	   cout << " vertex: " << *vi 
		<< " has degree: " << out_degree(*vi, g)
		<< endl;
	}
#endif

        int NE = num_edges(g);
        SEXP anslst, ncnt, ecnt, enlst;
        PROTECT(anslst = allocVector(VECSXP, 3));
        PROTECT(ncnt = NEW_INTEGER(1));
        PROTECT(ecnt = NEW_INTEGER(1));
        PROTECT(enlst = allocMatrix(INTSXP, 2, NE));

        INTEGER(ncnt)[0] = num_vertices(g);
        INTEGER(ecnt)[0] = NE;

	typedef graph_traits<Graph_ud>::edge_iterator edge_iterator;
	edge_iterator ei, e_end;
	for ( i = 0, tie(ei, e_end) = edges(g); ei != e_end ; ei++ )
        {
            INTEGER(enlst)[i++] = source(*ei, g);
            INTEGER(enlst)[i++] = target(*ei, g);
        }


	SET_VECTOR_ELT(anslst,0,ncnt);
	SET_VECTOR_ELT(anslst,1,ecnt);
        SET_VECTOR_ELT(anslst,2,enlst);
        UNPROTECT(4);
        return(anslst);
    }

}

