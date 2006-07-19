#include "RBGL.hpp"

extern "C"
{

#include <Rdefines.h>

    using namespace std;
    using namespace boost;

    static void delta_and_tau
    (SEXP num_verts_in, SEXP num_edges_in, SEXP R_edges_in,
     vector<int>& v_delta, vector<int>& v_tau)
    {
        Graph_ud g(num_verts_in, num_edges_in, R_edges_in);

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

        vector<int> v_delta, v_tau;
        delta_and_tau(num_verts_in, num_edges_in, R_edges_in, v_delta, v_tau);

        double nn = 0;		// count nodes w/ deg(v) >= 2
        double cG = 0;
        for ( i = 0; i < NV; i++ )
        {
            if ( v_tau[i] > 0 )
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

        vector<int> v_delta, v_tau;
        delta_and_tau(num_verts_in, num_edges_in, R_edges_in, v_delta, v_tau);

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

}

