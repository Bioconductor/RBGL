#include "RBGL.hpp"
#include <boost/graph/bandwidth.hpp>
#include <boost/graph/profile.hpp>
#include <boost/graph/wavefront.hpp>

extern "C"
{

	SEXP BGL_bandwidth(SEXP num_verts_in, SEXP num_edges_in, SEXP R_edges_in )
	{
		using namespace boost;

		Graph_ud g(num_verts_in, num_edges_in, R_edges_in);

		SEXP ansList, rbw;
		PROTECT(ansList = Rf_allocVector(VECSXP,1));
		PROTECT(rbw = Rf_allocVector(INTSXP, 1));

		INTEGER(rbw)[0] = bandwidth(g);

		SET_VECTOR_ELT(ansList,0,rbw);
		UNPROTECT(2);
		return(ansList);
	}

	SEXP BGL_profile(SEXP num_verts_in, SEXP num_edges_in, SEXP R_edges_in )
	{
		using namespace boost;

		Graph_ud g(num_verts_in, num_edges_in, R_edges_in);

		SEXP ansList, rbw;
		PROTECT(ansList = Rf_allocVector(VECSXP,1));
		PROTECT(rbw = Rf_allocVector(INTSXP, 1));

		INTEGER(rbw)[0] = profile(g);

		SET_VECTOR_ELT(ansList,0,rbw);
		UNPROTECT(2);
		return(ansList);
	}

	SEXP BGL_ith_wavefront(SEXP num_verts_in, SEXP num_edges_in,
                         SEXP R_edges_in,   SEXP init_ind)
	{
		using namespace boost;

		Graph_ud g(num_verts_in, num_edges_in, R_edges_in);

		SEXP ansList, rbw;
		PROTECT(ansList = Rf_allocVector(VECSXP,1));
		PROTECT(rbw = Rf_allocVector(INTSXP, 1));

		INTEGER(rbw)[0] = ith_wavefront(vertex((int)INTEGER(init_ind)[0], g), g);

		SET_VECTOR_ELT(ansList,0,rbw);
		UNPROTECT(2);
		return(ansList);
	}

	SEXP BGL_max_wavefront(SEXP num_verts_in, SEXP num_edges_in,
                         SEXP R_edges_in)
	{
		using namespace boost;

		Graph_ud g(num_verts_in, num_edges_in, R_edges_in);

		SEXP ansList, rbw;
		PROTECT(ansList = Rf_allocVector(VECSXP,1));
		PROTECT(rbw = Rf_allocVector(INTSXP, 1));

		INTEGER(rbw)[0] = max_wavefront(g);

		SET_VECTOR_ELT(ansList,0,rbw);
		UNPROTECT(2);
		return(ansList);
	}

	SEXP BGL_aver_wavefront(SEXP num_verts_in, SEXP num_edges_in,
                         SEXP R_edges_in)
	{
		using namespace boost;

		Graph_ud g(num_verts_in, num_edges_in, R_edges_in);

		SEXP ansList, rbw;
		PROTECT(ansList = Rf_allocVector(VECSXP,1));
		PROTECT(rbw = Rf_allocVector(REALSXP, 1));

		REAL(rbw)[0] = aver_wavefront(g);

		SET_VECTOR_ELT(ansList,0,rbw);
		UNPROTECT(2);
		return(ansList);
	}

	SEXP BGL_rms_wavefront(SEXP num_verts_in, SEXP num_edges_in,
                         SEXP R_edges_in)
	{
		using namespace boost;

		Graph_ud g(num_verts_in, num_edges_in, R_edges_in);

		SEXP ansList, rbw;
		PROTECT(ansList = Rf_allocVector(VECSXP,1));
		PROTECT(rbw = Rf_allocVector(REALSXP, 1));

		REAL(rbw)[0] = rms_wavefront(g);

		SET_VECTOR_ELT(ansList,0,rbw);
		UNPROTECT(2);
		return(ansList);
	}
}

