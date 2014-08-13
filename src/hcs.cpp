#include "RBGL.hpp"
#include "mincut.hpp"

// low-degree-vertice removal: external input of a decreasing sequence of int
std::vector<int> LDV;

int singleton_adoption_threshold=0; 

using namespace boost;
using namespace std;

typedef graph_traits < Graph_ud >::vertex_descriptor vertex_descriptor;
typedef graph_traits < Graph_ud >::vertex_iterator   vertex_iterator;
typedef graph_traits < Graph_ud >::edge_descriptor   edge_descriptor;
typedef graph_traits < Graph_ud >::edge_iterator     edge_iterator;
typedef graph_traits < Graph_ud >::edges_size_type   dst;

typedef std::vector<vertex_descriptor> Cluster;
typedef std::set<vertex_descriptor>    ClusterAsSet;
typedef std::vector<int>               V_Label;
typedef std::set<int>                  V_LabelAsSet;
typedef std::vector<V_Label>           ClusterResult;

#if DEBUG
static void output_graph_labels(const V_Label& v_label, const char* prompt)
{
    cout << prompt;
    for ( unsigned int j = 0; j < v_label.size(); j++ )
        cout << v_label[j] << " ";
    cout << endl << flush;
}

static void output_graph(const Graph_ud& g, const char* prompt)
{
    cout << prompt << " vertices in " << ": ";
    vertex_iterator vi, ve;
    for ( tie(vi, ve) = vertices(g); vi != ve; vi++ )
        cout << *vi+1 << " ";
    cout << endl;

    cout << prompt << " edges in " << ": ";
    edge_iterator ei, ee;
    vertex_descriptor u, v;
    for ( tie(ei, ee) = edges(g); ei != ee; ei++ )
    {
        u = source(*ei, g), v = target(*ei, g);
        cout << "(" << u+1 << "," << v+1 << ")" << " ";
    }
    cout << endl << flush;
}
#endif

static void build_subgraph(Cluster& vs_vec,
                           const Graph_ud& g_in, const V_Label& v_label_in,
                           Graph_ud& g_out, V_Label& v_label_out)
{
    // turning "vector" to "set" is just to ease lookup later
    Cluster::iterator vi;
    ClusterAsSet vs_set;
    for ( vi = vs_vec.begin(); vi != vs_vec.end(); vi++ )
        vs_set.insert(*vi);

    // keep track of vertex labels
    for ( unsigned int ii = 0; ii < v_label_in.size(); ii++ )
        if ( vs_set.count(ii) == 1 )
            v_label_out.push_back(v_label_in[ii]);

    // add edges of the subgraph
    vertex_descriptor u, v;
    edge_descriptor e1;
    bool inserted;
    edge_iterator ei, ei_end;
    for (tie(ei, ei_end) = edges(g_in); ei != ei_end; ++ei)
    {
        u = source(*ei, g_in), v = target(*ei, g_in);
        if ( vs_set.count(u) == 1 && vs_set.count(v) == 1 )
            tie(e1, inserted) = add_edge(u, v, g_out);
    }

    // remove vertices NOT in vs_set
    // complication comes from renumbering of vertices after removal
    vertex_iterator vvi, vvi_end;
    bool try_again = TRUE;

    V_Label flags(num_vertices(g_out), 0);
    V_Label::iterator fi;
    for ( tie(vvi, vvi_end) = vertices(g_out); vvi != vvi_end; vvi++ )
    {
        if ( vs_set.count(*vvi) == 1 ) flags[*vvi] = 1;
    }

    while ( try_again )
    {
        try_again = FALSE;
        fi = flags.begin();
        for ( tie(vvi, vvi_end) = vertices(g_out); vvi != vvi_end; vvi++, fi++ )
        {
            if ( *fi == 0 )
            {
                boost::remove_vertex(*vvi, g_out);
		flags.erase(fi);
                try_again = TRUE;
                break;
            }
        }
    }

#if DEBUG
    output_graph(g_out, "build_subgraph: ");
    output_graph_labels(v_label_out, "build_subgraph: ");
#endif
}

// repeatedly remove all vertices of degree < LDV[i] from G
// NOTE: "vertices w/ degree < LDV[i]" is implemented as dynamic
static void remove_vertices(const int d,
                            const Graph_ud& g_in, const V_Label& v_label_in,
                            Graph_ud& g_out, V_Label& v_label_out)
{
    g_out = g_in;
    v_label_out = v_label_in;

    vertex_iterator vvi, vvi_end;
    bool try_again = TRUE;
    V_Label::iterator i;

    while ( try_again )
    {
        try_again = FALSE;
        i = v_label_out.begin();
        for ( tie(vvi, vvi_end) = vertices(g_out); vvi != vvi_end; vvi++, i++ )
        {
            if ( out_degree(*vvi, g_out) < (unsigned int) d )
            {
                boost::clear_vertex(*vvi, g_out);
                boost::remove_vertex(*vvi, g_out);
                v_label_out.erase(i);
                try_again = TRUE;
                break;
            }
        }
    }
}

static dst cluster_distance(const Graph_ud& g, V_Label& v_label, int si, V_Label& ci)
{
    dst d = 0;
    unsigned int s_index;

    for ( s_index = 0; s_index < v_label.size(); s_index++ )
        if ( v_label[s_index] == si )
            break;

    if ( s_index >= v_label.size() ) return 0;  // ERROR

    V_LabelAsSet ci_set(ci.begin(), ci.end());
    Graph_ud::adjacency_iterator ai, ai_end;
    for ( tie(ai, ai_end) = adjacent_vertices(s_index, g);
            ai != ai_end; ai++ )
        if ( ci_set.count(v_label[*ai]) == 1 ) d++;

#if DEBUG
    cout << " cluster_distance between "
         << si << " and " << ci[0] << " is " << d
         << endl << flush;
#endif

    return d;
}

static void adopt_singleton(const Graph_ud& g, V_Label& v_label, ClusterResult& clusters)
{
    V_LabelAsSet singletons;
    V_LabelAsSet::iterator si, simax;

    ClusterResult::iterator ci, cimax;

    for ( ci = clusters.begin(); ci != clusters.end(); ci++ )
        if ( (*ci).size() == 1 )
            singletons.insert((*ci)[0]);

#if DEBUG
    cout << " found singletons: ";
    for ( si = singletons.begin(); si != singletons.end(); si++ )
        cout << (*si) << " ";
    cout << endl << flush;
#endif

    dst d, dmax;
    bool try_again = TRUE;

    while ( try_again )
    {
        try_again = FALSE;

        for ( si = singletons.begin(); si != singletons.end() && !try_again; )
        {
            dmax = 0;
            cimax = clusters.begin();
            for ( ci = clusters.begin(); ci < clusters.end(); ci++ )
            {
                d = cluster_distance(g, v_label, *si, *ci);
                if ( d >= dmax && (*ci).size() > 1 )
                {
                    cimax = ci;
                    dmax = d;
                }
            }

#if DEBUG
            cout << " max distance from " << *si
            << " to cluster " << (*cimax)[0] << " = " << dmax
            << endl << flush;
            if ( dmax > singleton_adoption_threshold && (*cimax).size() > 1 )
            {
                cout << " inserting " << *si
                << " to cluster " << (*cimax)[0]
                << endl << flush;
            }
#endif

            if ( (int) dmax > singleton_adoption_threshold && (*cimax).size() > (unsigned int) 1 )
            {
                // copy si to cur then increment si
                V_LabelAsSet::iterator cur = si++;
                (*cimax).push_back(*cur);
                singletons.erase(cur);
                try_again = TRUE;
            }
            else
            {
                ++si;
            }
        }
    }

    // remove singletons from clusters
    try_again = TRUE;
    while ( try_again )
    {
        try_again = FALSE;
        for ( ci = clusters.begin(); ci != clusters.end() && !try_again; ci++ )
        {
            if ( (*ci).size() == 1 && singletons.count((*ci)[0]) == 0 )
            {
                clusters.erase(ci);
                try_again = TRUE;
            }
        }
    }
}

static void remove_clusters(ClusterResult& clusters,
                            const Graph_ud& g_in,  const V_Label& v_label_in,
                            Graph_ud& g_out, V_Label& v_label_out)
{
    ClusterResult::iterator ci;
    V_Label::iterator       vi;
    V_LabelAsSet            v_set;

    for ( ci = clusters.begin(); ci < clusters.end(); ci++ )
    {
        for ( vi = (*ci).begin(); vi != (*ci).end(); vi++ )
            v_set.insert(*vi);
    }

    Cluster v_vec;

    // TODO: go over vertices instead of label indices??
    for ( unsigned int i = 0; i < v_label_in.size(); i++ )
        if ( v_set.count(v_label_in[i]) == 0 )
            v_vec.push_back(i);

    Graph_ud gt_out(num_vertices(g_in));
    V_Label vt_label_out;

    build_subgraph(v_vec, g_in, v_label_in, gt_out, vt_label_out);
 
    g_out = gt_out;
    v_label_out = vt_label_out;
}

// basic algorithm:
// HCS( G(V, E) )
// {
//   (H, HH, cut) <- min-cut(G)
//   if ( |cut| >= |V| / 2 )
//      return V
//   else
//   {
//      HCS(H)
//      HCS(HH)
//   }
// }
static void HCS_internal(const Graph_ud& g, V_Label& v_label, ClusterResult& clusters)
{
    Cluster s_vec, vs_vec;
    dst cut_capacity = ( num_vertices(g) <= 1 ) ? 0 :
                       min_cut(g, std::back_inserter(s_vec),
                               std::back_inserter(vs_vec) );

    if ( v_label.size() <= 0 )
    {
#if DEBUG
        cout << " 0 vertex in the graph: NEED TO LOOK AT THIS " << endl << flush;
#endif
    }
    else if ( v_label.size() <= 1 )
    {
        clusters.push_back(v_label);

#if DEBUG
        cout << " 1 vertex in the graph " << v_label[0] << endl << flush;
#endif
    }
    else if ( ((unsigned int) cut_capacity) >= num_vertices(g) / 2 )
    {
        clusters.push_back(v_label);

#if DEBUG
        cout << "found a highly connected subgraph with "
        << num_vertices(g) << " vertices: ";

        V_Label::iterator vi;
        for ( vi = v_label.begin(); vi != v_label.end(); ++vi)
            cout << *vi << " ";

        cout << endl << flush;
#endif

    }
    else
    {
        V_Label v_label_out1;
        V_Label v_label_out2;

        Graph_ud g1(num_vertices(g)), g2(num_vertices(g));

        build_subgraph(s_vec,  g, v_label, g1, v_label_out1);
        build_subgraph(vs_vec, g, v_label, g2, v_label_out2);

        ClusterResult cluster1, cluster2;
        HCS_internal(g1, v_label_out1, cluster1);
        HCS_internal(g2, v_label_out2, cluster2);

        ClusterResult::iterator ci;
        for ( ci = cluster1.begin(); ci != cluster1.end(); ci++ )
            clusters.push_back(*ci);
        for ( ci = cluster2.begin(); ci != cluster2.end(); ci++ )
            clusters.push_back(*ci);
    }
}

// refinements:
// (1) iterated HCS: handle cases w/ several mincuts in a graph
// (2) singleton adoption adoption
// (3) low degree heuristic: eliminate low-degree vertices before HCS
//
// final algorithm: (d[1..p]: input int array in decreasing order)
// HCS_LOOP( G(V, E) )
// {
//   for ( i = 1, i <= p, i++ )
//   {
//     H <- G;
//     repeatdly remove all vertices of degree < d_i from H
//
//     repeat
//     {
//        HCS(H)
//        perform singletons adoption
//        remove clustered vertices from H
//     } until ( no new cluster is found )
//
//     remove clusterd vertices from G
//   }
// }

static SEXP HCS_loop(Graph_ud& g, V_Label& v_label)
{
    unsigned int i, j, nn;

    Graph_ud g1, g2;
    V_Label v_label_orig1, v_label_orig2;

    ClusterResult clusterG, clusterH, clusterHCS;
    ClusterResult::iterator ci;

    V_Label::iterator vi;

    for ( i = 0; i < LDV.size(); i++ )
    {
        if ( i )
        {
            remove_clusters(clusterH, g, v_label, g1, v_label_orig1);
            g = g1;
            v_label = v_label_orig1;
        }

        // COMMENT OUT to test "adopt_singleton"
        remove_vertices(LDV[i], g, v_label, g1, v_label_orig1);
        
        if ( num_vertices(g1) == 0 ) 
        {
           g1 = g;
	   v_label_orig1 = v_label;
        }

#if DEBUG
        output_graph(g1, "after removing vertices: ");
        output_graph_labels(v_label_orig1, "after removing vertices: ");
#endif

        clusterH.clear();

        while ( num_vertices(g1) )
        {
            clusterHCS.clear();
            HCS_internal(g1, v_label_orig1, clusterHCS);

            adopt_singleton(g1, v_label_orig1, clusterHCS);

            remove_clusters(clusterHCS, g1, v_label_orig1, g2, v_label_orig2);
            g1 = g2;
            v_label_orig1 = v_label_orig2;

            for ( j = 0; j < clusterHCS.size(); j++ )
                clusterH.push_back(clusterHCS[j]);
        }

        for ( j = 0; j < clusterH.size(); j++ )
            clusterG.push_back(clusterH[j]);
    }

    // convert from ClusterResult to SEXP list
    SEXP ansList, cList;
    PROTECT(ansList = allocVector(VECSXP,(int)clusterG.size()));

    for ( i = 0, ci = clusterG.begin(); ci != clusterG.end(); i++, ci++ )
    {
        // convert one cluster to one SEXP vector
        nn = (*ci).size();
        PROTECT(cList = allocVector(INTSXP, nn));
        for ( j = 0, vi = (*ci).begin(); vi != (*ci).end(); j++, vi++ )
            INTEGER(cList)[j] = *vi;

        SET_VECTOR_ELT(ansList,i,cList);
        UNPROTECT(1);
    }

    UNPROTECT(1);
    return(ansList);
}

extern "C"
{

    // highly-connected-subgraphs: clustering based on graph connectivity
    SEXP BGL_highly_conn_sg (SEXP num_verts_in, SEXP num_edges_in,
                             SEXP R_edges_in,   SEXP R_weights_in,
                             SEXP sat, SEXP lldv, SEXP ldv)
    {
        int i, ne;
        Graph_ud g(num_verts_in, num_edges_in, R_edges_in, R_weights_in);

        // checking of parameters are done in R codes
        singleton_adoption_threshold = INTEGER(sat)[0];
        ne = INTEGER(lldv)[0];

        LDV.clear();
        for ( i = 0; i < ne; i++ ) LDV.push_back(INTEGER(ldv)[i]);

        // assign labels to vertices
        // "(i+1)*2" is to distinguish label from index (when testing)
        // "(i+1)" makes the vertex labels the same as those in R graph
        ne = num_vertices(g);
        V_Label v_label(ne);
        for ( i = 0; i < ne; i++ ) v_label[i] = (i+1);

        SEXP ansList = HCS_loop(g, v_label);

        return(ansList);
    }
}

