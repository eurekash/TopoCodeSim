//
// Created by eurek on 5/29/2019.
//

#ifndef COMPASS_V2_DECODER_H
#define COMPASS_V2_DECODER_H


#include <tuple>
#include <vector>
#include <set>
#include <tuple>
#include <queue>
#include <list>

using std::vector;
using std::pair;
using std::tuple;

struct list_node {
    int         value;
    list_node*  left;
    list_node*  right;
    list_node(int _value = -1): value(_value), left(NULL), right(NULL)  {}
};

class Decoder {
public:
    Decoder (const std::vector <bool>& _is_boundary);
    void    add_edge(int u, int v, double w);
    void    excite(int u);      
    void    reset();
    vector < pair<int, int> > decode();
	double	total_time;
	int     counter;

private:

    void    syndrome_validation();


    //information of each syndrome check (vertex)
    bool    odd(int u);
    int 	num_vertices;
    vector <bool>  	is_boundary;
    vector <bool>   defect;
    vector <int>	ndefects;
    vector <bool>	open;

    //data structure for the decoder graph
    int 	num_edges;
    vector <int>  	start;   	//for each edge e = (u, v), start[e] := u,  end[e] := v
    vector <int>  	end;		
    vector <double> 	weight;
    vector <double> 	cover;
    vector <int>  	oppo;		//opposite edge

    //disjoint sets for syndrome validation
    int     root_synd_valid(int u);
    vector <int>  	prt_synd_valid;
    vector <int>  	size_synd_valid;

    //disjoint sets for peeling
    int     root_peel(int u);
    vector <int>  	prt_peel;
    vector <int>  	size_peel;

    //clusters
    void            print_clusters();

    void            insert_edge (int e, int c);
    void 			remove_edge (int e);
    bool            empty(int c);
    int             merge_cluster (int c1, int c2);

	list_node*		pool;
    vector <list_node*>     edge_nodes;
    vector <list_node*>     cluster_heads;
    vector <int>            cluster_size;

    void            print_active_list();
    void            activate(int c);
    void            deactivate(int c);
    vector <list_node*>     vertex_nodes;
    vector <list_node*>     active_heads;
    vector <bool>   active;

    //peeling
    void    peel();
    bool    DFS(int, bool&);

    vector < tuple<double, int, int, int> > edge_list;
    vector < vector<int> >    MST_edges;
    vector < pair<int, int> > result;
    vector <bool> visited;

};

#endif //COMPASS_V2_DECODER_H
