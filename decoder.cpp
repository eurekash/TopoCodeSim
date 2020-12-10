// Created by eurek on 5/29/2019.
//

#include <algorithm>
#include <assert.h>
#include <time.h>
#include "decoder.h"

#define COUNT
//#define __DEBUG__

using std :: vector;
using std :: pair;
using std :: tuple;

#define INFINITY (1e9)


Decoder :: Decoder (const std::vector <bool>& _is_boundary) :
is_boundary ( _is_boundary ),
num_vertices ( _is_boundary.size() ),  num_edges(0),
defect ( num_vertices ),  ndefects ( num_vertices ),  open ( num_vertices ),
start(0), end(0), weight(0), cover(0), oppo(0),
prt_synd_valid( num_vertices ), size_synd_valid( num_vertices ),
prt_peel( num_vertices ),  size_peel ( num_vertices ),
cluster_size( num_vertices ),
edge_nodes(0), active_heads(0), active( num_vertices ),
edge_list(0),  MST_edges( num_vertices ),  visited( num_vertices ),
cluster_heads( num_vertices, NULL ),
vertex_nodes( num_vertices, NULL ),
dist( num_vertices, -1 ),
adjList( num_vertices, vector<int>() ),
total_time(0)
{
    counter = 0;

	pool = new list_node[num_vertices*2];
	list_node *pt = pool;

	for (int i = 0; i < num_vertices; i++, pt++) {
		cluster_heads[i] = pt;
		*pt = list_node();
	}

    for (int i = 0; i < num_vertices; i++, pt++) {
        vertex_nodes[i] = pt;
	    *pt	= list_node(i);
    }
}

bool Decoder :: odd(int u) {
    return ndefects[u] % 2 == 1 && open[u] == false;
}

void Decoder :: reset() {
    std::sort(edge_list.begin(), edge_list.end(), std::greater<std::tuple<double, int, int, int> > ());
    //reset lists of active clusters

    for (int i = 0; i < num_edges; i++) {
        list_node *h = active_heads[i];
        h -> left = h -> right = h;
    }

    for (int u = 0; u < num_vertices; u ++) {
        open[u] = is_boundary[u];
        defect[u] = false;
        ndefects[u] = 0;
        prt_synd_valid[u] = u; 
        size_synd_valid[u] = 1;
        prt_peel[u] = u;
        size_peel[u] = 1;
        active[u] = false;
        cluster_size[u] = 0;
        list_node *h = cluster_heads[u];
        h -> left = h -> right = h;
    }

    std :: fill(cover.begin(), cover.end(), 0.0);
    for (int e = 0; e < num_edges; e++) {
        insert_edge(e, start[e]);
    }

    result.resize(0);

}

void Decoder :: excite(int u) {
    defect[u] = true;
    ndefects[u] = 1;
    activate(u);
}

void Decoder :: add_edge(int u, int v, double w) {

    int e_uv = num_edges++;
    int e_vu = num_edges++;

    edge_list.emplace_back( std::make_tuple(w, e_uv, u, v) );

    start.emplace_back(u);
    end.emplace_back(v);
    weight.emplace_back(w);
    cover.emplace_back(0.);
    oppo.emplace_back(e_vu);
    edge_nodes.emplace_back(new list_node(e_uv));
    active_heads.emplace_back(new list_node());
	adjList[u].emplace_back(v);


    start.emplace_back(v);
    end.emplace_back(u);
    weight.emplace_back(w);
    cover.emplace_back(0.);
    oppo.emplace_back(e_uv);
    edge_nodes.emplace_back(new list_node(e_vu));
    active_heads.emplace_back(new list_node());
	adjList[v].emplace_back(u);
}

void Decoder :: insert_edge(int e, int c) {
    #ifdef __DEBUG__
    //printf("call Decoder::insert_edge(%d,%d)\n", e, c);
    #endif
    cluster_size[c] ++;
    list_node *h = cluster_heads[c];
    list_node *p = edge_nodes[e];
    p -> right = h;
    p -> left  = h -> left;
    h -> left -> right = p;
    h -> left  = p;
}

void Decoder :: remove_edge(int e) {
    int c = root_synd_valid(start[e]);
    counter += 4;
    cluster_size[c] --;

    list_node *p = edge_nodes[e];
    p -> left -> right = p -> right;
    p -> right -> left = p -> left;
}

void Decoder :: deactivate(int c) {

    counter += 5;
    if (active[c] == false)
        return ;

    active[c] = false;
    list_node *p = vertex_nodes[c];
    p -> left -> right = p -> right;
    p -> right -> left = p -> left ;
}

void Decoder :: activate(int c) {
    counter += 10;

    if (open[c] || (ndefects[c] & 1) == 0)    //should not be activated
        return ;
    if (active[c] == true)
        return ;

    active[c] = true;
    int s = cluster_size[c];
    list_node *p = vertex_nodes[c];
    list_node *h = active_heads[s];
    p -> left  = h;
    p -> right = h -> right;
    p -> right -> left = p;
    h -> right = p;
}

int Decoder :: merge_cluster(int c1, int c2) {   //merge c2 to c1
    c1 = root_synd_valid(c1);
    c2 = root_synd_valid(c2);

    if (c1 == c2)   
        return c1;

    counter += 20;

    if (size_synd_valid[c1] < size_synd_valid[c2])  
        std :: swap(c1, c2);

    //update the disjoint set structure
    prt_synd_valid[c2]  =  c1;
    size_synd_valid[c1] += size_synd_valid[c2];
    ndefects[c1] += ndefects[c2];
    open[c1] = open[c1] or open[c2];

    if (cluster_size[c2] == 0)    //c2 is empty
        return c1;  

    list_node *h1 = cluster_heads[c1];
    list_node *h2 = cluster_heads[c2];
    list_node *lh1 = h1 -> left;
    list_node *lh2 = h2 -> left;
    list_node *rh2 = h2 -> right;

    //concatenation
    h1 -> left      =   lh2;
    lh2 -> right    =   h1;
    lh1 -> right    =   rh2;
    rh2 -> left     =   lh1;

    //update the size of c1
    cluster_size[c1] += cluster_size[c2];

    //clear c2
    cluster_size[c2] = 0;
    h2 -> left = h2 -> right = h2;

    return c1;
}

int Decoder :: root_synd_valid(int u) {
    int &prt = prt_synd_valid[u];
    counter ++;
    return u == prt? u: prt = root_synd_valid(prt);
}

int Decoder :: root_peel(int u) {
    int &prt = prt_peel[u];
    return u == prt? u: prt = root_peel(prt);
}

vector< pair<int,int> > Decoder :: decode() {
	double start_syndrome_validation = clock();
    syndrome_validation();
    peel();
    return result;
}

void Decoder :: syndrome_validation() {

    #ifdef __DEBUG__
    printf("called Decoder::syndrome_validation()\n");
    #endif

    vector <int>  merge_list (num_vertices);
    vector <bool>  touched (num_vertices, false);

    for (int s = 1; s < num_edges; s++) {
        list_node *ahead = active_heads[s];
        while (ahead -> right != ahead) {   
            int c = ahead -> right -> value;
            //print_active_list();
            //print_clusters();
            deactivate(c);
            while (cluster_size[c] <= s) {
                //printf("current cluster = %d\n", c);
                double increment = INFINITY;   
                //find the increment value by enumerating all edges in c
                list_node *head = cluster_heads[c];
                for (list_node *p = head -> right; p != head; p = p -> right) {
                	counter ++;
                    int edgeId = p -> value;
                    if (root_synd_valid (end[edgeId]) == c) {
                        remove_edge(edgeId);
                    } else {
                        double delta = weight[edgeId] - cover[edgeId] - cover[oppo[edgeId]];
                        if (increment > delta)   increment = delta;
                    }
                }

                int total = 0;
                for (list_node *p = head -> right; p != head; p = p -> right) {
                	counter ++;
                    int edgeId = p -> value;
                    cover[edgeId] += increment;

                    if (cover[edgeId] + cover[ oppo[edgeId] ] >= weight[edgeId] - 1e-9) {
                        int c2 = root_synd_valid(end[edgeId]);

                        if (touched[c2] == false) {
                            merge_list[total ++] = root_synd_valid(end[edgeId]);
                            touched[c2] = true;
                        }

                        remove_edge(edgeId);
                        remove_edge(oppo[edgeId]);
                    }
                }

                //merge all clusters in the merge list
                while (total --) {
                    int &c2 = merge_list[total];
                    touched[c2] = false;
                    deactivate(c2);
                    c = merge_cluster(c, c2);
                }

                if (odd(c) == false)  break;
            }
            activate(c);
        }
    }
    //printf("counter = %d\n", counter);
}

void Decoder :: peel() {
    #ifdef __DEBUG__
    printf("call Decoder::peel()\n");
    #endif

    for (int i = 0; i < num_vertices; i++) {
        prt_peel[i] = i;
        size_peel[i] = 1;
        MST_edges[i].resize(0);
        visited[i] = false;
    }

    for (const auto &e: edge_list) {
        int id = std::get<1>(e);
        int u = std::get<2>(e);
        int v = std::get<3>(e);

        if (cover[id] + cover[oppo[id]] < weight[id] - 1e-9)  continue;
        if (root_synd_valid(u) == root_synd_valid(v) && root_peel(u) != root_peel(v)) {
            int root_u = root_peel(u);
            int root_v = root_peel(v);

            if (size_peel[root_u] >= size_peel[root_v]) {
                size_peel[root_u] += size_peel[root_v];
                prt_peel[root_v] = root_u;
            } else {
                size_peel[root_v] += size_peel[root_u];
                prt_peel[root_u] = root_v;
            }

            MST_edges[u].emplace_back(v);
            MST_edges[v].emplace_back(u);
        }
    }

    for (int i = 0; i < num_vertices ; i++) {
        int r = root_synd_valid(i);
        if (!visited[r]) {
            bool fuse_with_boundary = (ndefects[r] % 2 == 1);
            DFS(r, fuse_with_boundary);
        }
    }
}

bool Decoder::DFS(int u, bool& fuse_with_boundary) {
    visited[u] = true;

    bool flag = defect[u];
    if (is_boundary[u] && fuse_with_boundary) {
        flag = true;
        fuse_with_boundary = false;
    }
    for (auto &v: MST_edges[u]) {
        if (!visited[v] && DFS(v, fuse_with_boundary)) {
            result.emplace_back(std::make_pair(u, v));
            flag = !flag;
        }
    }
    return flag;
}

void Decoder :: init_shortest_path(int n2) {
	for (int i = 0; i < n2; i++) {
		dist[i] = 0;
		Q.push(i);
	}
}

bool Decoder :: shortest_path(int newNode, int dist_threshold) {
	while (!Q.empty() && Q.front() < newNode) {
		int u = Q.front();
		Q.pop();
		for (auto &v: adjList[u]) {
			if (dist[v] == -1) {
				dist[v] = dist[u] + 1;
				Q.push(v);
			}
		}
	}
	return Q.empty() || dist[Q.front()] >= dist_threshold;
}
