#include "toric_ancilla_block2.h"
#include "flags.h"
#include "random.h"
#include <cmath>
#include <cassert>
#include <algorithm>
#include <cstdlib>


inline double H(const double &p) {
	return log( (1-p) / p );
}

inline double B(const double &p1, const double &p2) {
	return p1 * (1-p2) + p2 * (1-p1);
}

inline int I(int n, int i, int j) {
	return n * i + j;
}

inline int L(int n, int i, int j) {
	return j > 0? I(n, i, j-1): I(n, i, n-1);
}

inline int R(int n, int i, int j) {
	return j < n-1? I(n, i, j+1): I(n, i, 0);
}

inline int U(int n, int i, int j) {
	return i > 0? I(n, i-1, j): I(n, n-1, j);
}

ToricAncBlk :: ToricAncBlk(int m, int n, int T, double p) {
	this->m = m;
	this->n = n;
	this->T = T;
	this->p = p;
	extractor = NULL;
	edges[0] = edges[1] = NULL;
	decoder[0] = decoder[1] = NULL;
}

void ToricAncBlk :: build_circuit() {
	num_hori = (m+1) * n;
	num_vert = m * (n+1);
	num_xstabs = (m+1) * (n+1);
	num_zstabs = m * n;

	extractor = new Extractor(num_hori + num_vert, num_xstabs, num_zstabs);

	for (int i = 0; i < num_xstabs; i++) {
		extractor->set_syndrome_bit_x(i);
	}

	for (int i = 0; i < num_zstabs; i++) {
		extractor->set_syndrome_bit_z(i);
	}

	//North
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			extractor->add_CNOT(0, I(n, i, j), 2, I(n, i, j));
			extractor->add_2qubit_depo(0, I(n, i, j), 2, I(n, i, j), p);
		}
		for (int j = 0; j <= n; j++) {
			extractor->add_CNOT(1, I(n+1, i+1, j), 0, I(n+1, i, j)+num_hori);
			extractor->add_2qubit_depo(1, I(n+1, i+1, j), 0, I(n+1, i, j)+num_hori, p);
		}
	}

	//West
	for (int j = 0; j < n; j++) {
		for (int i = 0; i < m; i++) {
			extractor->add_CNOT(0, I(n+1, i, j)+num_hori, 2, I(n, i, j));
			extractor->add_2qubit_depo(0, I(n+1, i, j)+num_hori, 2, I(n, i, j), p);
		}
		for (int i = 0; i <= m; i++) {
			extractor->add_CNOT(1, I(n+1, i, j+1), 0, I(n, i, j));
			extractor->add_2qubit_depo(1, I(n+1, i, j+1), 0, I(n, i, j), p);
		}
	}
	
	//East
	for (int j = 0; j < n; j++) {
		for (int i = 0; i < m; i++) {
			extractor->add_CNOT(0, I(n+1, i, j+1)+num_hori, 2, I(n, i, j));
			extractor->add_2qubit_depo(0, I(n+1, i, j+1)+num_hori, 2, I(n, i, j), p);
		}
		for (int i = 0; i <= m; i++) {
			extractor->add_CNOT(1, I(n+1, i, j), 0, I(n, i, j));
			extractor->add_2qubit_depo(1, I(n+1, i, j), 0, I(n, i, j), p);
		}
	}


	//South
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			extractor->add_CNOT(0, I(n, i+1, j), 2, I(n, i, j));
			extractor->add_2qubit_depo(0, I(n, i+1, j), 2, I(n, i, j), p);
		}
		for (int j = 0; j <= n; j++) {
			extractor->add_CNOT(1, I(n+1, i, j), 0, I(n+1, i, j)+num_hori);
			extractor->add_2qubit_depo(1, I(n+1, i, j), 0, I(n+1, i, j)+num_hori, p);
		}
	}

	//Measure
	for (int i = 0; i < num_xstabs; i++) {
		extractor->add_phase_flip(1, i, p);
	}

	for (int i = 0; i < num_zstabs; i++) {
		extractor->add_bit_flip(2, i, p);
	}

}

void ToricAncBlk :: build_decoder_graph_one_round() {
	std::vector<int> E0[2], E1[2];
	double p_err;

	extractor->init_enumerator();
	edges[0] = new EdgeTable[num_xstabs];
	edges[1] = new EdgeTable[num_zstabs];

	while (extractor->enumerate(E0[0], E1[0], E0[1], E1[1], p_err)) {
		int u, v;
		std::pair<char, int> w;

		for (int t = 0; t < 2; t++) {
			if (E0[t].size() == 1 && E1[t].size() == 1) {  //time-like excitation
				u = E0[t][0],  v = E1[t][0],  w = std::make_pair(1, v);
			} else if (E0[t].size() == 2) {  //space-like excitation
				u = E0[t][0],  v = E0[t][1],  w = std::make_pair(0, v);
			} else if (E1[t].size() == 2) {
				u = E1[t][0],  v = E1[t][1],  w = std::make_pair(0, v);
			} else if (E0[t].size() + E1[t].size() == 1) {
				u = E0[t].size() == 1? E0[t][0]: E1[t][0];
				v = -1;
				w = std::make_pair(0, v);
			} else  continue;
			edges[t][u][w] = B(edges[t][u][w], p_err);
		}
	}
}

void ToricAncBlk :: build_decoder_graph() {
	build_decoder_graph_one_round();
	
	//X Decoder, no space boundary

	std::vector <bool> boundary_x( (T+1) * (m+1) * (n+1), false );

	//set time boundary
	for (int i = 0; i < num_xstabs; i++) {
		boundary_x[ I(num_xstabs, T, i) ] = true;
	}

	decoder[0] = new Decoder( boundary_x );

	for (int t = 0; t < T; t++) {
		for (int u = 0; u < num_xstabs; u++) {
			for (auto &p: edges[0][u]) {
				int tp = p.first.first;
				int v  = p.first.second;
				double w = H(p.second);
				decoder[0]->add_edge( I(num_xstabs, t, u), I(num_xstabs, t+tp, v), w );
			}
		}
	}

	//Z Decoder, with space boundary
	int num_node_layer = num_zstabs + m+m+n+n;
	std::vector <bool> boundary_z( T * num_node_layer + num_zstabs, false );

	//space boundary
	for (int t = 0; t < T; t++) {
		for (int i = 0; i < m+m+n+n; i++) {
			boundary_z[I(num_node_layer, t, num_zstabs+i)] = true;
		}
	}

	//time boundary
	for (int i = 0; i < num_zstabs; i++) {
		boundary_z[I(num_node_layer, T, i)] = true;
	}

	decoder[1] = new Decoder( boundary_z );

	for (int t = 0; t < T; t++) {
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				int u = I(n, i, j);
				for (auto &p: edges[1][u]) {
					int tp = p.first.first;
					int v  = p.first.second;
					double w = H(p.second);
					if (v == -1) {
						if (i == 0) {  //north boundary
							decoder[1]->add_edge(I(num_node_layer, t, u), I(num_node_layer, t, m*n+j), w);
						} else if (j == 0) {  //west 
							decoder[1]->add_edge(I(num_node_layer, t, u), I(num_node_layer, t, m*n+n+i), w);
						} else if (j == n-1) {  //east
							decoder[1]->add_edge(I(num_node_layer, t, u), I(num_node_layer, t, m*n+n+m+i), w);
						} else if (i == m-1) {  //south
							decoder[1]->add_edge(I(num_node_layer, t, u), I(num_node_layer, t, m*n+n+m+m+j), w);
						}
					} else {
						decoder[1]->add_edge(I(num_node_layer, t, u), I(num_node_layer, t+tp, v), w);
					}
				}
			}
		}
	}
}

bool ToricAncBlk :: simulate() {
	std::vector<int> data[2]; 
	std::vector<int> anc[2];
	std::vector<int> anc_last_round[2];

	decoder[0]->reset();
	decoder[1]->reset();
	
	data[0].resize(num_hori + num_vert, 0);
	data[1].resize(num_hori + num_vert, 0);
	anc[0].resize(num_xstabs, 0);
	anc_last_round[0].resize(num_xstabs, 0);
	anc[1].resize(num_zstabs, 0);
	anc_last_round[1].resize(num_zstabs, 0);

	int num_node_layer[2] = {num_xstabs, num_zstabs + m+m+n+n};

	for (int i = 0; i < data[1].size(); i++)  {
		if (Random::get_instance()->random_number() < 0.5)
			data[1][i] = 1;
	}

	for (int t = 0; t < T; t++) {
		extractor->execute(data[0], data[1], anc[0], anc[1], true);
		for (int o = 0; o < 2; o++) {
			for (int i = 0; i < anc[o].size(); i++) {
				if (anc_last_round[o][i] ^ anc[o][i]) {
					decoder[o]->excite(I(num_node_layer[o], t, i));
				}
			}
			anc_last_round[o] = anc[o];
		}
	}

	//correct Z errors
	std::vector< std::pair<int, int> > correct = decoder[0]->decode();
	for (auto &p: correct) {
		int u = p.first % num_node_layer[0];
		int v = p.second % num_node_layer[0];
		int x1 = u / (n+1);
		int y1 = u % (n+1);
		int x2 = v / (n+1);
		int y2 = v % (n+1);

		if (y1 != y2) {
			data[1][I(n, x1, std::min(y1, y2))] ^= 1;
		}
		if (x1 != x2) {
			data[1][I(n+1, std::min(x1, x2), y2) + num_hori] ^= 1;
		}
	}

	//correct X errors
	correct = decoder[1]->decode();
	for (auto &p: correct) {
		int u = p.first % num_node_layer[1];
		int v = p.second % num_node_layer[1];
		int x1, y1, x2, y2;
		if (u < num_zstabs) {
			x1 = u / n,  y1 = u % n;
		} else if (u < num_zstabs + n) { //north boundary
			x1 = -1,  y1 = u - num_zstabs;
		} else if (u < num_zstabs + n+m) {  //west boundary
			x1 = u - num_zstabs - n,  y1 = -1;
		} else if (u < num_zstabs + n+m+m) {  //east boundary
			x1 = u - num_zstabs - n - m,  y1 = n;
		} else {
			x1 = m,  y1 = u - num_zstabs - n - m - m;
		}


		if (v < num_zstabs) {
			x2 = v / n,  y2 = v % n;
		} else if (v < num_zstabs + n) { //north boundary
			x2 = -1,  y2 = v - num_zstabs;
		} else if (v < num_zstabs + n+m) {  //west boundary
			x2 = v - num_zstabs - n,  y2 = -1;
		} else if (v < num_zstabs + n+m+m) {  //east boundary
			x2 = v - num_zstabs - n - m,  y2 = n;
		} else {
			x2 = m,  y2 = v - num_zstabs - n - m - m;
		}

		if (y1 != y2) {
			data[0][I(n+1, x1, std::min(y1, y2)) + num_hori] ^= 1;
		}
		if (x1 != x2) {
			data[0][I(n, std::max(x1, x2), y2)] ^= 1;
		}
	}

	//detect errors
	extractor->execute(data[0], data[1], anc[0], anc[1], true);  //false for perfect detector
	for (int o = 0; o < 2; o++) {
		for (int i = 0; i < anc[o].size(); i++) {
			if (anc[o][i])  return false;
		}
	}

	/*
	for (int i = 0; i < n*3; i++)  putchar('-');
	puts("");
	for (int i = 0; i < m; i++) {
		//horizontal
		for (int j = 0; j < n; j++) {
			putchar('+');
			if (data[0][I(n, i, j)] && data[1][I(n,i,j)])  putchar('Y');
			else if (data[0][I(n,i,j)])  putchar('X');
			else if (data[1][I(n,i,j)])  putchar('Z');
			else  putchar('.');
		}
		puts("+");
		//vertical
		for (int j = 0; j <= n; j++) {
			int index = I(n+1,i,j) + num_hori;
			if (data[0][index] && data[1][index])  putchar('Y');
			else if (data[0][index])  putchar('X');
			else if (data[1][index])  putchar('Z');
			else  putchar('.');
			putchar(' ');
		}
		puts("");
	}

	for (int j = 0; j < n; j++) {
		putchar('+');
		int index = I(n, n, j);
		if (data[0][index] && data[1][index])  putchar('Y');
		else if (data[0][index])  putchar('X');
		else if (data[1][index])  putchar('Z');
		else  putchar('.');
	}
	puts("");
	*/
	return true;
}


