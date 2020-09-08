#include "toric_ancilla_block.h"
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
	extractor_x = NULL;
	extractor_z = NULL;
	edges = NULL;
	decoder = NULL;
}

void ToricAncBlk :: build_circuit() {
	num_hori = (m+1) * n;
	num_vert = m * (n+1);
	num_xstabs = (m+1) * (n+1);
	num_zstabs = m * n;

	extractor_x = new Extractor(num_hori + num_vert, num_xstabs, 0);
	extractor_z = new Extractor(num_hori + num_vert, 0, num_zstabs);

	for (int i = 0; i < num_xstabs; i++) {
		extractor_x->set_syndrome_bit_x(i);
	}

	for (int i = 0; i < num_zstabs; i++) {
		extractor_z->set_syndrome_bit_z(i);
	}

	//North
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			extractor_z->add_CNOT(0, I(n, i, j), 2, I(n, i, j));
			extractor_z->add_2qubit_depo(0, I(n, i, j), 2, I(n, i, j), p);
		}
		for (int j = 0; j <= n; j++) {
			extractor_x->add_CNOT(1, I(n+1, i+1, j), 0, I(n+1, i, j)+num_hori);
			extractor_x->add_2qubit_depo(1, I(n+1, i+1, j), 0, I(n+1, i, j)+num_hori, p);
		}
	}

	//West
	for (int j = 0; j < n; j++) {
		for (int i = 0; i < m; i++) {
			extractor_z->add_CNOT(0, I(n+1, i, j)+num_hori, 2, I(n, i, j));
			extractor_z->add_2qubit_depo(0, I(n+1, i, j)+num_hori, 2, I(n, i, j), p);
		}
		for (int i = 0; i <= m; i++) {
			extractor_x->add_CNOT(1, I(n+1, i, j+1), 0, I(n, i, j));
			extractor_x->add_2qubit_depo(1, I(n+1, i, j+1), 0, I(n, i, j), p);
		}
	}
	
	//East
	for (int j = 0; j < n; j++) {
		for (int i = 0; i < m; i++) {
			extractor_z->add_CNOT(0, I(n+1, i, j+1)+num_hori, 2, I(n, i, j));
			extractor_z->add_2qubit_depo(0, I(n+1, i, j+1)+num_hori, 2, I(n, i, j), p);
		}
		for (int i = 0; i <= m; i++) {
			extractor_x->add_CNOT(1, I(n+1, i, j), 0, I(n, i, j));
			extractor_x->add_2qubit_depo(1, I(n+1, i, j), 0, I(n, i, j), p);
		}
	}


	//South
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			extractor_z->add_CNOT(0, I(n, i+1, j), 2, I(n, i, j));
			extractor_z->add_2qubit_depo(0, I(n, i+1, j), 2, I(n, i, j), p);
		}
		for (int j = 0; j <= n; j++) {
			extractor_x->add_CNOT(1, I(n+1, i, j), 0, I(n+1, i, j)+num_hori);
			extractor_x->add_2qubit_depo(1, I(n+1, i, j), 0, I(n+1, i, j)+num_hori, p);
		}
	}

	//Measure
	for (int i = 0; i < num_xstabs; i++) {
		extractor_x->add_phase_flip(1, i, p);
	}

	for (int i = 0; i < num_zstabs; i++) {
		extractor_z->add_bit_flip(2, i, p);
	}

}

void ToricAncBlk :: build_decoder_graph_one_round() {
	std::vector<int> E0, E1, E0p, E1p;
	double p_err;

	extractor_x->init_enumerator();
	edges = new EdgeTable[num_xstabs];

	while (extractor_x->enumerate(E0, E1, E0p, E1p, p_err)) {
		int u, v;
		std::pair<char, int> w;

		if (E0.size() == 1 && E1.size() == 1) {  //time-like excitation
			u = E0[0],  v = E1[0],  w = std::make_pair(1, v);
		} else if (E0.size() == 2) {  //space-like excitation
			u = E0[0],  v = E0[1],  w = std::make_pair(0, v);
		} else if (E1.size() == 2) {
			u = E1[0],  v = E1[1],  w = std::make_pair(0, v);
		} else if (E0.size() + E1.size() == 1) {
			u = E0.size() == 1? E0[0]: E1[0];
			v = -1;
			w = std::make_pair(0, v);
		} else  continue;
		edges[u][w] = B(edges[u][w], p_err);
	}
}

void ToricAncBlk :: build_decoder_graph() {
	build_decoder_graph_one_round();
	
	decoder = new Decoder( std::vector<bool> ( T*(m+1)*(n+1), false ) );

	for (int t = 0; t < T; t++) {
		for (int u = 0; u < num_xstabs; u++) {
			for (auto &p: edges[u]) {
				int tp = p.first.first;
				int v  = p.first.second;
				double w = H(p.second);
				if (t == 0 && tp == 0) {
					decoder->add_edge( I(num_xstabs, t, u), I(num_xstabs, t + tp, v), 0 );
				} else if (t + tp < T) {
					decoder->add_edge( I(num_xstabs, t, u), I(num_xstabs, t + tp, v), w );
				}
			}
		}
	}
}

bool ToricAncBlk :: simulate() {
	std::vector<int> data[2]; 
	std::vector<int> anc[2];
	std::vector<int> anc_last_round[2];

	decoder->reset();
	
	data[0].resize(num_hori + num_vert, 0);
	data[1].resize(num_hori + num_vert, 0);
	anc[0].resize(num_xstabs, 0);
	anc_last_round[0].resize(num_xstabs, 0);
	anc[1].resize(num_zstabs, 0);
	anc_last_round[1].resize(num_zstabs, 0);

	for (int i = 0; i < data[1].size(); i++)  {
		if (Random::get_instance()->random_number() < 0.5)
			data[1][i] = 1;
	}

	int nexcitation = 0;
	for (int t = 0; t < T; t++) {
		extractor_x->execute(data[0], data[1], anc[0], anc[1], true);
		for (int i = 0; i < anc[0].size(); i++) {
			if (anc_last_round[0][i] ^ anc[0][i]) {
				decoder->excite(I(num_xstabs, t, i));
				nexcitation++;
			}
		}
		anc_last_round[0] = anc[0];
	}

	if (nexcitation % 2 == 1)  return false;
	//correct Z errors
	std::vector< std::pair<int, int> > correct = decoder->decode();
	for (auto &p: correct) {
		int u = p.first % num_xstabs;
		int v = p.second % num_xstabs;
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

	//detect errors

	extractor_z->execute(data[0], data[1], anc[0], anc[1], true);  //false for perfect detector
	for (int i = 0; i < anc[1].size(); i++) {
		if (anc[1][i])  return false;
	}

	extractor_x->execute(data[0], data[1], anc[0], anc[1], true);  //false for perfect detector
	for (int i = 0; i < anc[0].size(); i++) {
		if (anc[0][i])  return false;
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


